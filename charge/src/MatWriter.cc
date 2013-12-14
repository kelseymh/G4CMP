#include "MatWriter.hh" 
#include "inttypes.h"
#include <fstream>
#include <cstdio>

using std::vector;

MatWriter::MatWriter(G4String outFile, G4String varName, G4int numRows, G4int numColumns, vector<G4double> data, G4bool append)
{

    int32_t matType = 14; //Matlab "matrix" flag.
                         
    // Now I will typecast the input variables into size-definite types.
    int32_t rows = int32_t(numRows);
    int32_t columns = int32_t(numColumns);
    int32_t nameLength = int32_t(varName.length());
    int32_t dataSize = int32_t(columns*rows*8); //8 bytes per double
    // These are Matlab tag ID's
    int32_t sizeOfMatrix = dataSize + 8*7; //seven tags, eight bytes per tag.
    int8_t mxDOUBLE_CLASS = 6; //double precision
    int32_t miUINT32 = 6;
    int32_t miINT32 = 5;
    int32_t miDOUBLE = 9;
    int32_t miINT8 = 1;
    // And here are some numbers
    int8_t zero = 0;
    int32_t eight = 8;

    FILE* fid;
    if(!append)
    {
        fid = fopen(outFile.c_str(),"w");

        int16_t version = 0x0100; //Matlab needs this to know it's a mat file
        char header[42] = {"Matlab MAT file written by CDMS Geant DMC"};
        int8_t space = 32; //ASCII space
        int16_t MI = 0x4D49; //The letters "MI" in hex. 
                             //If Matlab reads this and sees "IM," it knows 
                             //to byte-swap the file.
        // Write text header (124 bytes)
        fwrite(header,1,42,fid);
        for(int i=0; i<124-42; ++i)
            fwrite(&space,1,1,fid);

        //Version and Endianness checker
        fwrite(&version,2,1,fid);
        fwrite(&MI,2,1,fid);
    }
    else
    {
        fid = fopen(outFile.c_str(),"a");
    }

    fwrite(&matType,4,1,fid);
    fwrite(&sizeOfMatrix,4,1,fid);

    fwrite(&miUINT32,4,1,fid);
    fwrite(&eight,4,1,fid);
    fwrite(&mxDOUBLE_CLASS,1,1,fid);
    for(int i=0; i<7; ++i)
        fwrite(&zero,1,1,fid);

    fwrite(&miINT32,4,1,fid);
    fwrite(&eight,4,1,fid);
    fwrite(&rows,4,1,fid);
    fwrite(&columns,4,1,fid);

    fwrite(&miINT8,4,1,fid);
    fwrite(&nameLength,4,1,fid);
    fwrite(varName,1,nameLength,fid);

    fwrite(&miDOUBLE,4,1,fid);
    fwrite(&dataSize,4,1,fid);
    for(int j=0; j<columns; ++j)
        for(int i=0; i<rows; ++i)
            fwrite(&data[i+(j*rows)],8,1,fid);

    fclose(fid);
}
