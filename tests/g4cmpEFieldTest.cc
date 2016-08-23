/***********************************************************************\
 * This software is licensed under the terms of the GNU General Public *
 * License version 3 or later. See G4CMP/LICENSE for the full license. *
\***********************************************************************/
/*
 * A test for checking changes to the EField code.
 * Takes in the x, y, and z dimensions of a rectangular box, the
 * number of points, and an EPot file name. Uniformly distributes
 * points in the rectangular box centered at the origin using the
 * given potential from the EPot file.
 * Outputs a file with the voltage and E field components at each position.
 * Also prints how long it takes to run.
 */
#include "G4CMPMeshElectricField.hh"
#include <ctime>
#include<fstream>

int main(int argc, char** argv) {

    G4double lx, ly, lz;
    G4int N;
    G4String filename;
    time_t start, end;
    if (argc == 6) {
        lx = atof(argv[1]);
        ly = atof(argv[2]);
        lz = atof(argv[3]);
        N = atoi(argv[4]);
        filename = argv[5];
    }
    if (argc != 6) {
        G4cout << "g4cmpEFieldTest x_length y_length z_length number_points EPotfilename" << G4endl;
        return 0;
    }

    std::fstream outputFile;
    outputFile.open("EFieldTestdata.txt", std::fstream::out);
    if (!outputFile.is_open()) {
        std::cerr << "Cannot create file, error was: " << strerror(errno) << std::endl;
    }
    outputFile << "Data formatted as x, y, z, V, Ex, Ey, Ez" << "\r\n";

    std::time(&start);

    G4CMPMeshElectricField EField(filename);
    G4double EMfield[6];
    //G4double normE;
    G4int n = cbrt(N);
    G4double Deltax = lx/n;
    G4double Deltay = ly/n;
    G4double Deltaz = lz/n;

    for (G4int i = 0; i < n; ++i)
        for (G4int j = 0; j < n; ++j)
            for (G4int k = 0; k < n; ++k) {
                G4double x = -(Deltax * n/2) + (i * Deltax);
                G4double y = -(Deltay * n/2) + (j * Deltay);
                G4double z = -(Deltaz * n/2) + (k * Deltaz);
                outputFile << x << " " << y << " " << z << " ";
                const G4double point[4] = {x, y, z, 0};
                outputFile << EField.GetPotential(point) << " ";
                EField.GetFieldValue(point, EMfield);
                outputFile << EMfield[3] << " " << EMfield[4] << " " << EMfield[5] << " ";
                //normE = EMfield[3]*EMfield[3] + EMfield[4]*EMfield[4] + EMfield[5]*EMfield[5];
                //outputFile << normE;
                outputFile << "\r\n";
            }

    outputFile.close();

    std::time(&end);
    G4cout << "This took " << difftime(end, start) << " s" << G4endl;

    return 0;
}
