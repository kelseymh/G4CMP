#ifndef CDMS_FET1_h
#define CDMS_FET1_h
#include <vector>
#include "globals.hh"
#include "G4Run.hh"
#include "G4ThreeVector.hh"
#include "FETMessenger1.hh"

using std::vector;

class FET1
{
    public:
        FET1();
        ~FET1();
        void Run();

    private:
        //FET constants
        G4int channels, ramoModel, timeBins;
        G4double decayTime, deltaZ, dt, preTrig, templateEnergy;
        vector<G4double> radIn;
        vector<G4double> radOut;
        vector<G4double> Z;
        vector<G4int> quad;
        vector<G4String> channelNames;
        //Phonon constants
        G4bool Fano;
        G4int atomicNumber, massNumber;
        G4double elecHoleGap, chargeEnergy, fanoFactor, elecHoleCreation;
        //Calculated stuff
        G4double fetStart;
        //vector< vector< vector<G4double> > > FETTemplate;
        vector< vector<G4double> > fetTrace;
        vector<vector<G4double> > elecPos;
        vector<G4double> elecT;
        vector<vector<G4double> > holePos;
        vector<G4double> holeT;
        vector<G4double> charBias;
        FETMessenger1* messenger;

    private:
        void ReadFETConstantsFile(G4String filename);
        void ReadPhononConstantsFile(G4String filename);
        //G4double CalculateChargeSupressionFactor(const G4Run* run);
        void BuildFETTemplate(vector<vector<vector<G4double> > >& FETTemplate);
        void CalculateTrace();
        void ReadRamoInputFile(const G4String& filename, vector<vector<G4double> >& X, vector<G4double>& V);
        void SaveResults(const G4String& matFile, const G4String& varName, G4int numRows, G4int numColumns,const vector<G4double>& colData, G4bool append);
        G4bool QuadrantFlip(const G4int& quadrant, vector< vector<G4double> >& elecPos, vector< vector<G4double> >& holePos);
};
#endif
