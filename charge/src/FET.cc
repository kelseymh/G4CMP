#include <vector>
#include <math.h>
#include <stdlib.h>
#include "FET.hh"
#include "inttypes.h"
#include <fstream>
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "AlminumElectrodeHit.hh"
#include "AlminumElectrodeSensitivity.hh"
#include <ctime>
#include <algorithm>
#include "Interpolation3D.hh"
#include <sstream>
#include <string>
#include "G4SDManager.hh"
#include "MatWriter.hh"
#include "FETMessenger.hh"

using std::string;
using std::vector;
using std::min;
using std::max;

FET::FET(G4RunManager* runMan):RunManager(runMan)
{
    messenger = new FETMessenger(this);
}

FET::~FET()
{}

void FET::Run()
{
    ReadFETConstantsFile(G4String("ConstantsFET.txt"));
    ReadPhononConstantsFile(G4String("ConstantsPhonon.txt"));
    CalculateTrace(RunManager->GetCurrentRun());
}

void FET::ReadFETConstantsFile(G4String filename)
{
    std::ifstream constants;
    constants.open(filename);
    G4String buffer;
    G4String varName;
    G4String varVal;
    G4bool comment = false;
    size_t pos;

    while(!constants.eof())
    {
        getline(constants, buffer);
        if((!comment) && (buffer.find("//") != string::npos))
        {
            if(buffer.find("//") == 0)
                continue;
            else
                buffer = buffer.substr(0, buffer.find("//") + 2);
        }
        else if((!comment) && (buffer.find("/*") != string::npos))
        {
            buffer.erase(buffer.find("/*"), buffer.find("*/") + 2);
            if(buffer.find("*/") == string::npos)
                comment = true;
        }
        else if(buffer.find("*/") != string::npos)
        {
            comment = false;
            buffer.erase(buffer.find("*/"), 2);
        }

        if(!comment && !buffer.empty())
        {
            while(buffer.find(' ') != string::npos) //Erase all spaces
                buffer.erase(buffer.find(' '), 1);

            pos = buffer.find('='); //To split buffer into lhs and rhs
            if(pos == string::npos)
                continue;
            else
            {
                varName = buffer.substr(0, pos++);
                varVal = buffer.substr(pos);

                if(*(--varVal.end()) == ';') //Strip ; from end
                    varVal.erase(--varVal.end());

                G4String num;
                if(varName == "channels")
                    channels = atoi(varVal);
                else if(varName == "ramoModel")
                    ramoModel = atoi(varVal);
                else if(varName == "timeBins")
                    timeBins = atoi(varVal);
                else if(varName == "decayTime")
                    decayTime = atof(varVal)*s;
                else if(varName == "deltaZ")
                    deltaZ = atof(varVal)*m;
                else if(varName == "dt")
                    dt = atof(varVal)*s;
                else if(varName == "preTrig")
                    preTrig = atof(varVal)*s;
                else if(varName == "templateEnergy")
                    templateEnergy = atof(varVal)*eV;
                else if(varName == "radIn")
                {
                    varVal.erase(varVal.begin()); //Strip { and } from line
                    varVal.erase(--varVal.end());
                    while(!varVal.empty())
                    {
                        radIn.push_back(atof(varVal.substr(0, varVal.find(',')).c_str())*m);
                        if(varVal.find(',') == string::npos)
                            varVal.clear();
                        else
                            varVal.erase(0, varVal.find(',') + 1);
                    }
                }
                else if(varName == "radOut")
                {
                    varVal.erase(varVal.begin()); //Strip { and } from line
                    varVal.erase(--varVal.end()); 
                    while(!varVal.empty())
                    {
                        radOut.push_back(atof(varVal.substr(0, varVal.find(',')).c_str())*m);
                        if(varVal.find(',') == string::npos)
                            varVal.clear();
                        else
                            varVal.erase(0, varVal.find(',') + 1);
                    }
                }
                else if(varName == "Z")
                {
                    varVal.erase(varVal.begin()); //Strip { and } from line
                    varVal.erase(--varVal.end());
                    while(!varVal.empty())
                    {
                        Z.push_back(atof(varVal.substr(0, varVal.find(',')).c_str())*m);
                        if(varVal.find(',') == string::npos)
                            varVal.clear();
                        else
                            varVal.erase(0, varVal.find(',') + 1);
                    }
                }
                else if(varName == "quad")
                {
                    varVal.erase(varVal.begin()); //Strip { and } from line
                    varVal.erase(--varVal.end());
                    while(!varVal.empty())
                    {
                        quad.push_back(atoi(varVal.substr(0, varVal.find(',')).c_str()));
                        if(varVal.find(',') == string::npos)
                            varVal.clear();
                        else
                            varVal.erase(0, varVal.find(',') + 1);
                    }
                }
                else if(varName == "channelNames")
                {
                    varVal.erase(varVal.begin()); //Strip { and } from line
                    varVal.erase(--varVal.end());
                    while(!varVal.empty())
                    {
                        if(varVal.find(',') == string::npos)
                        {
                            channelNames.push_back(varVal.substr(1, varVal.length()-2));
                            varVal.clear();
                        }
                        else
                        {
                            channelNames.push_back(varVal.substr(1, varVal.find(',')-2));
                            varVal.erase(0, varVal.find(',') + 1);
                        }
                    }
                }
                else
                    G4cout << "FETSim Warning: Variable " << varName 
                           << " is unused." << G4endl;
            }
        }
    }
}

void FET::ReadPhononConstantsFile(G4String filename)
{
    std::fstream constants(filename);
    G4String buffer;
    G4String varName;
    G4String varVal;
    G4bool comment = false;
    size_t pos;

    while(!constants.eof())
    {
        getline(constants, buffer);
        if((!comment) && (buffer.find("//") != string::npos))
        {
            if(buffer.find("//") == 0)
                continue;
            else
                buffer = buffer.substr(0, buffer.find("//") + 2);
        }
        else if((!comment) && (buffer.find("/*") != string::npos))
        {
            buffer.erase(buffer.find("/*"), buffer.find("*/") + 2);
            if(buffer.find("*/") == string::npos)
                comment = true;
        }
        else if(buffer.find("*/") != string::npos)
        {
            comment = false;
            buffer.erase(buffer.find("*/"), 2);
        }

        if(!comment && !buffer.empty())
        {
            while(buffer.find(' ') != string::npos) //Erase all spaces
                buffer.erase(buffer.find(' '), 1);

            pos = buffer.find('=');
            if(pos == string::npos)
                continue;
            else
            {
                varName = buffer.substr(0, pos++);
                varVal = buffer.substr(pos);

                if(varName != "elecHoleGap" && varName != "chargeEnergy"
                        && varName != "massNumber" && varName != "atomicNumber"
                        && varName != "Fano" && varName != "fanoFactor"
                        && varName != "elecHoleCreation")
                    continue;

                if(*(--varVal.end()) == ';') //Strip ; from end
                    varVal.erase(--varVal.end());

                if(varName == "elecHoleGap")
                    elecHoleGap = atof(varVal)*eV;
                else if(varName == "chargeEnergy")
                    chargeEnergy = atof(varVal)*eV;
                else if(varName == "atomicNumber")
                    atomicNumber = atoi(varVal);
                else if(varName == "massNumber")
                    massNumber = atoi(varVal);
                else if(varName == "Fano")
                    Fano = atoi(varVal);
                else if(varName == "fanoFactor")
                    fanoFactor = atof(varVal);
                else if(varName == "elecHoleCreation")
                    elecHoleCreation = atof(varVal)*eV;
                else
                    G4cout << "FETSim Warning: Variable " << varName 
                           << " is unused." << G4endl;;
            }
        }
    }
}

G4double FET::CalculateChargeSupressionFactor(const G4Run* run)
{
    const vector<const G4Event*>* eventVec = run->GetEventVector();
    G4int numEvents = eventVec->size();
    G4HCofThisEvent* pHCofEvent;
    AlminumElectrodeHitsCollection* pHitColl;
    G4double hitEnergy;
    vector<G4double> eventEnergy(numEvents,0);
    G4double totalEnergy = 0;
    G4SDManager* fSDM = G4SDManager::GetSDMpointer();
    G4int collectionID = fSDM->GetCollectionID("AlminumElectrodeHit");
    for(int i=0; i<numEvents; ++i)
    {
        pHCofEvent = (eventVec->at(i))->GetHCofThisEvent();
        pHitColl = static_cast<AlminumElectrodeHitsCollection*>(pHCofEvent->GetHC(collectionID));
        for(int k=0; k<pHitColl->entries(); ++k)
        {
            hitEnergy = (static_cast<AlminumElectrodeHit*>(pHitColl->GetVector()->at(k)))->GetEDep();
            eventEnergy[i] += hitEnergy;
            totalEnergy += hitEnergy;
        }
    }

    G4double simChargeEnergy = 0;
    G4double fullChargeEnergy = 0;
    G4double simChargeEnergyEvent;
    G4double fullChargeEnergyEvent;
    G4double chFraction;
    G4int ncc; //Number of charge carriers
    for(int i=0; i<numEvents; ++i) 
    {
        chFraction = eventEnergy[i]/totalEnergy; //Assume elec recoil
        ncc = 4*floor(chargeEnergy*chFraction/elecHoleGap/4 + .5); //use floor(x+.5) to round.
        simChargeEnergyEvent = ncc * elecHoleGap;
        simChargeEnergy += simChargeEnergyEvent;

        srand(time(NULL));
        if(Fano)
            fullChargeEnergyEvent = elecHoleGap*floor(eventEnergy[i]/elecHoleCreation +
                    rand()/RAND_MAX *sqrt(fanoFactor*eventEnergy[i]/elecHoleCreation)); //Assume elec recoil.
        else
            fullChargeEnergyEvent = elecHoleGap*floor(eventEnergy[i]/elecHoleCreation);

        fullChargeEnergy += fullChargeEnergyEvent;
    }

    return(simChargeEnergy/fullChargeEnergy);
}

void FET::BuildFETTemplate(vector<vector<vector<G4double> > >& FETTemplate)
{
    std::fstream templateFile("FETTemplates.txt");
    if(!templateFile.fail())
    {
        for(int i=0; i<4; ++i)
            for(int j=0; j<4; ++j)
                for(int k=0; k<4096; ++k)
                {
                    templateFile >> FETTemplate[i][j][k];
                }
    }
    else
    {
        for(int i=0; i<channels; ++i)
        {
            for(int j=0; j<preTrig/dt; ++j)
                FETTemplate[i][i][j] = 0;

            for(int j=1; j<timeBins-preTrig/dt+1; ++j)
                FETTemplate[i][i][j+preTrig/dt-1] = exp(-j*dt/decayTime);
        }
    }
    templateFile.close();
}

void FET::CalculateTrace(const G4Run* run)
{
    vector<G4int> potChan2FETChan(channels,0);

    std::fstream ramoFile("RamoChannels.txt");
    vector<G4String> ramoChannelNames;
    G4String buffer;
    while(!ramoFile.eof())
    {
        ramoFile >> buffer;
        ramoChannelNames.push_back(buffer);
    }
    ramoFile.close();

    for(int i=0; i<channels; ++i)
        for(int j=0; j<channels; ++j)
            if(ramoChannelNames[i] == channelNames[j])
            {
                potChan2FETChan[i] = j;
                break;
            }

    G4double chargeSupressionFactor = CalculateChargeSupressionFactor(run);

    fetTrace = vector< vector<G4double> >(channels, vector<G4double>(timeBins,0));
    charBias.resize(channels);

    const vector<const G4Event*>* eventVec = run->GetEventVector();
    G4int numEvents = eventVec->size();

    G4HCofThisEvent* pHCofEvent;
    AlminumElectrodeHitsCollection* pHitColl;
    AlminumElectrodeHit hit;
    G4SDManager* fSDM = G4SDManager::GetSDMpointer();
    G4int collectionID = fSDM->GetCollectionID("AlminumElectrodeHit");

    G4ThreeVector positionVec;
    vector<G4double> position(3,0);
    for(int q=0; q<quad.size(); ++q)
    {
        for(int i=0; i<numEvents; ++i)
        {
            pHCofEvent = (eventVec->at(i))->GetHCofThisEvent();
            pHitColl = static_cast <AlminumElectrodeHitsCollection*>(pHCofEvent->GetHC(collectionID));
            G4cout << "pHitColl->GetSize() = " << pHitColl->entries() << G4endl;
            for(int j=0; j<pHitColl->entries(); ++j)
            {
                hit = *(static_cast<AlminumElectrodeHit*>(pHitColl->GetVector()->at(j)));
                if(hit.GetCharge() < 0)
                {
                    positionVec = hit.GetWorldPos();
                    for(int k=0; k<3; ++k)
                        position[k] = positionVec[k];
                    elecPos.push_back(position);
                    elecT.push_back(hit.GetTime());
                }
                else
                {
                    positionVec = hit.GetWorldPos();
                    for(int k=0; k<3; ++k)
                        position[k] = positionVec[k];
                    holePos.push_back(position);
                    holeT.push_back(hit.GetTime());
                }
            }
        }
        QuadrantFlip(quad[q],elecPos,holePos);

        vector<vector<vector<G4double> > > FETTemplate(channels, 
                vector<vector<G4double> >(channels, vector<G4double>(timeBins,0)));
        BuildFETTemplate(FETTemplate);

        vector<Interpolation3D> ramo;
        vector<vector<G4double> > X;
        vector<vector<G4double> > V;
        vector<G4double> tempV;

        ReadRamoInputFile("EpotRamoChan1.txt",X,tempV);
        V.push_back(tempV);
        ramo.push_back(Interpolation3D(X,tempV));

        ReadRamoInputFile("EpotRamoChan2.txt",X,tempV);
        V.push_back(tempV);
        ramo.push_back(Interpolation3D(X,tempV));

        ReadRamoInputFile("EpotRamoChan3.txt",X,tempV);
        V.push_back(tempV);
        ramo.push_back(Interpolation3D(X,tempV));

        ReadRamoInputFile("EpotRamoChan4.txt",X,tempV);
        V.push_back(tempV);
        ramo.push_back(Interpolation3D(X,tempV));

        vector<G4double> elecPot(elecPos.size(),0);
        vector<G4double> holePot(holePos.size(),0);
        for(int channelRamo = 0; channelRamo < channels; ++channelRamo)
        {
            G4double minV = 0;
            G4double maxV = 0;
            for(int i=0; i<V[channelRamo].size(); ++i)
            {
                if(V[channelRamo][i] > maxV)
                    maxV = V[channelRamo][i];
                else if(V[channelRamo][i] < minV)
                    minV = V[channelRamo][i];
            }
            if(maxV > 0)
                charBias[potChan2FETChan[channelRamo]] = maxV;
            else
                charBias[potChan2FETChan[channelRamo]] = minV;

            if(ramoModel == 1)
            {
                for(int i=0; i<elecPos.size(); ++i)
                {
                    if( (abs(sqrt(elecPos[i][0]*elecPos[i][0]+elecPos[i][1]*elecPos[i][1])) 
                                >= radIn[channelRamo]) 
                        && (abs(sqrt(elecPos[i][0]*elecPos[i][0]+elecPos[i][1]*elecPos[i][1])
                                <= radOut[channelRamo]) 
                        && elecPos[i][2] >= Z[channelRamo]-abs(deltaZ)
                        && elecPos[i][2] <= Z[channelRamo]+abs(deltaZ)))
                        elecPot[i] = 1;
                    else
                        elecPot[i] = 0;
                }
            }
            else if(ramoModel == 2)
            {
                for(int i=0; i<elecPos.size(); ++i)
                {
                    elecPot[i] = ramo[potChan2FETChan[channelRamo]].GetValue(&elecPos[i][0]);
                    holePot[i] = ramo[potChan2FETChan[channelRamo]].GetValue(&holePos[i][0]);
                }
            }
            else
                G4cout << "This is bad- no Ramo Model specified" << G4endl;

            fetStart=DBL_MAX;
            for(int i=0; i<holeT.size(); ++i)
                if(min(holeT[i],elecT[i])<fetStart)
                    fetStart = min(holeT[i],elecT[i]);

            for(int chargeNum=0; chargeNum<holeT.size(); ++chargeNum)
                for(int channelDAQ=0; channelDAQ<channels; ++channelDAQ)
                    for(int traceIdx=0; traceIdx<timeBins; ++traceIdx)
                    {
                        fetTrace[channelDAQ][traceIdx] += 
                            elecPot[chargeNum]*FETTemplate[channelDAQ][channelRamo][traceIdx]
                            /chargeSupressionFactor;
                        fetTrace[channelDAQ][traceIdx] += 
                            -1*holePot[chargeNum]*FETTemplate[channelDAQ][channelRamo][traceIdx]
                            /chargeSupressionFactor;
                    }
        }
        std::stringstream out;
        out << quad[q];
        vector<G4double> fetOutput(4096*4,0);
        for(int i=0; i<4096; ++i)
            for(int j=0; j<4; ++j)
                fetOutput[j+(i*4)] = fetTrace[j][i]/volt;
        G4String name = "FETTrace";
        SaveResults("FETSimOutput_Quad" + out.str() + ".mat", name, 4, 4096, fetOutput, 0);
        name = "fetStart";
        vector<G4double> fetstartoutput(1,fetStart/s);
        SaveResults("FETSimOutput_Quad" + out.str() + ".mat", name, 1, 1, fetstartoutput, 1);
        name = "charBias";
        for(int i=0; i<4; ++i)
            charBias[i] = charBias[i]/volt;
        SaveResults("FETSimOutput_Quad" + out.str() + ".mat", name, 4, 1, charBias, 1);
    }
}

void FET::ReadRamoInputFile(const G4String& filename, vector<vector<G4double> >& X, vector<G4double>& V)
{
    X.clear();
    V.clear();
    // Read in File that has x,y,z,V
    std::fstream input(filename.c_str());
	    
    G4double x,y,z,v;
    G4int lineNum=0;
    G4int n=0;
    vector<G4double> tempVec(3,0);
    while(!input.eof())
    {
        input >> x >> y >> z >> v;

        if(X.size()>=100000*n)
        {
            ++n;
            X.reserve(100000*n);
            V.reserve(100000*n);
        }
            
        tempVec[0] = x*m;
        tempVec[1] = y*m;
        tempVec[2] = z*m - 1.27*cm;
        X.push_back(tempVec);
        V.push_back(v*volt);

        ++lineNum;
    }
    X.reserve(lineNum);
    V.reserve(lineNum);

    input.close();
}

void FET::SaveResults(const G4String& matFile, const G4String& varName, G4int numRows, G4int numColumns,const vector<G4double>& colData, G4bool append)
{
    MatWriter(matFile, varName, numRows, numColumns, colData, append);
}

G4bool FET::QuadrantFlip(const G4int &quadrant, vector< vector<G4double> > &elecPos, vector< vector<G4double> > &holePos)
{
    if (quadrant == 2 || quadrant == 3)
    {
        for(int i = 0; i < elecPos.size(); ++i)
            elecPos[i][1] = -1*elecPos[i][1];
        for(int i = 0; i < holePos.size(); ++i)
            holePos[i][1] = -1*holePos[i][1];
    }
    if (quadrant == 3 || quadrant == 4)
    {
        for(int i = 0; i < elecPos.size(); ++i)
            elecPos[i][2] = -1*elecPos[i][2];
        for(int i = 0; i < holePos.size(); ++i)
            holePos[i][2] = -1*holePos[i][2];
    }
    return(true);
}
