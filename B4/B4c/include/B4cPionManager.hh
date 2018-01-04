#ifndef B4cPionManager_hh
#define B4cPionManager_hh


#include "globals.hh"
#include "G4ThreeVector.hh"


class PionManager {

public:


~PionManager();

static PionManager& Instance(){

        static PionManager EMan;
        return EMan;

}

void SetEnergyPh1(G4double eph1);
G4double GetEnergyPh1();

void SetEnergyPh2(G4double eph2);
G4double GetEnergyPh2();

void SetMomPh1(G4ThreeVector mp1);
G4ThreeVector GetMomPh1();

void SetMomPh2(G4ThreeVector mp2);
G4ThreeVector GetMomPh2();

void SetTD();
void ReSetTD();
bool GetTD();


private:

PionManager(){
}

G4double EnergyPh1;
G4double EnergyPh2;

G4ThreeVector MomPh1;
G4ThreeVector MomPh2;

bool td;

};


inline PionManager& GetEInst(){
        return PionManager::Instance();
}




#endif
