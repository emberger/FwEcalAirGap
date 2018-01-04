#include "B4cPionManager.hh"


PionManager::~PionManager(){
}

void PionManager::SetEnergyPh1(G4double eph1){
        EnergyPh1=eph1;
}
G4double PionManager::GetEnergyPh1(){
        return EnergyPh1;
}

void PionManager::SetEnergyPh2(G4double eph2){
        EnergyPh2=eph2;
}
G4double PionManager::GetEnergyPh2(){
        return EnergyPh2;
}

void PionManager::SetMomPh1(G4ThreeVector mp1){
        MomPh1=mp1;
}
G4ThreeVector PionManager::GetMomPh1(){
        return MomPh1;
}

void PionManager::SetMomPh2(G4ThreeVector mp2){
        MomPh2=mp2;
}
G4ThreeVector PionManager::GetMomPh2(){
        return MomPh2;
}

void PionManager::SetTD(){
        td=true;
}

void PionManager::ReSetTD(){
        td=false;
}


bool PionManager::GetTD(){
        return td;
}
