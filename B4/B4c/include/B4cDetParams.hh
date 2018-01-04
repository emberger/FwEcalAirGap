#ifndef B4cDetParams
#define B4cDetParams

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "globals.hh"



//Specifyes the geometrical information

class DetParams {

public:

static DetParams& Instance(){

								static DetParams Det;
								return Det;

}



~DetParams();

void InitDet();

void SetfNofOuterLayers(G4double nla);
G4double GetfNofOuterLayers();

void SetOutertileLenX(G4double tx);
G4double GetOutertileLenX();

void SetOutertileLenY(G4double ty);
G4double GetOutertileLenY();

void SetOuterabsoThickness(G4double abs);
G4double GetOuterabsoThickness();

void SetOutergapThickness(G4double gap);
G4double GetOutergapThickness();



void SetfNofInnerLayers(G4double nla);
G4double GetfNofInnerLayers();

void SetInnertileLenX(G4double tx);
G4double GetInnertileLenX();

void SetInnertileLenY(G4double ty);
G4double GetInnertileLenY();

void SetInnerabsoThickness(G4double abs);
G4double GetInnerabsoThickness();

void SetInnergapThickness(G4double gap);
G4double GetInnergapThickness();

void SetPvesselThickness(G4double pvt);
G4double GetPvesselThickness();


void SetcalorSizeXY(G4double cs);
G4double GetcalorSizeXY();

//void SetairgapThickness(G4double at);
//G4double GetairgapThickness();

void SetWorldMult(G4double wm);

G4double GetnofInnerTilesX();
G4double GetnofInnerTilesY();

G4double GetnofOuterTilesX();
G4double GetnofOuterTilesY();

G4double GetDetectorThickness();

G4double GetInnerlayerThickness();
G4double GetOuterlayerThickness();

G4double GetInnercalorThickness();
G4double GetOutercalorThickness();

G4double GetInnertilesPerLayer();
G4double GetOutertilesPerLayer();

G4double GetWorldSizeXY();
G4double GetWorldSizeZ();

void SetInnerAbsFirst(G4bool Icfg);
G4bool GetInnerAbsFirst();

void SetOuterAbsFirst(G4bool Ocfg);
G4bool GetOuterAbsFirst();

//G4double GetGunPos();

private:

G4double fNofInnerLayers;
G4double fNofOuterLayers;

G4double OutertileLenX;
G4double OutertileLenY;
G4double nofOuterTilesX;
G4double nofOuterTilesY;
G4double OutertilesPerLayer;

G4double OuterabsoThickness;
G4double OutergapThickness;
G4double OuterlayerThickness;

G4double InnertileLenX;
G4double InnertileLenY;
G4double nofInnerTilesX;
G4double nofInnerTilesY;
G4double InnertilesPerLayer;

G4double InnerabsoThickness;
G4double InnergapThickness;
G4double InnerlayerThickness;



G4double calorSizeXY;
//G4double airgapThickness;
G4double DetectorThickness;

G4double PvesselThickness;

G4double InnercalorThickness;
G4double OutercalorThickness;

G4double WorldSizeXY;
G4double WorldSizeZ;
G4double WorldMult;

G4double GunPos;

G4bool Innerabsofirst;
G4bool Outerabsofirst;
DetParams(){
}
DetParams(const DetParams& other){
}

};


inline DetParams& GetInst(){
								return DetParams::Instance();
}

#endif
