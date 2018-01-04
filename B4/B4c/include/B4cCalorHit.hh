//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B4cCalorHit.hh 100946 2016-11-03 11:28:08Z gcosmo $
//
/// \file B4cCalorHit.hh
/// \brief Definition of the B4cCalorHit class

#ifndef B4cCalorHit_h
#define B4cCalorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4Threading.hh"

/// Calorimeter hit class
///
/// It defines data members to store the the energy deposit and track lengths
/// of charged particles in a selected volume:
/// - fEdep, fTrackLength

class B4cCalorHit : public G4VHit
{
public:
B4cCalorHit();
B4cCalorHit(const B4cCalorHit&);
virtual ~B4cCalorHit();

// operators
const B4cCalorHit& operator=(const B4cCalorHit&);
G4int operator==(const B4cCalorHit&) const;

inline void* operator new(size_t);
inline void operator delete(void*);

// methods from base class
virtual void Draw() {
}
virtual void Print();

// methods to handle data
void Add(G4double de, G4double dl);

void SetTouch();
void SetCellInfo();
// void SetPhotonNumber(G4int nr);
void SetX(G4int x);
void SetY(G4int y);
void SetZ(G4int z);
void SetCalorPart(std::string p);



// get methods
G4double GetEdep() const;
G4double GetTrackLength() const;

std::string GetCalorPart() const;

G4bool GetTouch();
G4bool GetCellInfo();


//G4int GetPhotonNumber();

G4int GetX();
G4int GetY();
G4int GetZ();
//G4int GetCalorSeg();


private:
G4double fEdep;            ///< Energy deposit in the sensitive volume
G4double fTrackLength;    ///< Track length in the  sensitive volume

G4bool fTouched;      // Bool to select interesting cells
G4bool fCellInfo;     // Bool to distinguish cell info from layered info and total accounting
//G4int Photon;

std::string calpart;

G4int Xpos;     //
G4int Ypos;     // Coordinates to identify Calorimeter Cell
G4int Zpos;     //

//G4int CalSeg; //calorimeter Segment

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using B4cCalorHitsCollection = G4THitsCollection<B4cCalorHit>;

extern G4ThreadLocal G4Allocator<B4cCalorHit>* B4cCalorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* B4cCalorHit::operator new(size_t)
{
        if (!B4cCalorHitAllocator) {
                B4cCalorHitAllocator = new G4Allocator<B4cCalorHit>;
        }
        void *hit;
        hit = (void *) B4cCalorHitAllocator->MallocSingle();
        return hit;
}

inline void B4cCalorHit::operator delete(void *hit)
{
        if (!B4cCalorHitAllocator) {
                B4cCalorHitAllocator = new G4Allocator<B4cCalorHit>;
        }
        B4cCalorHitAllocator->FreeSingle((B4cCalorHit*) hit);
}

inline void B4cCalorHit::Add(G4double de, G4double dl) {
        fEdep += de;
        fTrackLength += dl;
}

inline G4double B4cCalorHit::GetEdep() const {
        return fEdep;
}

inline std::string B4cCalorHit::GetCalorPart() const{
  return calpart;
}

inline G4double B4cCalorHit::GetTrackLength() const {
        return fTrackLength;
}

inline G4bool B4cCalorHit::GetTouch(){
        return fTouched;
}

inline void B4cCalorHit::SetTouch(){
        fTouched=true;
}

inline void B4cCalorHit::SetCalorPart(std::string p){
  calpart=p;
}

inline void B4cCalorHit::SetCellInfo(){
        fCellInfo=true;
}

inline G4bool B4cCalorHit::GetCellInfo(){
        return fCellInfo;
}

inline void B4cCalorHit::SetX(G4int x){
        Xpos=x;
}
inline void B4cCalorHit::SetY(G4int y){
        Ypos=y;
}
inline void B4cCalorHit::SetZ(G4int z){
        Zpos=z;
}

inline G4int B4cCalorHit::GetX(){
        return Xpos;
}
inline G4int B4cCalorHit::GetY(){
        return Ypos;
}
inline G4int B4cCalorHit::GetZ(){
        return Zpos;
}

// inline void B4cCalorHit::SetCalorSeg(G4int seg){
//         CalSeg=seg;
// }
//
// inline G4int B4cCalorHit::GetCalorSeg(){
//         return CalSeg;
// }


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
