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
// $Id: B4PrimaryGeneratorAction.cc 100946 2016-11-03 11:28:08Z gcosmo $
//
/// \file B4PrimaryGeneratorAction.cc
/// \brief Implementation of the B4PrimaryGeneratorAction class

#include "B4PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include <time.h>
#include "TStyle.h"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4PrimaryGeneratorAction::B4PrimaryGeneratorAction()
        : G4VUserPrimaryGeneratorAction(),
        fParticleGun(nullptr)
{
        G4int nofParticles = 1;
        fParticleGun = new G4ParticleGun(nofParticles);

        // default particle kinematic
        //
        auto particleDefinition
                = G4ParticleTable::GetParticleTable()->FindParticle("pi0");
        fParticleGun->SetParticleDefinition(particleDefinition);



        G4double worldZHalfLength = 0.;
        auto worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");

        // Check that the world volume has box shape
        G4Box* worldBox = nullptr;
        if (  worldLV ) {
                worldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
        }

        if ( worldBox ) {
                worldZHalfLength = worldBox->GetZHalfLength();
        }
        else  {
                G4ExceptionDescription msg;
                msg << "World volume of box shape not found." << G4endl;
                msg << "Perhaps you have changed geometry." << G4endl;
                msg << "The gun will be place in the center.";
                G4Exception("B4PrimaryGeneratorAction::GeneratePrimaries()",
                            "MyCode0002", JustWarning, msg);
        }

        TCanvas * dummy = new TCanvas();
        f = new TFile("/home/iwsatlas1/emberger/Geant4/Current/4piSensitiveDetector_mod/B4-build/pi0_kinematics.root");
        h1 = (TH2D*)f->Get("pi0_energyAngle");
        std::cout<<"file opened"<<std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4PrimaryGeneratorAction::~B4PrimaryGeneratorAction()
{
        delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
        // // This function is called at the begining of event
        //
        // // In order to avoid dependence of PrimaryGeneratorAction
        // // on DetectorConstruction class we get world volume
        // // from G4LogicalVolumeStore
        // //
        // //Beam Direction
        // TVector3 BeamDirection(0.,0.,1.);
        // BeamDirection=BeamDirection.Unit();
        //
        //
        // // Set gun position
        // G4ThreeVector PartPos(0.,0.,0.);   // 0 for middle of calorimeter
        //
        // TVector3 MomentumDirection;
        //
        //
        //
        //
        // //f.ls();
        // Double_t energy;
        // Double_t cos_theta;
        //
        // h1->GetRandom2(energy, cos_theta);
        // //std::cout<<cos_theta<<std::endl;
        //
        // fParticleGun->SetParticleEnergy(energy*GeV);
        //
        //
        // if(cos_theta>0) {
        //         Double_t theta = TMath::ACos(cos_theta);
        //
        //         Double_t radius = TMath::Tan(theta);
        //
        //         Double_t a=G4UniformRand()*TMath::TwoPi();
        //
        //         std::cout<<a<<std::endl;
        //         Double_t x=TMath::Cos(a)*TMath::Abs(radius);
        //         Double_t y=TMath::Sin(a)*TMath::Abs(radius);
        //         std::cout<<x<<":"<<y<<std::endl;
        //
        //         MomentumDirection.SetXYZ(x+BeamDirection.X(), y+BeamDirection.X(), BeamDirection.Z());
        //         MomentumDirection=MomentumDirection.Unit();
        //         std::cout<<">0"<<std::endl;
        //
        //
        //
        // }
        // else if(cos_theta==0) {
        //         Double_t theta = TMath::ACos(cos_theta);
        //
        //         Double_t a=G4UniformRand()*TMath::TwoPi();
        //
        //         std::cout<<a<<std::endl;
        //         Double_t x=TMath::Cos(a)*1;
        //         Double_t y=TMath::Sin(a)*1;
        //         std::cout<<x<<":"<<y<<std::endl;
        //
        //         MomentumDirection.SetXYZ(x+BeamDirection.X(), y+BeamDirection.Y(), 0);
        //         MomentumDirection=MomentumDirection.Unit();
        //         std::cout<<"==0"<<std::endl;
        //
        // }
        // else if(cos_theta<0) {
        //         Double_t theta = TMath::ACos(cos_theta);
        //
        //         Double_t radius = TMath::Tan(theta);
        //
        //         Double_t a=G4UniformRand()*TMath::TwoPi();
        //
        //         Double_t x=TMath::Cos(a)*TMath::Abs(radius);
        //         Double_t y=TMath::Sin(a)*TMath::Abs(radius);
        //         std::cout<<x<<":"<<y<<std::endl;
        //
        //         MomentumDirection.SetXYZ(x+BeamDirection.X(), y+BeamDirection.Y(), -BeamDirection.Z());
        //         MomentumDirection=MomentumDirection.Unit();
        //         std::cout<<"<0"<<std::endl;
        //
        // }
        //
        //
        // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(MomentumDirection.X(), MomentumDirection.Y(), MomentumDirection.Z()));
        //
        // fParticleGun->SetParticlePosition(PartPos);

        fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
