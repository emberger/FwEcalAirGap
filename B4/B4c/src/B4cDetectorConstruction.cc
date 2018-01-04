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
// $Id: B4cDetectorConstruction.cc 101905 2016-12-07 11:34:39Z gunter $
//
/// \file B4cDetectorConstruction.cc
/// \brief Implementation of the B4cDetectorConstruction class

#include "B4cDetectorConstruction.hh"
#include "B4cCalorimeterSD.hh"
#include "B4cReadoutGeometry.hh"



#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4RotationMatrix.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "B4cDetParams.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//instance of Readout Geometry
MyRO * ROGeom = new MyRO("readout");


G4ThreadLocal
G4GlobalMagFieldMessenger* B4cDetectorConstruction::fMagFieldMessenger = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cDetectorConstruction::B4cDetectorConstruction()
        : G4VUserDetectorConstruction(),
        fStepLimit(NULL),
        fCheckOverlaps(true),
        fNofLayers(-1)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cDetectorConstruction::~B4cDetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4cDetectorConstruction::Construct()
{

        //build Readout Geometry
        ROGeom->BuildROGeometry();

        // Define materials
        DefineMaterials();

        // Define volumes
        return DefineVolumes();




}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cDetectorConstruction::DefineMaterials()
{
        // Lead material defined using NIST Manager
        auto nistManager = G4NistManager::Instance();
        nistManager->FindOrBuildMaterial("G4_Pb");
        nistManager->FindOrBuildMaterial("G4_Cu");
        nistManager->FindOrBuildMaterial("G4_POLYSTYRENE");
        nistManager->FindOrBuildMaterial("G4_Ti");
        nistManager->FindOrBuildMaterial("G4_AIR");


        //BGO material
        G4double BGOdensity=7.13*g/cm3;
        G4int BGOcomponents=3;
        G4Material *  BGO = new G4Material("BGO", BGOdensity, BGOcomponents);

        BGO->AddElement(nistManager->FindOrBuildElement(83),4);
        BGO->AddElement(nistManager->FindOrBuildElement(32),3);
        BGO->AddElement(nistManager->FindOrBuildElement(8),12);

        // PEN Material
        G4double Pdensity;
        Pdensity=1.36*g/cm3;
        G4int Pcomponents=3;
        G4Material* pen=new G4Material("pen",Pdensity,Pcomponents);

        pen->AddElement(nistManager->FindOrBuildElement(6),14);
        pen->AddElement(nistManager->FindOrBuildElement(1),10);
        pen->AddElement(nistManager->FindOrBuildElement(8),4);

        // LYSO Crystal

        G4double LYSOdensity = 7.1*g/cm3;
        G4int LYSOcomponents = 4;
        G4Material * LYSO = new G4Material("LYSO", LYSOdensity, LYSOcomponents);

        LYSO->AddElement(nistManager->FindOrBuildElement(71), 9);
        LYSO->AddElement(nistManager->FindOrBuildElement(39), 1);
        LYSO->AddElement(nistManager->FindOrBuildElement(14), 5);
        LYSO->AddElement(nistManager->FindOrBuildElement(8), 25);

        //TBBPA material
        G4double TBBPAdensity=2.12*g/cm3;
        G4int TBBPAcomponents=4;
        G4Material * TBBPA = new G4Material("TBBPA", TBBPAdensity, TBBPAcomponents);

        TBBPA->AddElement(nistManager->FindOrBuildElement(6),15);
        TBBPA->AddElement(nistManager->FindOrBuildElement(1),12);
        TBBPA->AddElement(nistManager->FindOrBuildElement(35),4);
        TBBPA->AddElement(nistManager->FindOrBuildElement(8),2);

        //Fiberglass Material
        G4double FIBERdensity=2.5*g/cm3;
        G4int FIBERcomponents=3;
        G4Material * fiberglass= new G4Material("fiberglass", FIBERdensity, FIBERcomponents);

        fiberglass->AddElement(nistManager->FindOrBuildElement(8),5);
        fiberglass->AddElement(nistManager->FindOrBuildElement(14),1);
        fiberglass->AddElement(nistManager->FindOrBuildElement(5),2);

        //FR4 Material
        G4double FR4density=1.8*g/cm3;
        G4int FR4components=3;
        G4Material * FR4= new G4Material("FR4", FR4density, FR4components);

        FR4->AddMaterial(fiberglass, 0.92);
        FR4->AddMaterial(TBBPA, 0.07);
        FR4->AddMaterial(nistManager->FindOrBuildMaterial("G4_Cu"),0.01);

        //Steel Material
        G4double SteelDensity=8.027*g/cm3;
        G4int Steelcomponents=5;
        G4Material * Steel=new G4Material("Steel", SteelDensity, Steelcomponents);

        Steel->AddMaterial(nistManager->FindOrBuildMaterial("G4_Fe"),0.709 );
        Steel->AddMaterial(nistManager->FindOrBuildMaterial("G4_Cr"),0.18 );
        Steel->AddMaterial(nistManager->FindOrBuildMaterial("G4_Ni"), 0.10);
        Steel->AddMaterial(nistManager->FindOrBuildMaterial("G4_Mn"), 0.01);
        Steel->AddMaterial(nistManager->FindOrBuildMaterial("G4_C"), 0.001);




        // Liquid argon material
        G4double a; // mass of a mole;
        G4double z; // z=mean number of protons;
        G4double density;
//  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
//         // The argon by NIST Manager is a gas with a different density

        // Vacuum
        new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                       kStateGas, 2.73*kelvin, 3.e-18*pascal);

        // Print materials
        G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4cDetectorConstruction::DefineVolumes()
{

        auto defaultMaterial = G4Material::GetMaterial("Galactic");
        auto absorberMaterial = G4Material::GetMaterial("G4_Pb");
        auto gapMaterial = G4Material::GetMaterial("G4_POLYSTYRENE");
        auto airMaterial = G4Material::GetMaterial("G4_AIR");
        auto PVesselMaterial = G4Material::GetMaterial("Steel");

        if ( !defaultMaterial || !absorberMaterial || !gapMaterial ) {
                G4ExceptionDescription msg;
                msg << "Cannot retrieve materials already defined.";
                G4Exception("B4DetectorConstruction::DefineVolumes()",
                            "MyCode0001", FatalException, msg);
        }

        //
        // World
        //
        auto worldS
                = new G4Box("World", // its name
                            GetInst().GetWorldSizeXY()/2, GetInst().GetWorldSizeXY()/2, GetInst().GetWorldSizeZ()/2); // its size

        auto worldLV
                = new G4LogicalVolume(
                worldS,            // its solid
                defaultMaterial,   // its material
                "World");          // its name

        auto worldPV
                = new G4PVPlacement(
                0,                 // no rotation
                G4ThreeVector(),   // at (0,0,0)
                worldLV,           // its logical volume
                "World",           // its name
                0,                 // its mother  volume
                false,             // no boolean operation
                0,                 // copy number
                fCheckOverlaps);   // checking overlaps

        //
        // Detector
        //
        auto DetectorS
                = new G4Box("Detector", // its name
                            GetInst().GetcalorSizeXY()/2, GetInst().GetcalorSizeXY()/2, GetInst().GetDetectorThickness()/2); // its size

        auto DetectorLV
                = new G4LogicalVolume(
                DetectorS,      // its solid
                defaultMaterial,   // its material
                "Detector");    // its name

        G4ThreeVector segD(0*mm,0*mm,(GetInst().GetDetectorThickness()/2)*mm);

        new G4PVPlacement(
                0,                 // no rotation
                segD,
                DetectorLV,           // its logical volume
                "Detector",     // its name
                worldLV,           // its mother  volume
                false,             // no boolean operation
                1,                 // copy number
                fCheckOverlaps);// checking overlaps

        //
        // CalorimeterInside
        //
        auto calorimeterInsideS
                = new G4Box("CalorimeterInside",         // its name
                            GetInst().GetcalorSizeXY()/2, GetInst().GetcalorSizeXY()/2, GetInst().GetInnercalorThickness()/2);         // its size

        auto calorInsideLV
                = new G4LogicalVolume(
                calorimeterInsideS,              // its solid
                defaultMaterial,           // its material
                "CalorimeterInside");            // its name

        G4ThreeVector segI(0*mm,0*mm,(-GetInst().GetDetectorThickness()/2+GetInst().GetInnercalorThickness()/2)*mm);

        new G4PVPlacement(
                0,                         // no rotation
                segI,
                calorInsideLV,                   // its logical volume
                "CalorimeterInside",             // its name
                DetectorLV,                   // its mother  volume
                false,                     // no boolean operation
                1,                         // copy number
                fCheckOverlaps);        // checking overlaps
        //
        // PressureVessel
        //
        auto PressureVesselS
                = new G4Box("PressureVessel",                         // its name
                            GetInst().GetcalorSizeXY()/2, GetInst().GetcalorSizeXY()/2, GetInst().GetPvesselThickness()/2);                         // its size

        auto PressureVesselLV
                = new G4LogicalVolume(
                PressureVesselS,                              // its solid
                PVesselMaterial,                           // its material
                "PressureVessel");                            // its name

        G4ThreeVector segP(0*mm,0*mm,(-GetInst().GetDetectorThickness()/2+GetInst().GetInnercalorThickness()+GetInst().GetPvesselThickness()/2)*mm);

        new G4PVPlacement(
                0,                                         // no rotation
                segP,
                PressureVesselLV,                                   // its logical volume
                "PressureVessel",                             // its name
                DetectorLV,                                   // its mother  volume
                false,                                     // no boolean operation
                1,                                         // copy number
                fCheckOverlaps);                        // checking overlaps

        //
        // CalorimeterOutside
        //
        auto calorimeterOutsideS
                = new G4Box("CalorimeterOutside",         // its name
                            GetInst().GetcalorSizeXY()/2, GetInst().GetcalorSizeXY()/2, GetInst().GetOutercalorThickness()/2);         // its size

        auto calorOutsideLV
                = new G4LogicalVolume(
                calorimeterOutsideS,              // its solid
                defaultMaterial,           // its material
                "CalorimeterOutside");            // its name

        G4ThreeVector segO(0*mm,0*mm,(-GetInst().GetDetectorThickness()/2+GetInst().GetInnercalorThickness()+GetInst().GetPvesselThickness()+GetInst().GetOutercalorThickness()/2)*mm);

        new G4PVPlacement(
                0,                         // no rotation
                segO,
                calorOutsideLV,                   // its logical volume
                "CalorimeterOutside",             // its name
                DetectorLV,                   // its mother  volume
                false,                     // no boolean operation
                1,                         // copy number
                fCheckOverlaps);        // checking overlaps
        //
        //AirGap
        //

        // auto AirGapS
        //         = new G4Box("AirGap", // its name
        //                     GetInst().GetcalorSizeXY()/2, GetInst().GetcalorSizeXY()/2, GetInst().GetairgapThickness()/2); // its size
        //
        // auto AirGapLV
        //         = new G4LogicalVolume(
        //         AirGapS,         // its solid
        //         airMaterial, // its material
        //         "AirGapLV");         // its name
        //
        //
        // new G4PVPlacement(
        //         0,                 // no rotation
        //         G4ThreeVector(0., 0., -(GetInst().GetInnergapThickness()/2)),         // its position
        //         AirGapLV,         // its logical volume
        //         "AirGap",            // its name
        //         AirGapLV,           // its mother  volume
        //         false,             // no boolean operation
        //         0,                 // copy number
        //         fCheckOverlaps);         // checking overlaps






        //
        // InnerLayer
        //

        auto InnerlayerS
                = new G4Box("InnerLayer", // its name
                            GetInst().GetcalorSizeXY()/2, GetInst().GetcalorSizeXY()/2, GetInst().GetInnerlayerThickness()/2); //its size

        auto InnerlayerLV
                = new G4LogicalVolume(
                InnerlayerS,            // its solid
                defaultMaterial,   // its material
                "InnerLayer");          // its name

        new G4PVReplica(
                "InnerLayer",           // its name
                InnerlayerLV,           // its logical volume
                calorInsideLV,           // its mother
                kZAxis,            // axis of replication
                GetInst().GetfNofInnerLayers(),         // number of replica
                GetInst().GetInnerlayerThickness());   // witdth of replica



        //
        // InnerAbsorber
        //
        auto InnerabsorberS
                = new G4Box("InnerAbso", // its name
                            GetInst().GetcalorSizeXY()/2, GetInst().GetcalorSizeXY()/2, GetInst().GetInnerabsoThickness()/2); // its size

        auto InnerabsorberLV
                = new G4LogicalVolume(
                InnerabsorberS,         // its solid
                absorberMaterial,  // its material
                "InnerAbsoLV");         // its name

        if(GetInst().GetInnerAbsFirst()) {
                new G4PVPlacement(
                        0,         // no rotation
                        G4ThreeVector(0., 0., -(GetInst().GetInnergapThickness()/2)), // its position
                        InnerabsorberLV, // its logical volume
                        "InnerAbso",    // its name
                        InnerlayerLV,   // its mother  volume
                        false,     // no boolean operation
                        0,         // copy number
                        fCheckOverlaps); // checking overlaps
        }
        else{
                new G4PVPlacement(
                        0,   // no rotation
                        G4ThreeVector(0., 0., (GetInst().GetInnergapThickness()/2)), // its position
                        InnerabsorberLV, // its logical volume
                        "InnerAbso", // its name
                        InnerlayerLV, // its mother  volume
                        false, // no boolean operation
                        0,   // copy number
                        fCheckOverlaps); // checking overlaps
        }
        //
        // InnerGap
        //
        auto InnergapS
                = new G4Box("InnerGap", // its name
                            GetInst().GetcalorSizeXY()/2, GetInst().GetcalorSizeXY()/2, GetInst().GetInnergapThickness()/2); // its size

        auto InnergapLV
                = new G4LogicalVolume(
                InnergapS,              // its solid
                gapMaterial,       // its material
                "InnerGapLV");          // its name

        if(GetInst().GetInnerAbsFirst()) {
                new G4PVPlacement(
                        0,         // no rotation
                        G4ThreeVector(0., 0., GetInst().GetInnerabsoThickness()/2), // its position
                        InnergapLV,     // its logical volume
                        "InnerGap",     // its name
                        InnerlayerLV,   // its mother  volume
                        false,     // no boolean operation
                        0,         // copy number
                        fCheckOverlaps); // checking overlaps
        }
        else{
                new G4PVPlacement(
                        0,   // no rotation
                        G4ThreeVector(0., 0., -GetInst().GetInnerabsoThickness()/2), // its position
                        InnergapLV, // its logical volume
                        "InnerGap", // its name
                        InnerlayerLV, // its mother  volume
                        false, // no boolean operation
                        0,   // copy number
                        fCheckOverlaps); // checking overlaps

        }

        //
        // OuterLayer
        //

        auto OuterlayerS
                = new G4Box("OuterLayer", // its name
                            GetInst().GetcalorSizeXY()/2, GetInst().GetcalorSizeXY()/2, GetInst().GetOuterlayerThickness()/2); //its size

        auto OuterlayerLV
                = new G4LogicalVolume(
                OuterlayerS,            // its solid
                defaultMaterial,   // its material
                "OuterLayer");          // its name

        new G4PVReplica(
                "OuterLayer",           // its name
                OuterlayerLV,           // its logical volume
                calorOutsideLV,           // its mother
                kZAxis,            // axis of replication
                GetInst().GetfNofOuterLayers(),         // number of replica
                GetInst().GetOuterlayerThickness());   // witdth of replica



        //
        // OuterAbsorber
        //

        auto OuterabsorberS
                = new G4Box("OuterAbso",                 // its name
                            GetInst().GetcalorSizeXY()/2, GetInst().GetcalorSizeXY()/2, GetInst().GetOuterabsoThickness()/2);                 // its size

        auto OuterabsorberLV
                = new G4LogicalVolume(
                OuterabsorberS,                         // its solid
                absorberMaterial,                  // its material
                "OuterAbsoLV");                         // its name

        if(GetInst().GetOuterAbsFirst()) {
                new G4PVPlacement(
                        0,                         // no rotation
                        G4ThreeVector(0., 0., -(GetInst().GetOutergapThickness()/2)),                 // its position
                        OuterabsorberLV,                 // its logical volume
                        "Abso",                    // its name
                        OuterlayerLV,                   // its mother  volume
                        false,                     // no boolean operation
                        0,                         // copy number
                        fCheckOverlaps);                 // checking overlaps
        }
        else{
                new G4PVPlacement(
                        0,                   // no rotation
                        G4ThreeVector(0., 0., (GetInst().GetOutergapThickness()/2)),                 // its position
                        OuterabsorberLV,                 // its logical volume
                        "OuterAbso",                 // its name
                        OuterlayerLV,                 // its mother  volume
                        false,                 // no boolean operation
                        0,                   // copy number
                        fCheckOverlaps);                 // checking overlaps
        }

        //
        // OuterGap
        //

        auto OutergapS
                = new G4Box("OuterGap",                 // its name
                            GetInst().GetcalorSizeXY()/2, GetInst().GetcalorSizeXY()/2, GetInst().GetOutergapThickness()/2);                 // its size

        auto OutergapLV
                = new G4LogicalVolume(
                OutergapS,                              // its solid
                gapMaterial,                       // its material
                "OuterGapLV");                          // its name

        if(GetInst().GetOuterAbsFirst()) {
                new G4PVPlacement(
                        0,                         // no rotation
                        G4ThreeVector(0., 0., GetInst().GetOuterabsoThickness()/2),                 // its position
                        OutergapLV,                     // its logical volume
                        "OuterGap",                     // its name
                        OuterlayerLV,                   // its mother  volume
                        false,                     // no boolean operation
                        0,                         // copy number
                        fCheckOverlaps);                 // checking overlaps
        }
        else{
                new G4PVPlacement(
                        0,                   // no rotation
                        G4ThreeVector(0., 0., -GetInst().GetOuterabsoThickness()/2),                 // its position
                        OutergapLV,                 // its logical volume
                        "OuterGap",                 // its name
                        OuterlayerLV,                 // its mother  volume
                        false,                 // no boolean operation
                        0,                   // copy number
                        fCheckOverlaps);                 // checking overlaps

        }



        //
        // print parameters
        //

        // G4cout
        //         << G4endl
        //         << "------------------------------------------------------------" << G4endl
        //         << "---> The calorimeter is " << GetInst().GetfNofLayers() << " layers of: [ "
        //         << GetInst().GetabsoThickness()/mm << "mm of " << absorberMaterial->GetName()
        //         << " + "
        //         << GetInst().GetgapThickness()/mm << "mm of " << gapMaterial->GetName() << " ] " << G4endl
        //         << "------------------------------------------------------------" << G4endl;

        //
        // Visualization attributes
        //
        worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

        auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
        simpleBoxVisAtt->SetVisibility(true);
        DetectorLV->SetVisAttributes(simpleBoxVisAtt);

        //
        // Always return the physical World
        //
        return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cDetectorConstruction::ConstructSDandField()
{
        G4SDManager::GetSDMpointer()->SetVerboseLevel(2);

        //
        // Sensitive detectors
        //
//  auto absoSD
//    = new B4cCalorimeterSD("AbsorberSD", "AbsorberHitsCollection", GetInst().GetfNofLayers(), GetInst().GettilesPerLayer(), GetInst().GetnofTilesX());
//  G4SDManager::GetSDMpointer()->AddNewDetector(absoSD);
//  SetSensitiveDetector("AbsoLV",absoSD);

        auto gapSD
                = new B4cCalorimeterSD("GapSD", "GapHitsCollection");

        G4SDManager::GetSDMpointer()->AddNewDetector(gapSD);

        SetSensitiveDetector("InnerGapLV",gapSD);
        SetSensitiveDetector("OuterGapLV",gapSD);

        gapSD->SetROgeometry(ROGeom);

        //G4cout<<ROGeom->GetName()<<G4endl;

        //
        // Magnetic field
        //
        // Create global magnetic field messenger.
        // Uniform magnetic field is then created automatically if
        // the field value is not zero.
        G4ThreeVector fieldValue;
        fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
        fMagFieldMessenger->SetVerboseLevel(1);

        // Register the field messenger for deleting
        G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
