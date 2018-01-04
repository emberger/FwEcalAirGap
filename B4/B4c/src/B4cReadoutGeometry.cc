#include "B4cReadoutGeometry.hh"
#include "B4cCalorimeterSD.hh"


#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SDManager.hh"
#include "G4Box.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Material.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"
#include "G4Colour.hh"
#include <iostream>

#include "B4cDetParams.hh"

MyRO::MyRO() : G4VReadOutGeometry(){
}

MyRO::MyRO(G4String aString) : G4VReadOutGeometry(aString){



}

MyRO::~MyRO(){
}





G4VPhysicalVolume* MyRO::Build(){

								// Geometry parameters also for Detector Construction
								GetInst().SetfNofInnerLayers(40);  // in #
								GetInst().SetfNofOuterLayers(50);

								GetInst().SetcalorSizeXY(2000); // in mm

								GetInst().SetInnertileLenX(20); // in mm
								GetInst().SetInnertileLenY(20); // in mm

								GetInst().SetOutertileLenX(20);         // in mm
								GetInst().SetOutertileLenY(20);         // in mm

								GetInst().SetInnerabsoThickness(2);     // in mm
								GetInst().SetInnergapThickness(5);      // in mm

								GetInst().SetOuterabsoThickness(2); // in mm
								GetInst().SetOutergapThickness(5); // in mm

								GetInst().SetPvesselThickness(20); //in mm

								GetInst().SetInnerAbsFirst(true);    // true if structure should be |Abso-Gap|Abso-Gap|... false if structure should be |Gap-Abso|Gap-Abso|...
								GetInst().SetOuterAbsFirst(false);

								GetInst().SetWorldMult(10.);//in mm
								GetInst().InitDet(); // dont forget!!!




//	   auto ROworldSizeXY = 1.2 * GetInst().GetcalorSizeXY();
//	    auto ROworldSizeZ  = 1.2 * GetInst().GetcalorThickness();

								auto fCheckOverlaps=true;

								// A dummy material is used to fill the volumes of the readout geometry.
								// ( It will be allowed to set a NULL pointer in volumes of such virtual
								// division in future, since this material is irrelevant for tracking.)
								G4Material * dummyMat  = new G4Material("dummyMat", 1., 1.*g/mole, 1.*g/cm3);



								//Build Readout World

								auto ROWorldS = new G4Box("ROWorld",
																																		GetInst().GetWorldSizeXY()/2,
																																		GetInst().GetWorldSizeXY()/2,
																																		GetInst().GetWorldSizeZ()/2);

								auto ROWorldLog = new G4LogicalVolume(ROWorldS,
																																														dummyMat,
																																														"ROWorldLogical",
																																														0,
																																														0,
																																														0);


								G4PVPlacement *ROWorld = new G4PVPlacement(0,
																																																			G4ThreeVector(),
																																																			"ROWorldPhysical",
																																																			ROWorldLog,
																																																			0,
																																																			false,
																																																			0);

								//
								// Detector
								//
								auto RODetectorS
																= new G4Box("RODetector",                                        // its name
																												GetInst().GetcalorSizeXY()/2, GetInst().GetcalorSizeXY()/2, GetInst().GetDetectorThickness()/2);                                  // its size

								auto RODetectorLV
																= new G4LogicalVolume(
																RODetectorS,                                             // its solid
																dummyMat,                                          // its material
																"RODetector");                                           // its name

								G4ThreeVector ROsegD(0*mm,0*mm,(GetInst().GetDetectorThickness()/2)*mm);

								new G4PVPlacement(
																0,                                                        // no rotation
																ROsegD,                                          // at (0,0,0)
																RODetectorLV,                                                  // its logical volume
																"RODetector",                                            // its name
																ROWorldLog,                                                  // its mother  volume
																false,                                                    // no boolean operation
																1,                                                        // copy number
																fCheckOverlaps);                                       // checking overlaps

								//
								// CalorimeterInside
								//
								auto ROcalorimeterInsideS
																= new G4Box("ROCalorimeterInside",                                                // its name
																												GetInst().GetcalorSizeXY()/2, GetInst().GetcalorSizeXY()/2, GetInst().GetInnercalorThickness()/2);                                          // its size

								auto ROcalorInsideLV
																= new G4LogicalVolume(
																ROcalorimeterInsideS,                                                     // its solid
																dummyMat,                                                  // its material
																"ROCalorimeterInside");                                                   // its name

								G4ThreeVector ROsegI(0*mm,0*mm,(-GetInst().GetDetectorThickness()/2+GetInst().GetInnercalorThickness()/2)*mm);

								new G4PVPlacement(
																0,                                                                // no rotation
																ROsegI,                                                  // at (0,0,0)
																ROcalorInsideLV,                                                          // its logical volume
																"ROCalorimeterInside",                                                    // its name
																RODetectorLV,                                                          // its mother  volume
																false,                                                            // no boolean operation
																1,                                                                // copy number
																fCheckOverlaps);                                               // checking overlaps
								// checking overlaps

								//
								// ROCalorimeterOutside
								//
								auto ROcalorimeterOutsideS
																= new G4Box("ROCalorimeterOutside",                                                // its name
																												GetInst().GetcalorSizeXY()/2, GetInst().GetcalorSizeXY()/2, GetInst().GetOutercalorThickness()/2);                                          // its size

								auto ROcalorOutsideLV
																= new G4LogicalVolume(
																ROcalorimeterOutsideS,                                                     // its solid
																dummyMat,                                                  // its material
																"ROCalorimeterOutside");                                                   // its name

								G4ThreeVector ROsegO(0*mm,0*mm,(-GetInst().GetDetectorThickness()/2+GetInst().GetInnercalorThickness()+GetInst().GetPvesselThickness()+GetInst().GetOutercalorThickness()/2)*mm);

								new G4PVPlacement(
																0,                                                                // no rotation
																ROsegO,                                                  // at (0,0,0)
																ROcalorOutsideLV,                                                          // its logical volume
																"ROCalorimeterOutside",                                                    // its name
																RODetectorLV,                                                          // its mother  volume
																false,                                                            // no boolean operation
																1,                                                                // copy number
																fCheckOverlaps);                                               // checking overlaps







								//
								// InnerLayer
								//

								auto ROInnerlayerS
																= new G4Box("ROInnerLayer",             // its name
																												GetInst().GetcalorSizeXY()/2, GetInst().GetcalorSizeXY()/2, GetInst().GetInnerlayerThickness()/2);             //its size

								auto ROInnerlayerLV
																= new G4LogicalVolume(
																ROInnerlayerS,                        // its solid
																dummyMat,               // its material
																"ROInnerLayer");                      // its name

								new G4PVReplica(
																"ROInnerLayer",                       // its name
																ROInnerlayerLV,                       // its logical volume
																ROcalorInsideLV,                       // its mother
																kZAxis,                        // axis of replication
																GetInst().GetfNofInnerLayers(),                     // number of replica
																GetInst().GetInnerlayerThickness());               // witdth of replica


								//-------------------------------
								//build Inner calorimeter readout cells
								//-------------------------------

								//build Gap
								auto InnerGapS
																=new G4Box("InnerGap", GetInst().GetcalorSizeXY()/2, GetInst().GetcalorSizeXY()/2, GetInst().GetInnergapThickness()/2);

								auto InnerGapLV
																=new G4LogicalVolume(InnerGapS,
																																					dummyMat,
																																					"InnerGap");
								if(GetInst().GetInnerAbsFirst()) {

																new G4PVPlacement(0,
																																		G4ThreeVector(0., 0., GetInst().GetInnerabsoThickness()/2 ),
																																		InnerGapLV,
																																		"Gap",
																																		ROInnerlayerLV,
																																		false,
																																		0,
																																		true);


								}
								else{
																new G4PVPlacement(0,
																																		G4ThreeVector(0., 0., -GetInst().GetInnerabsoThickness()/2 ),
																																		InnerGapLV,
																																		"InnerGap",
																																		ROInnerlayerLV,
																																		false,
																																		0,
																																		true);

								}
								//build strips parallel to X axis

								auto InnerStripS
																=new G4Box("InnerStrip", GetInst().GetcalorSizeXY()/2,GetInst().GetInnertileLenY()/2,GetInst().GetInnergapThickness()/2);

								auto InnerStripLV
																=new G4LogicalVolume(InnerStripS,
																																					dummyMat,
																																					"InnerStrip");

								new G4PVReplica("InnerStrip",
																								InnerStripLV, //logical volume
																								InnerGapLV, //mother volume
																								kYAxis,
																								GetInst().GetnofInnerTilesY(),
																								GetInst().GetInnertileLenY());


								// build cells in Xaxis strip

								auto InnerCellS
																=new G4Box("InnerCell", GetInst().GetInnertileLenX()/2, GetInst().GetInnertileLenY()/2, GetInst().GetInnergapThickness()/2);

								auto InnerCellLV
																=new G4LogicalVolume(InnerCellS,
																																					dummyMat,
																																					"InnerCell");

								new G4PVReplica("InnerCell",
																								InnerCellLV, //logical volume
																								InnerStripLV, //mother volume
																								kXAxis,
																								GetInst().GetnofInnerTilesX(),
																								GetInst().GetInnertileLenX());





								//
								// OuterLayer
								//

								auto ROOuterlayerS
																= new G4Box("ROOuterLayer",                             // its name
																												GetInst().GetcalorSizeXY()/2, GetInst().GetcalorSizeXY()/2, GetInst().GetOuterlayerThickness()/2);                             //its size

								auto ROOuterlayerLV
																= new G4LogicalVolume(
																ROOuterlayerS,                                        // its solid
																dummyMat,                               // its material
																"ROOuterLayer");                                      // its name

								new G4PVReplica(
																"ROOuterLayer",                                       // its name
																ROOuterlayerLV,                                       // its logical volume
																ROcalorOutsideLV,                                       // its mother
																kZAxis,                                        // axis of replication
																GetInst().GetfNofOuterLayers(),                                     // number of replica
																GetInst().GetOuterlayerThickness());                               // witdth of replica
								//-------------------------------
								//build Outer calorimeter readout cells
								//-------------------------------

								//build Gap
								auto OuterGapS
																=new G4Box("OuterGap", GetInst().GetcalorSizeXY()/2, GetInst().GetcalorSizeXY()/2, GetInst().GetOutergapThickness()/2);

								auto OuterGapLV
																=new G4LogicalVolume(OuterGapS,
																																					dummyMat,
																																					"OuterGap");
								if(GetInst().GetOuterAbsFirst()) {

																new G4PVPlacement(0,
																																		G4ThreeVector(0., 0., GetInst().GetOuterabsoThickness()/2 ),
																																		OuterGapLV,
																																		"Gap",
																																		ROOuterlayerLV,
																																		false,
																																		0,
																																		true);


								}
								else{
																new G4PVPlacement(0,
																																		G4ThreeVector(0., 0., -GetInst().GetOuterabsoThickness()/2 ),
																																		OuterGapLV,
																																		"OuterGap",
																																		ROOuterlayerLV,
																																		false,
																																		0,
																																		true);

								}
								//build strips parallel to X axis

								auto OuterStripS
																=new G4Box("OuterStrip", GetInst().GetcalorSizeXY()/2,GetInst().GetOutertileLenY()/2,GetInst().GetOutergapThickness()/2);

								auto OuterStripLV
																=new G4LogicalVolume(OuterStripS,
																																					dummyMat,
																																					"OuterStrip");

								new G4PVReplica("OuterStrip",
																								OuterStripLV,                 //logical volume
																								OuterGapLV,                 //mother volume
																								kYAxis,
																								GetInst().GetnofOuterTilesY(),
																								GetInst().GetOutertileLenY());


								// build cells in Xaxis strip

								auto OuterCellS
																=new G4Box("OuterCell", GetInst().GetOutertileLenX()/2, GetInst().GetOutertileLenY()/2, GetInst().GetOutergapThickness()/2);

								auto OuterCellLV
																=new G4LogicalVolume(OuterCellS,
																																					dummyMat,
																																					"OuterCell");

								new G4PVReplica("OuterCell",
																								OuterCellLV,                 //logical volume
																								OuterStripLV,                 //mother volume
																								kXAxis,
																								GetInst().GetnofOuterTilesX(),
																								GetInst().GetOutertileLenX());











								auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,0.0,1.0));
								simpleBoxVisAtt->SetVisibility(true);
								InnerCellLV->SetVisAttributes(simpleBoxVisAtt);
								OuterCellLV->SetVisAttributes(simpleBoxVisAtt);




								auto dummy
																= new B4cCalorimeterSD("dummy", "dummyCollection");

								InnerCellLV->SetSensitiveDetector(dummy);
								OuterCellLV->SetSensitiveDetector(dummy);


								return ROWorld;




}
