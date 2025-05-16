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
//
/// \file B1/src/DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4CutTubs.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

namespace B1
{

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  G4VPhysicalVolume *DetectorConstruction::Construct()
  {
    // Get nist material manager
    G4NistManager *nist = G4NistManager::Instance();

    // Tube material
    //
    G4Material *tube_mat = nist->FindOrBuildMaterial("G4_Cu");

    // Option to switch on/off checking of volumes overlaps
    //
    G4bool checkOverlaps = true;

    //
    // World
    //
    G4double world_sizeXY = 400 * cm;
    G4double world_sizeZ = 400 * cm;
    G4Material *world_mat = nist->FindOrBuildMaterial("G4_AIR");

    auto solidWorld = new G4Box("World",                                                    // its name
                                0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ); // its size

    auto logicWorld = new G4LogicalVolume(solidWorld, // its solid
                                          world_mat,  // its material
                                          "World");   // its name

    auto physWorld = new G4PVPlacement(nullptr,         // no rotation
                                       G4ThreeVector(), // at (0,0,0)
                                       logicWorld,      // its logical volume
                                       "World",         // its name
                                       nullptr,         // its mother  volume
                                       false,           // no boolean operation
                                       0,               // copy number
                                       checkOverlaps);  // overlaps checking

    //
    // Tube
    //
    G4double R_inner = 30 * mm;
    G4double R_outer = 36 * mm;
    G4double L_tube = 51.5 * mm;
    G4double phi_0 = 0;
    G4double phi_1 = 2 * M_PI;
    G4ThreeVector n_bot = G4ThreeVector(0., 0., -1.);
    G4ThreeVector n_top = G4ThreeVector(0., 1., 1.);

    auto solidTube = new G4CutTubs("MyTube",                                              // its name
                                   R_inner, R_outer, L_tube, phi_0, phi_1, n_bot, n_top); // its size

    auto logicTube = new G4LogicalVolume(solidTube, // its solid
                                         tube_mat,  // its material
                                         "MyTube"); // its name

    new G4PVPlacement(nullptr,         // no rotation
                      G4ThreeVector(), // at (0,0,0)
                      logicTube,       // its logical volume
                      "MyTube",        // its name
                      logicWorld,      // its mother  volume
                      false,           // no boolean operation
                      0,               // copy number
                      checkOverlaps);  // overlaps checking

    G4VisAttributes *copperVisAttributes = new G4VisAttributes(G4Colour(0.7, 0.4, 0.1));
    copperVisAttributes->SetForceWireframe(true); // Translucenty
    logicTube->SetVisAttributes(copperVisAttributes); // Setting copper colour for copper

    //
    // Target 1 layer
    //
    G4Material *target1_mat = nist->FindOrBuildMaterial("G4_W");
    G4ThreeVector pos1 = G4ThreeVector(0, 0, 50.5 * mm );

    // Blinchik wolfram
    G4double R_inner_target = 0 * mm;
    G4double R_outer_target = 30 * mm;
    G4double L_tube_target = 0.5 * mm;
    G4double phi_0_target = 0;
    G4double phi_1_target = 2 * M_PI;
    G4ThreeVector n_bot_target = G4ThreeVector(0., -1., -1.);
    G4ThreeVector n_top_target = G4ThreeVector(0., 1., 1.);

    auto solidTarget1 = new G4CutTubs("WolfTarget", R_inner_target, R_outer_target, L_tube_target, phi_0_target, phi_1_target, n_bot_target, n_top_target);

    auto logicTarget1 = new G4LogicalVolume(solidTarget1,  // its solid
                                            target1_mat,   // its material
                                            "WolfTarget"); // its name

    new G4PVPlacement(nullptr,        // no rotation
                      pos1,           // at position
                      logicTarget1,   // its logical volume
                      "WolfTarget",   // its name
                      logicWorld,     // its mother  volume
                      false,          // no boolean operation
                      0,              // copy number
                      checkOverlaps); // overlaps checking

    G4VisAttributes *wolframVisAttributes = new G4VisAttributes(G4Colour(0.75, 0.75, 0.75));
    // copperVisAttributes->SetForceWireframe(true); // Translucenty
    logicTarget1->SetVisAttributes(wolframVisAttributes); // Setting copper colour for copper

    //
    // Target 2 layer
    //
    G4Material *target2_mat = nist->FindOrBuildMaterial("G4_Cu");
    G4ThreeVector pos2 = G4ThreeVector(0, 0, 51 * mm);

    // Blinchik cuprum

    auto solidTarget2 = new G4CutTubs("CuTarget", R_inner_target, R_outer_target, L_tube_target, phi_0_target, phi_1_target, n_bot_target, n_top_target);

    auto logicTarget2 = new G4LogicalVolume(solidTarget2, // its solid
                                            target2_mat,  // its material
                                            "CuTarget");  // its name

    new G4PVPlacement(nullptr,        // no rotation
                      pos2,           // at position
                      logicTarget2,   // its logical volume
                      "CuTarget",     // its name
                      logicWorld,     // its mother  volume
                      false,          // no boolean operation
                      0,              // copy number
                      checkOverlaps); // overlaps checking

    // G4VisAttributes * copperVisAttributes = new G4VisAttributes(G4Colour(0.7, 0.4, 0.1));
    // copperVisAttributes->SetForceWireframe(true); // Translucenty
    logicTarget2->SetVisAttributes(copperVisAttributes); // Setting copper colour for copper

    //
    // Detector volume
    //
    G4ThreeVector pos3 = G4ThreeVector(0, 0, L_tube);
    // Rotation
    G4RotationMatrix *rot3 = new G4RotationMatrix();
    rot3->rotateX(-45 * deg);
    rot3->rotateY(0 * deg);
    rot3->rotateZ(0 * deg);

    G4double Rmin_det = 36 * mm;
    G4double Rmax_det = 50 * mm;
    G4double SPhi_det = 0 * deg;
    G4double DPhi_det = 180 * deg;
    G4double STheta_det = 0 * deg;
    G4double DTheta_det = 180 * deg;

    auto solidDetector = new G4Sphere("AirDetector", Rmin_det, Rmax_det, SPhi_det, DPhi_det, STheta_det, DTheta_det);

    auto logicDetector = new G4LogicalVolume(solidDetector,  // its solid
                                             world_mat,      // its material
                                             "AirDetector"); // its name

    new G4PVPlacement(rot3,           // rotation
                      pos3,           // at position
                      logicDetector,  // its logical volume
                      "AirDetector",  // its name
                      logicWorld,     // its mother  volume
                      false,          // no boolean operation
                      0,              // copy number
                      checkOverlaps); // overlaps checking

    G4VisAttributes *detectorVisAttributes = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    detectorVisAttributes->SetForceWireframe(true); // Translucenty
    logicDetector->SetVisAttributes(detectorVisAttributes);

    //
    // Detector volume 0.7 * m
    //
    G4ThreeVector pos4 = G4ThreeVector(0, 0, L_tube + 0.7 * m);
    // No Rotation

    G4double R_inner_det07 = 0 * mm;
    G4double R_outer_det07 = 200 * cm;
    G4double L_tube_det07 = 1 * mm;
    G4double phi_0_det07 = 0;
    G4double phi_1_det07 = 2 * M_PI;

    auto solidDetector07 = new G4Tubs("AirDetector_07", R_inner_det07, R_outer_det07, L_tube_det07, phi_0_det07, phi_1_det07);

    auto logicDetector07 = new G4LogicalVolume(solidDetector07,   // its solid
                                               world_mat,         // its material
                                               "AirDetector_07"); // its name

    new G4PVPlacement(nullptr,          // rotation
                      pos4,             // at position
                      logicDetector07,  // its logical volume
                      "AirDetector_07", // its name
                      logicWorld,       // its mother  volume
                      false,            // no boolean operation
                      0,                // copy number
                      checkOverlaps);   // overlaps checking

    G4VisAttributes *detectorVisAttributes_07 = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    detectorVisAttributes_07->SetForceWireframe(true); // Translucenty
    logicDetector07->SetVisAttributes(detectorVisAttributes_07);

    //
    // Detector volume 1 * m
    //
    G4ThreeVector pos5 = G4ThreeVector(0, 0, L_tube + 1. * m);
    // No Rotation

    G4double R_inner_det10 = 0 * mm;
    G4double R_outer_det10 = 200 * cm;
    G4double L_tube_det10 = 1 * mm;
    G4double phi_0_det10 = 0;
    G4double phi_1_det10 = 2 * M_PI;

    auto solidDetector10 = new G4Tubs("AirDetector_10", R_inner_det10, R_outer_det10, L_tube_det10, phi_0_det10, phi_1_det10);

    auto logicDetector10 = new G4LogicalVolume(solidDetector10,   // its solid
                                               world_mat,         // its material
                                               "AirDetector_10"); // its name

    new G4PVPlacement(nullptr,          // rotation
                      pos5,             // at position
                      logicDetector10,  // its logical volume
                      "AirDetector_10", // its name
                      logicWorld,       // its mother  volume
                      false,            // no boolean operation
                      0,                // copy number
                      checkOverlaps);   // overlaps checking

    G4VisAttributes *detectorVisAttributes_10 = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    detectorVisAttributes_10->SetForceWireframe(true); // Translucenty
    logicDetector10->SetVisAttributes(detectorVisAttributes_10);

    // Set AirDetector as scoring volume
    //
    fScoringVolume = logicDetector;
    fScoringPan07 = logicDetector07;
    fScoringPan10 = logicDetector10;

    //
    // always return the physical World
    //
    return physWorld;
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
