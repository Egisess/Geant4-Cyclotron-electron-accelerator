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

    G4double L_tube = 50 * mm + 30 * mm;
    G4VisAttributes *copperVisAttributes = new G4VisAttributes(G4Colour(0.7, 0.4, 0.1));

    //
    // Target 1 layer
    //
    G4Material *target1_mat = nist->FindOrBuildMaterial("G4_W");
    G4ThreeVector pos1 = G4ThreeVector(0, 0, 50.0 * mm);

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
    // Sphere detector volume 1
    //
    G4ThreeVector pos4 = pos1;
    // Rotation
    G4RotationMatrix *rot4 = new G4RotationMatrix();
    rot4->rotateX(0 * deg);
    rot4->rotateY(0 * deg);
    rot4->rotateZ(0 * deg);

    G4double Rmin_sph1 = 100 * mm;
    G4double Rmax_sph1 = 125 * mm;
    G4double SPhi_sph1 = 0 * deg;
    G4double DPhi_sph1 = 360 * deg;
    G4double STheta_sph1 = 0 * deg;
    G4double DTheta_sph1 = 360 * deg;

    auto solidSphere1 = new G4Sphere("AirDetector", Rmin_sph1, Rmax_sph1, SPhi_sph1, DPhi_sph1, STheta_sph1, DTheta_sph1);

    auto logicSphere1 = new G4LogicalVolume(solidSphere1,  // its solid
                                            world_mat,     // its material
                                            "AirSphere1"); // its name

    new G4PVPlacement(rot4,           // rotation
                      pos4,           // at position
                      logicSphere1,   // its logical volume
                      "AirSphere1",   // its name
                      logicWorld,     // its mother  volume
                      false,          // no boolean operation
                      0,              // copy number
                      checkOverlaps); // overlaps checking

    G4VisAttributes *sphericalVisAttributes = new G4VisAttributes(G4Colour(0.25, 0.25, 0.25));
    sphericalVisAttributes->SetForceWireframe(true); // Translucenty
    logicSphere1->SetVisAttributes(sphericalVisAttributes);

    //
    // Sphere detector volume 2
    //

    G4double Rmin_sph2 = 126 * mm;
    G4double Rmax_sph2 = 150 * mm;

    auto solidSphere2 = new G4Sphere("AirDetector", Rmin_sph2, Rmax_sph2, SPhi_sph1, DPhi_sph1, STheta_sph1, DTheta_sph1);

    auto logicSphere2 = new G4LogicalVolume(solidSphere2,  // its solid
                                            world_mat,     // its material
                                            "AirSphere2"); // its name

    new G4PVPlacement(rot4,           // rotation
                      pos4,           // at position
                      logicSphere2,   // its logical volume
                      "AirSphere2",   // its name
                      logicWorld,     // its mother  volume
                      false,          // no boolean operation
                      0,              // copy number
                      checkOverlaps); // overlaps checking

    logicSphere2->SetVisAttributes(sphericalVisAttributes);

    //
    // Sphere detector volume 3
    //

    G4double Rmin_sph3 = 151 * mm;
    G4double Rmax_sph3 = 175 * mm;

    auto solidSphere3 = new G4Sphere("AirDetector", Rmin_sph3, Rmax_sph3, SPhi_sph1, DPhi_sph1, STheta_sph1, DTheta_sph1);

    auto logicSphere3 = new G4LogicalVolume(solidSphere3,  // its solid
                                            world_mat,     // its material
                                            "AirSphere3"); // its name

    new G4PVPlacement(rot4,           // rotation
                      pos4,           // at position
                      logicSphere3,   // its logical volume
                      "AirSphere3",   // its name
                      logicWorld,     // its mother  volume
                      false,          // no boolean operation
                      0,              // copy number
                      checkOverlaps); // overlaps checking

    logicSphere3->SetVisAttributes(sphericalVisAttributes);

    //
    // Sphere detector volume 4
    //

    G4double Rmin_sph4 = 176 * mm;
    G4double Rmax_sph4 = 200 * mm;

    auto solidSphere4 = new G4Sphere("AirDetector", Rmin_sph4, Rmax_sph4, SPhi_sph1, DPhi_sph1, STheta_sph1, DTheta_sph1);

    auto logicSphere4 = new G4LogicalVolume(solidSphere4,  // its solid
                                            world_mat,     // its material
                                            "AirSphere4"); // its name

    new G4PVPlacement(rot4,           // rotation
                      pos4,           // at position
                      logicSphere4,   // its logical volume
                      "AirSphere4",   // its name
                      logicWorld,     // its mother  volume
                      false,          // no boolean operation
                      0,              // copy number
                      checkOverlaps); // overlaps checking

    logicSphere4->SetVisAttributes(sphericalVisAttributes);

    //
    // Sphere detector volume 5
    //

    G4double Rmin_sph5 = 201 * mm;
    G4double Rmax_sph5 = 225 * mm;

    auto solidSphere5 = new G4Sphere("AirDetector", Rmin_sph5, Rmax_sph5, SPhi_sph1, DPhi_sph1, STheta_sph1, DTheta_sph1);

    auto logicSphere5 = new G4LogicalVolume(solidSphere5,  // its solid
                                            world_mat,     // its material
                                            "AirSphere5"); // its name

    new G4PVPlacement(rot4,           // rotation
                      pos4,           // at position
                      logicSphere5,   // its logical volume
                      "AirSphere5",   // its name
                      logicWorld,     // its mother  volume
                      false,          // no boolean operation
                      0,              // copy number
                      checkOverlaps); // overlaps checking

    logicSphere5->SetVisAttributes(sphericalVisAttributes);

    //
    // Sphere detector volume 6
    //

    G4double Rmin_sph6 = 226 * mm;
    G4double Rmax_sph6 = 250 * mm;

    auto solidSphere6 = new G4Sphere("AirDetector", Rmin_sph6, Rmax_sph6, SPhi_sph1, DPhi_sph1, STheta_sph1, DTheta_sph1);

    auto logicSphere6 = new G4LogicalVolume(solidSphere6,  // its solid
                                            world_mat,     // its material
                                            "AirSphere5"); // its name

    new G4PVPlacement(rot4,           // rotation
                      pos4,           // at position
                      logicSphere6,   // its logical volume
                      "AirSphere6",   // its name
                      logicWorld,     // its mother  volume
                      false,          // no boolean operation
                      0,              // copy number
                      checkOverlaps); // overlaps checking

    logicSphere6->SetVisAttributes(sphericalVisAttributes);

    // Set AirDetector as scoring volume
    //
    fScoringVolume1 = logicSphere1;
    fScoringVolume2 = logicSphere2;
    fScoringVolume3 = logicSphere3;
    fScoringVolume4 = logicSphere4;
    fScoringVolume5 = logicSphere5;
    fScoringVolume6 = logicSphere6;

    //
    // always return the physical World
    //
    return physWorld;
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
