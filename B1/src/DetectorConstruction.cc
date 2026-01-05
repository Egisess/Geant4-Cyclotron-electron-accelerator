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
// * regarding  this software system or assume any liability for its  *
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
    G4Material *tube_mat = nist->FindOrBuildMaterial("G4_Cu");

    // Option to switch on/off checking of volumes overlaps
    G4bool checkOverlaps = true;

    //
    // World
    //
    G4double world_sizeXY = 400 * cm;
    G4double world_sizeZ = 400 * cm;
    G4Material *world_mat = nist->FindOrBuildMaterial("G4_AIR");

    auto solidWorld = new G4Box("World",
                                0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ);

    auto logicWorld = new G4LogicalVolume(solidWorld,
                                          world_mat,
                                          "World");

    auto physWorld = new G4PVPlacement(nullptr,
                                       G4ThreeVector(),
                                       logicWorld,
                                       "World",
                                       nullptr,
                                       false,
                                       0,
                                       checkOverlaps);

    G4double L_tube = 50 * mm + 30 * mm;
    G4VisAttributes *copperVisAttributes = new G4VisAttributes(G4Colour(0.7, 0.4, 0.1));

    //
    // Target 1 layer
    //
    G4Material *target1_mat = nist->FindOrBuildMaterial("G4_W");
    G4ThreeVector pos1 = G4ThreeVector(0, 0, 50.0 * mm);

    // Blinchik wolfram
    G4double R_inner_target = 0 * mm;
    G4double R_outer_target = 60 * mm;
    G4double L_tube_target = 0.5 * mm;
    G4double phi_0_target = 0;
    G4double phi_1_target = 2 * M_PI;
    G4ThreeVector n_bot_target = G4ThreeVector(0., -1, -std::sqrt(3));
    G4ThreeVector n_top_target = G4ThreeVector(0., 1, std::sqrt(3));

    auto solidTarget1 = new G4CutTubs("WolfTarget", R_inner_target, R_outer_target, L_tube_target, phi_0_target, phi_1_target, n_bot_target, n_top_target);
    auto logicTarget1 = new G4LogicalVolume(solidTarget1,
                                            target1_mat,
                                            "WolfTarget");

    new G4PVPlacement(nullptr,
                      pos1,
                      logicTarget1,
                      "WolfTarget",
                      logicWorld,
                      false,
                      0,
                      checkOverlaps);

    G4VisAttributes *wolframVisAttributes = new G4VisAttributes(G4Colour(0.75, 0.75, 0.75));
    logicTarget1->SetVisAttributes(wolframVisAttributes);

    // //
    // // Sphere detector volume 1
    // //
    // G4ThreeVector pos4 = pos1;
    // G4RotationMatrix *rot4 = new G4RotationMatrix();
    // rot4->rotateX(0 * deg);
    // rot4->rotateY(0 * deg);
    // rot4->rotateZ(0 * deg);

    // G4double Rmin_sph1 = 100 * mm;
    // G4double Rmax_sph1 = 125 * mm;
    // G4double SPhi_sph1 = 0 * deg;
    // G4double DPhi_sph1 = 360 * deg;
    // G4double STheta_sph1 = 0 * deg;
    // G4double DTheta_sph1 = 360 * deg;

    // auto solidSphere1 = new G4Sphere("AirDetector", Rmin_sph1, Rmax_sph1, SPhi_sph1, DPhi_sph1, STheta_sph1, DTheta_sph1);
    // auto logicSphere1 = new G4LogicalVolume(solidSphere1,
    //                                         world_mat,
    //                                         "AirSphere1");

    // new G4PVPlacement(rot4,
    //                   pos4,
    //                   logicSphere1,
    //                   "AirSphere1",
    //                   logicWorld,
    //                   false,
    //                   0,
    //                   checkOverlaps);

    // G4VisAttributes *sphericalVisAttributes = new G4VisAttributes(G4Colour(0.25, 0.25, 0.25));
    // sphericalVisAttributes->SetForceWireframe(true);
    // logicSphere1->SetVisAttributes(sphericalVisAttributes);

    // //
    // // Sphere detector volume 2
    // //
    // G4double Rmin_sph2 = 126 * mm;
    // G4double Rmax_sph2 = 150 * mm;
    // auto solidSphere2 = new G4Sphere("AirDetector", Rmin_sph2, Rmax_sph2, SPhi_sph1, DPhi_sph1, STheta_sph1, DTheta_sph1);
    // auto logicSphere2 = new G4LogicalVolume(solidSphere2,
    //                                         world_mat,
    //                                         "AirSphere2");

    // new G4PVPlacement(rot4,
    //                   pos4,
    //                   logicSphere2,
    //                   "AirSphere2",
    //                   logicWorld,
    //                   false,
    //                   0,
    //                   checkOverlaps);
    // logicSphere2->SetVisAttributes(sphericalVisAttributes);

    // //
    // // Sphere detector volume 3
    // //
    // G4double Rmin_sph3 = 151 * mm;
    // G4double Rmax_sph3 = 175 * mm;
    // auto solidSphere3 = new G4Sphere("AirDetector", Rmin_sph3, Rmax_sph3, SPhi_sph1, DPhi_sph1, STheta_sph1, DTheta_sph1);
    // auto logicSphere3 = new G4LogicalVolume(solidSphere3,
    //                                         world_mat,
    //                                         "AirSphere3");

    // new G4PVPlacement(rot4,
    //                   pos4,
    //                   logicSphere3,
    //                   "AirSphere3",
    //                   logicWorld,
    //                   false,
    //                   0,
    //                   checkOverlaps);
    // logicSphere3->SetVisAttributes(sphericalVisAttributes);

    // //
    // // Sphere detector volume 4
    // //
    // G4double Rmin_sph4 = 176 * mm;
    // G4double Rmax_sph4 = 200 * mm;
    // auto solidSphere4 = new G4Sphere("AirDetector", Rmin_sph4, Rmax_sph4, SPhi_sph1, DPhi_sph1, STheta_sph1, DTheta_sph1);
    // auto logicSphere4 = new G4LogicalVolume(solidSphere4,
    //                                         world_mat,
    //                                         "AirSphere4");

    // new G4PVPlacement(rot4,
    //                   pos4,
    //                   logicSphere4,
    //                   "AirSphere4",
    //                   logicWorld,
    //                   false,
    //                   0,
    //                   checkOverlaps);
    // logicSphere4->SetVisAttributes(sphericalVisAttributes);

    // //
    // // Sphere detector volume 5
    // //
    // G4double Rmin_sph5 = 201 * mm;
    // G4double Rmax_sph5 = 225 * mm;
    // auto solidSphere5 = new G4Sphere("AirDetector", Rmin_sph5, Rmax_sph5, SPhi_sph1, DPhi_sph1, STheta_sph1, DTheta_sph1);
    // auto logicSphere5 = new G4LogicalVolume(solidSphere5,
    //                                         world_mat,
    //                                         "AirSphere5");

    // new G4PVPlacement(rot4,
    //                   pos4,
    //                   logicSphere5,
    //                   "AirSphere5",
    //                   logicWorld,
    //                   false,
    //                   0,
    //                   checkOverlaps);
    // logicSphere5->SetVisAttributes(sphericalVisAttributes);

    // //
    // // Sphere detector volume 6
    // //
    // G4double Rmin_sph6 = 226 * mm;
    // G4double Rmax_sph6 = 250 * mm;
    // auto solidSphere6 = new G4Sphere("AirDetector", Rmin_sph6, Rmax_sph6, SPhi_sph1, DPhi_sph1, STheta_sph1, DTheta_sph1);
    // auto logicSphere6 = new G4LogicalVolume(solidSphere6,
    //                                         world_mat,
    //                                         "AirSphere6");

    // new G4PVPlacement(rot4,
    //                   pos4,
    //                   logicSphere6,
    //                   "AirSphere6",
    //                   logicWorld,
    //                   false,
    //                   0,
    //                   checkOverlaps);
    // logicSphere6->// Radiation capacitors - Water detector grid

    // Radiation capacitors - Water detector grid

    G4Material *water_material = nist->FindOrBuildMaterial("G4_WATER");

    double start_x = -500 * mm;
    double start_z = -500 * mm;
    double start_y = -100 * mm;

    double end_x = 500 * mm;
    double end_z = 500 * mm;
    double end_y = -600 * mm;

    double dx = 50 * mm;
    double dz = 50 * mm;
    double dy = 20 * mm;

    // Calculate number of cells in each direction
    fNX = (int)((end_x - start_x) / dx);
    fNZ = (int)((end_z - start_z) / dz);
    fNY = (int)((start_y - end_y) / dy);

    G4cout << "Creating water detector grid: " << fNX << "x" << fNZ << "x" << fNY << " cells" << G4endl;

    // Reserve space for water detectors
    fWaterDetectors.reserve(fNX * fNY * fNZ);

    // Create box shape for water detector cell
    auto solidWaterCell = new G4Box("WaterCell",
                                    0.5 * dx,
                                    0.5 * dy,
                                    0.5 * dz);

    // Visualization attributes for water
    G4VisAttributes *waterVisAttributes = new G4VisAttributes(G4Colour(0.0, 0.5, 1.0, 0.3));
    waterVisAttributes->SetForceSolid(true);

    // Create grid of water detectors
    int cellCounter = 0;
    for (int ix = 0; ix < fNX; ix++)
    {
      for (int iz = 0; iz < fNZ; iz++)
      {
        for (int iy = 0; iy < fNY; iy++)
        {
          // Calculate position of cell center
          G4double x_pos = start_x + (ix + 0.5) * dx;
          G4double z_pos = start_z + (iz + 0.5) * dz;
          G4double y_pos = start_y - (iy + 0.5) * dy;

          G4ThreeVector cellPosition = G4ThreeVector(x_pos, y_pos, z_pos);

          // Create logical volume for this cell
          G4String cellName = "WaterCell_" + std::to_string(ix) + "_" +
                              std::to_string(iy) + "_" + std::to_string(iz);

          auto logicWaterCell = new G4LogicalVolume(solidWaterCell,
                                                    water_material,
                                                    cellName);

          logicWaterCell->SetVisAttributes(waterVisAttributes);

          // Add to vector of water detectors
          fWaterDetectors.push_back(logicWaterCell);

          // Place the cell
          new G4PVPlacement(nullptr,
                            cellPosition,
                            logicWaterCell,
                            cellName,
                            logicWorld,
                            false,
                            cellCounter,
                            checkOverlaps);

          cellCounter++;
        }
      }
    }

    G4cout << "Total water detector cells created: " << cellCounter << G4endl;

    // Set scoring volumes
    // fScoringVolume1 = logicSphere1;
    // fScoringVolume2 = logicSphere2;
    // fScoringVolume3 = logicSphere3;
    // fScoringVolume4 = logicSphere4;
    // fScoringVolume5 = logicSphere5;
    // fScoringVolume6 = logicSphere6;

    return physWorld;
  }
  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
}