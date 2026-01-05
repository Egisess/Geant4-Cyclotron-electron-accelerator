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
/// \file B1_300/include/DetectorConstruction.hh
/// \brief Definition of the B1_300::DetectorConstruction class

#ifndef B1_300DetectorConstruction_h
#define B1_300DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include <vector>

class G4VPhysicalVolume;
class G4LogicalVolume;

namespace B1_300
{

  /// Detector construction class to define materials and geometry.

  class DetectorConstruction : public G4VUserDetectorConstruction
  {
  public:
    DetectorConstruction() = default;
    ~DetectorConstruction() override = default;

    G4VPhysicalVolume *Construct() override;

    // Get all water detector logical volumes
    const std::vector<G4LogicalVolume *> &GetWaterDetectors() const
    {
      return fWaterDetectors;
    }

    // Get grid dimensions for analysis
    G4int GetNX() const { return fNX; }
    G4int GetNY() const { return fNY; }
    G4int GetNZ() const { return fNZ; }

    // Get total number of water cells
    G4int GetTotalWaterCells() const { return fWaterDetectors.size(); }

    // Convert cell indices to copyNo and vice versa
    G4int GetCopyNo(G4int ix, G4int iy, G4int iz) const
    {
      return ix * fNY * fNZ + iz * fNY + iy;
    }

    void GetIndices(G4int copyNo, G4int &ix, G4int &iy, G4int &iz) const
    {
      ix = copyNo / (fNY * fNZ);
      G4int remainder = copyNo % (fNY * fNZ);
      iz = remainder / fNY;
      iy = remainder % fNY;
    }

    // Old scoring volumes (if still needed)
    std::vector<G4LogicalVolume *> GetScoringVolumes() const
    {
      std::vector<G4LogicalVolume *> volumes;
      if (fScoringVolume1)
        volumes.push_back(fScoringVolume1);
      if (fScoringVolume2)
        volumes.push_back(fScoringVolume2);
      if (fScoringVolume3)
        volumes.push_back(fScoringVolume3);
      if (fScoringVolume4)
        volumes.push_back(fScoringVolume4);
      if (fScoringVolume5)
        volumes.push_back(fScoringVolume5);
      if (fScoringVolume6)
        volumes.push_back(fScoringVolume6);
      return volumes;
    }

  protected:
    // Water detector grid
    std::vector<G4LogicalVolume *> fWaterDetectors;
    G4int fNX = 0; // Number of cells in X
    G4int fNY = 0; // Number of cells in Y
    G4int fNZ = 0; // Number of cells in Z

    // Old scoring volumes (spheres - commented out in .cc)
    G4LogicalVolume *fScoringVolume1 = nullptr;
    G4LogicalVolume *fScoringVolume2 = nullptr;
    G4LogicalVolume *fScoringVolume3 = nullptr;
    G4LogicalVolume *fScoringVolume4 = nullptr;
    G4LogicalVolume *fScoringVolume5 = nullptr;
    G4LogicalVolume *fScoringVolume6 = nullptr;
  };
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif