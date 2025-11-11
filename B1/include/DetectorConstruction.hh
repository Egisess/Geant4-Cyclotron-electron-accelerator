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
/// \file B1/include/DetectorConstruction.hh
/// \brief Definition of the B1::DetectorConstruction class

#ifndef B1DetectorConstruction_h
#define B1DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include <vector>

class G4VPhysicalVolume;
class G4LogicalVolume;

namespace B1
{

  /// Detector construction class to define materials and geometry.

  class DetectorConstruction : public G4VUserDetectorConstruction
  {
  public:
    DetectorConstruction() = default;
    ~DetectorConstruction() override = default;

    G4VPhysicalVolume *Construct() override;

    // Change the return type to a vector of logical volumes
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
