//
// ********************************************************************
// * License and Disclaimer *
// 
// * The Geant4 software is copyright of the Copyright Holders of *
// * the Geant4 Collaboration. It is provided under the terms and *
// * conditions of the Geant4 Software License, included in the file *
// * LICENSE and available at http://cern.ch/geant4/license . These *
// * include a list of copyright holders. *
// 
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work make any representation or warranty, express or implied, *
// * regarding this software system or assume any liability for its *
// * use. Please see the license in the file LICENSE and URL above *
// * for the full disclaimer and the limitation of liability. *
// 
// * This code implementation is the result of the scientific and *
// * technical work of the GEANT4 collaboration. *
// * By using, copying, modifying or distributing the software (or *
// * any work based on the software) you agree to acknowledge its *
// * use in resulting scientific publications, and indicate your *
// * acceptance of all terms of the Geant4 Software license. *
// ********************************************************************
//
//
/// \file B1/src/SteppingAction.cc
/// \brief Implementation of the B1::SteppingAction class
#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"
#include "stdio.h"
#include "stdlib.h"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4VProcess.hh"
#include "G4SteppingManager.hh"
#include "G4SystemOfUnits.hh"
#include <vector>

extern G4String output_name;
G4double x, y, z, kinEn;
G4int track_num;
G4String particleName;
G4bool isFirstStepInVolume;

// Track the total energy deposited in each volume
static G4double totalEdepMain = 0.0;
static G4double totalEdepPan07 = 0.0;
static G4double totalEdepPan10 = 0.0;

namespace B1
{
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
SteppingAction::SteppingAction(EventAction* eventAction)
: fEventAction(eventAction)
{}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void SteppingAction::UserSteppingAction(const G4Step* step)
{
  // Initialize the scoring volumes if they haven't been already
  if (!fScoringVolume || !fScoringPan07 || !fScoringPan10) {
    const auto detConstruction = static_cast<const DetectorConstruction*>(
      G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    
    // Get all scoring volumes from the detector construction
    std::vector<G4LogicalVolume*> scoringVolumes = detConstruction->GetScoringVolumes();
    
    // Assign the volumes to the appropriate pointers
    // Assuming the volumes are returned in the same order they are defined in DetectorConstruction
    if (scoringVolumes.size() > 0) fScoringVolume = scoringVolumes[0];
    if (scoringVolumes.size() > 1) fScoringPan07 = scoringVolumes[1];
    if (scoringVolumes.size() > 2) fScoringPan10 = scoringVolumes[2];
  }
  
  // Get volume of the current step
  G4LogicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()
                              ->GetVolume()->GetLogicalVolume();
  
  // Calculate energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit();
  
  // Process the step if it's in any of our scoring volumes
  bool isInScoringVolume = false;
  G4String volumeName;
  G4double currentTotalEdep = 0.0;  // Total energy deposited in the current volume
  
  if (volume == fScoringVolume) {
    isInScoringVolume = true;
    volumeName = "ScoringVolume";
    totalEdepMain += edepStep;
    currentTotalEdep = totalEdepMain;
    
    // Add energy deposition for the main scoring volume
    fEventAction->AddEdep(edepStep);
  }
  else if (volume == fScoringPan07) {
    isInScoringVolume = true;
    volumeName = "ScoringPan07";
    totalEdepPan07 += edepStep;
    currentTotalEdep = totalEdepPan07;
    
    // If you have a method to track energy in Pan07
    // fEventAction->AddEdepPan07(edepStep);
  }
  else if (volume == fScoringPan10) {
    isInScoringVolume = true;
    volumeName = "ScoringPan10";
    totalEdepPan10 += edepStep;
    currentTotalEdep = totalEdepPan10;
    
    // If you have a method to track energy in Pan10
    // fEventAction->AddEdepPan10(edepStep);
  }
  
  if (!isInScoringVolume) return;
  
  // Track information
  track_num = step->GetTrack()->GetTrackID();
  x = step->GetPreStepPoint()->GetPosition()[0];
  y = step->GetPreStepPoint()->GetPosition()[1];
  z = step->GetPreStepPoint()->GetPosition()[2];
  kinEn = step->GetPreStepPoint()->GetKineticEnergy();
  isFirstStepInVolume = step->IsFirstStepInVolume();
  particleName = step->GetTrack()->GetParticleDefinition()->GetParticleName();
  
  // Calculate dose in Gray (J/kg)
  G4double dose = 0.0;
  if (volume == fScoringVolume && fScoringVolume) {
    G4double mass = fScoringVolume->GetMass();
    if (mass > 0) {
      dose = currentTotalEdep / mass / gray;  // Convert to Gray units
    }
  }
  else if (volume == fScoringPan07 && fScoringPan07) {
    G4double mass = fScoringPan07->GetMass();
    if (mass > 0) {
      dose = currentTotalEdep / mass / gray;  // Convert to Gray units
    }
  }
  else if (volume == fScoringPan10 && fScoringPan10) {
    G4double mass = fScoringPan10->GetMass();
    if (mass > 0) {
      dose = currentTotalEdep / mass / gray;  // Convert to Gray units
    }
  }
  
  FILE *f_out;
  f_out = fopen(output_name, "a");
  fprintf(f_out, "%d %e %e %e %e %d %s %s %e\n", 
          track_num, x / mm, y / mm, z / mm, kinEn / keV, 
          isFirstStepInVolume, static_cast<char const *>(particleName),
          static_cast<char const *>(volumeName), dose);
  fclose(f_out);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
}