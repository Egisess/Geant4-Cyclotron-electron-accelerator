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
/// \file B1_300/src/SteppingAction.cc
/// \brief Implementation of the B1_300::SteppingAction class

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
#include <algorithm>

extern G4String output_name;

namespace B1_300
{

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  SteppingAction::SteppingAction(EventAction *eventAction)
      : fEventAction(eventAction), fVolumesInitialized(false)
  {
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void SteppingAction::UserSteppingAction(const G4Step *step)
  {
    // Инициализация детекторов при первом вызове
    if (!fVolumesInitialized)
    {
      const auto detConstruction = static_cast<const DetectorConstruction *>(
          G4RunManager::GetRunManager()->GetUserDetectorConstruction());

      fScoringVolumes = detConstruction->GetScoringVolumes();
      fWaterDetectors = detConstruction->GetWaterDetectors();

      for (size_t i = 0; i < fWaterDetectors.size(); ++i)
      {
        fWaterDetectorMap[fWaterDetectors[i]] = i;
      }

      G4cout << "SteppingAction initialized with " << fWaterDetectors.size()
             << " water detectors" << G4endl;

      fVolumesInitialized = true;
    }

    // Получаем текущий объем и физический объем
    G4LogicalVolume *volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
    G4VPhysicalVolume *physVol = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

    // ============================================================
    // ВОДЯНОЙ ДЕТЕКТОР
    // Проверяем ДО фильтра по edepStep — нам нужны все шаги
    // гамма-кванта внутри ячейки, чтобы суммировать дозу.
    // ============================================================
    auto it = fWaterDetectorMap.find(volume);
    if (it != fWaterDetectorMap.end())
    {
      G4int trackID = step->GetTrack()->GetTrackID();
      G4String particleName = step->GetTrack()->GetParticleDefinition()->GetParticleName();
      G4int copyNo = physVol->GetCopyNo();

      G4ThreeVector position = step->GetPreStepPoint()->GetPosition();
      G4double x = position.x();
      G4double y = position.y();
      G4double z = position.z();

      G4double kinEn = step->GetPreStepPoint()->GetKineticEnergy();
      G4double edepStep = step->GetTotalEnergyDeposit();
      G4int isFirst = step->IsFirstStepInVolume() ? 1 : 0;

      // Индексы ячейки в сетке детектора
      const auto detConstruction = static_cast<const DetectorConstruction *>(
          G4RunManager::GetRunManager()->GetUserDetectorConstruction());
      G4int ix, iy, iz;
      detConstruction->GetIndices(copyNo, ix, iy, iz);

      // Записываем в файл только гамма-кванты
      if (particleName == "gamma")
      {
        FILE *f_out = fopen(output_name, "a");
        if (f_out)
        {
          fprintf(f_out, "%d %d %d %d %d %e %e %e %e %s %e %d\n",
                  trackID,                                 // Номер трека
                  copyNo,                                  // Номер ячейки (уникальный)
                  ix, iy, iz,                              // Индексы ячейки в решётке
                  x / mm, y / mm, z / mm,                  // Координаты в мм
                  kinEn / keV,                             // Кинетическая энергия в кэВ
                  static_cast<char const *>(particleName), // Тип частицы
                  edepStep / keV,                          // Отложенная энергия в кэВ
                  isFirst);                                // 1 = первый шаг в ячейке
          fclose(f_out);
        }
      }

      return; // Выходим — водяной детектор обработан
    }

    // ============================================================
    // СТАРЫЕ ДЕТЕКТОРЫ (сферы)
    // Для них фильтр по edepStep остаётся на своём месте.
    // ============================================================
    G4double edepStep = step->GetTotalEnergyDeposit();
    if (edepStep == 0.)
      return;

    fEventAction->AddEdep(edepStep);

    int detectorNumber = -1;
    for (size_t i = 0; i < fScoringVolumes.size(); ++i)
    {
      if (volume == fScoringVolumes[i])
      {
        detectorNumber = i + 1;
        break;
      }
    }

    if (detectorNumber == -1)
      return;

    if (!step->IsFirstStepInVolume())
      return;

    G4int trackID = step->GetTrack()->GetTrackID();
    G4String particleName = step->GetTrack()->GetParticleDefinition()->GetParticleName();

    G4ThreeVector position = step->GetPreStepPoint()->GetPosition();
    G4double x = position.x();
    G4double y = position.y();
    G4double z = position.z();

    G4double kinEn = step->GetPreStepPoint()->GetKineticEnergy();

    FILE *f_out = fopen(output_name, "a");
    if (f_out)
    {
      fprintf(f_out, "SPHERE %d %d %e %e %e %e %s %e\n",
              trackID,
              detectorNumber,
              x / mm, y / mm, z / mm,
              kinEn / keV,
              static_cast<char const *>(particleName),
              edepStep / keV);
      fclose(f_out);
    }
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}