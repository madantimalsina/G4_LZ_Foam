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
/// \file GeSteppingAction.cc
/// \brief Implementation of the GeSteppingAction class

#include "GeSteppingAction.hh"
#include "GeEventAction.hh"
#include "GeDetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GeSteppingAction::GeSteppingAction(GeEventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringVolume(0),
  fScoringVolume1(0),
  fScoringVolume2(0),
  fScoringVolume3(0)

{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GeSteppingAction::~GeSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GeSteppingAction::UserSteppingAction(const G4Step* step)
{
  
  if (!fScoringVolume) { 
    const GeDetectorConstruction* detectorConstruction
      = static_cast<const GeDetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detectorConstruction->GetScoringVolume();   
  }
  
  
  if (!fScoringVolume1) { 
    const GeDetectorConstruction* detectorConstruction
      = static_cast<const GeDetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume1 = detectorConstruction->GetScoringVolume1();   
  }
  
  if (!fScoringVolume2) { 
    const GeDetectorConstruction* detectorConstruction
      = static_cast<const GeDetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume2 = detectorConstruction->GetScoringVolume2();   
  }

    if (!fScoringVolume3) { 
    const GeDetectorConstruction* detectorConstruction
      = static_cast<const GeDetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume3 = detectorConstruction->GetScoringVolume3();   
  }

  // get volume of the current step
  G4LogicalVolume* volume 
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();
      
  // check if we are in scoring volume
  //if (volume != fScoringVolume) return;
  
  // collect energy deposited in this step
  //G4double edepStep = step->GetTotalEnergyDeposit();
  //fEventAction->AddEdep(edepStep);
  
  
  if ((volume != fScoringVolume) && (volume != fScoringVolume1) && (volume != fScoringVolume2) && (volume != fScoringVolume3)) return;

  if ((volume != fScoringVolume1) && (volume != fScoringVolume2) && (volume != fScoringVolume3)){
  // collect energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit();
  fEventAction->AddEdep(edepStep);
  return;
  }
  
  if ((volume != fScoringVolume) && (volume != fScoringVolume2) && (volume != fScoringVolume3)){
  // collect energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit();
  fEventAction->AddEdep1(edepStep);
  return;
  }
  
 if (volume == fScoringVolume3){
    // collect energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit();
  G4Track* track = step->GetTrack();
  G4ThreeVector position = track->GetPosition();
  fEventAction->AddEdepPE(edepStep, position);

 }

  G4double edepStep = step->GetTotalEnergyDeposit();
  fEventAction->AddEdep4(edepStep);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

