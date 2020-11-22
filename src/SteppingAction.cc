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
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "Run.hh"
#include "EventAction.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, EventAction* event)
: G4UserSteppingAction(), fDetector(det), fEventAction(event)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
    int fdebug = 0;
    if (fdebug>1) std::cout << "SteppingAction::UserSteppingAction(const G4Step* aStep)" << std::endl;
    Run* run = static_cast<Run*>(
                                 G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    
    if (fdebug>1) std::cout << "//which volume ?" << std::endl;
    //which volume ?
    //
    G4LogicalVolume* lVolume = aStep->GetPreStepPoint()->GetTouchableHandle()
    ->GetVolume()->GetLogicalVolume();
    G4int iVol = 0;
    if (lVolume == fDetector->GetLogicScint_1()) iVol = 1;
    if (lVolume == fDetector->GetLogicScint_2()) iVol = 2;

    if (fdebug>1) std::cout << "//count processes" << std::endl;
    // count processes
    //
    const G4StepPoint* startPoint = aStep->GetPreStepPoint();
    const G4StepPoint* endPoint = aStep->GetPostStepPoint();
    
    if (fdebug>1) std::cout << "startPoint: "      << "("
    << endPoint->GetPosition().x() << ","
    << endPoint->GetPosition().y() << ","
    << endPoint->GetPosition().z() << ")"
    << std::endl;
    
    if (fdebug>1) std::cout << "endPoint: "      << "("
    << endPoint->GetPosition().x() << ","
    << endPoint->GetPosition().y() << ","
    << endPoint->GetPosition().z() << ")"
    << std::endl;
    
    G4double stepLength = aStep -> GetStepLength();
    if (fdebug>1) std::cout << "stepLength: " << stepLength / CLHEP::nm << "nm" << std::endl;
    if (fabs(stepLength / CLHEP::nm)<1) return;
    
    const G4VProcess* process   = endPoint->GetProcessDefinedStep();
    if (fdebug>1) std::cout << "process: " << process->GetProcessName() << std::endl;
    run->CountProcesses(process, iVol);
    
    if (fdebug>1) std::cout << "//energy deposit" << std::endl;
    // energy deposit
    //
    G4double edepStep = aStep->GetTotalEnergyDeposit();
    if (edepStep <= 0.) return;
    G4double time   = aStep->GetPreStepPoint()->GetGlobalTime();
    G4double weight = aStep->GetPreStepPoint()->GetWeight();
    fEventAction->AddEdep(iVol, edepStep, time, weight);
    
    if (fdebug>1) std::cout << "//fill ntuple id = 2" << std::endl;
    //fill ntuple id = 2
    G4int id = 2;
    analysisManager->FillNtupleDColumn(id,0, edepStep);
    analysisManager->FillNtupleDColumn(id,1, time/s);
    analysisManager->FillNtupleDColumn(id,2, weight);
    analysisManager->AddNtupleRow(id);
    if (fdebug>1) std::cout << "done SteppingAction::UserSteppingAction(const G4Step* aStep)" << std::endl;
    
    // continue here:
    // stream step energy deposit into output csv:
    // for each event, sum energy deposition in each detector,
    // and record it to output file:
    //
    // event, detector, particle type, track id, parent id, Edep, weight, time, process name, start x,y,z, end x,y,z
    
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
