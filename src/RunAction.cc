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
/// \file RunAction.cc
/// \brief Implementation of the RunAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "Run.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* prim)
: G4UserRunAction(),
fDetector(det), fPrimary(prim), fRun(0), fHistoManager(0)
{
    // Book predefined histograms
    fHistoManager = new HistoManager();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
    delete fHistoManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* RunAction::GenerateRun()
{ 
    fRun = new Run(fDetector);
    return fRun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
    int fdebug = 0 ;
    if (fdebug>1) std::cout << "RunAction::BeginOfRunAction(const G4Run*)" << std::endl;
    // save Rndm status
    G4RunManager::GetRunManager()->SetRandomNumberStore(false);
    if (isMaster) G4Random::showEngineStatus();
    
    // keep run condition
    if (fPrimary) {
        G4ParticleDefinition* particle
        = fPrimary->GetParticleGun()->GetParticleDefinition();
        G4double energy = fPrimary->GetParticleGun()->GetParticleEnergy();
        fRun->SetPrimary(particle, energy);
    }
    
    //    //histograms
    //    //
    //    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    //    if ( analysisManager->IsActive() ) {
    //        analysisManager->OpenFile();
    //    }
    
    // output csv file
    // copy from EventAction::EndOfEventAction(const G4Event*evt)
    csvfile.open("particles.csv", std::ios_base::out);
    
    csvfile
    << "eventId"          << ","
    << "trackId"          << ","
    << "parentId"         << ","
    << "PrtclName"        << ","
    << "|p_init|/MeV"       << ","
    << "CreatorProcessName" << ",";

    // track hit-position and energy deposition in scintillators
    for (int iVol=0; iVol<4; iVol++){
        csvfile
        << "EdepVol"    <<iVol<<"/MeV"      << ",";
    }
     
    csvfile
    << "PDGcode"            << ","
    << "p_init_x/MeV"       << ","
    << "p_init_y/MeV"       << ","
    << "p_init_z/MeV"       << ",";
    
    for (int iVol=0; iVol<4; iVol++){
        csvfile
        << "FirstPtVol" <<iVol<<"_Ek/MeV"   << ","
        << "FirstPtVol" <<iVol<<"_time/ns"  << ","
        << "FirstPtVol" <<iVol<<"_x/mm"     << ","
        << "FirstPtVol" <<iVol<<"_y/mm"     << ","
        << "FirstPtVol" <<iVol<<"_z/mm"     << ","
        << "LastPtVol"  <<iVol<<"_x/mm"     << ","
        << "LastPtVol"  <<iVol<<"_y/mm"     << ","
        << "LastPtVol"  <<iVol<<"_z/mm"     << ",";
    }
    
    // end line
    csvfile << std::endl;
    
    // close file
    csvfile.close();
    
    if (fdebug>1) std::cout << "done RunAction::BeginOfRunAction(const G4Run*)" << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
    int fdebug = 0;
    if (fdebug>1) std::cout << "RunAction::EndOfRunAction(const G4Run*)" << std::endl;
    if (isMaster) fRun->EndOfRun();
    
    //    //save histograms
    //    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    //    if ( analysisManager->IsActive() ) {
    //        analysisManager->Write();
    //        analysisManager->CloseFile();
    //    }
    
    // show Rndm status
    if (isMaster) G4Random::showEngineStatus();
    if (fdebug>1) std::cout << "done RunAction::EndOfRunAction(const G4Run*)" << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
