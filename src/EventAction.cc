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
/// \file EventAction.cc
/// \brief Implementation of the EventAction class
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"

#include "Run.hh"
#include "HistoManager.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
:G4UserEventAction(),
fEdep1(0.), fEdep2(0.), fWeight1(0.), fWeight2(0.),
fTime0(-1*s)
{ } 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
    int fdebug = 0;
    if (fdebug>1) std::cout << "EventAction::BeginOfEventAction(const G4Event*)" << std::endl;
    fEdep1 = fEdep2 = fWeight1 = fWeight2 = 0.;
    fTime0 = -1*s;
    if (fdebug>1) std::cout << "done EventAction::BeginOfEventAction(const G4Event*)" << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::AddEdep(G4int iVol, G4double edep,
                          G4double time, G4double weight)
{
    int fdebug = 0;
    if (fdebug>1) std::cout << "EventAction::AddEdep()" << std::endl;
    // initialize t0
    if (fTime0 < 0.) fTime0 = time;
    
    // out of time window ?
    const G4double TimeWindow (1*microsecond);
    if (std::fabs(time - fTime0) > TimeWindow) return;
    
    if (iVol == 1) { fEdep1 += edep; fWeight1 += edep*weight;}
    if (iVol == 2) { fEdep2 += edep; fWeight2 += edep*weight;}
    if (fdebug>1) std::cout << "done EventAction::AddEdep()" << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event*evt)
{
    int fdebug = 2;
    if (fdebug>0) std::cout << "EventAction::EndOfEventAction(const G4Event*evt)" << std::endl;
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    
    G4double Etot = fEdep1 + fEdep2;
    G4double Wtot = (fWeight1 + fWeight2)/Etot;
    // pulse height in target
    if (fEdep1 > 0.) {
        fWeight1 /= fEdep1;
        analysisManager->FillH1(0, fEdep1, fWeight1);
    }
    // pulse height in detector
    if (fEdep2 > 0.) {
        fWeight2 /= fEdep2;
        analysisManager->FillH1(1, fEdep2, fWeight2);
    }
    // total
    analysisManager->FillH1(2, Etot, Wtot);
    // threshold in target and detector
    const G4double Threshold1(10*keV), Threshold2(10*keV);
    //coincidence, anti-coincidences
    G4bool coincidence       = ((fEdep1 >= Threshold1) && (fEdep2 >= Threshold2));
    G4bool anti_coincidence1 = ((fEdep1 >= Threshold1) && (fEdep2 <  Threshold2));
    G4bool anti_coincidence2 = ((fEdep1 <  Threshold1) && (fEdep2 >= Threshold2));
    if (coincidence)       analysisManager->FillH1(3, fEdep2, fWeight2);
    if (anti_coincidence1) analysisManager->FillH1(4, fEdep1, fWeight1);
    if (anti_coincidence2) analysisManager->FillH1(5, fEdep2, fWeight2);
    // pass energies to Run
    Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    run->AddEdep (fEdep1, fEdep2);
    
    
    if (fdebug>0) std::cout << "event " << evt -> GetEventID() << ",opening output csv and write data" << std::endl;
    // open output csv and write data
    csvfile.open("particles.csv", std::ios_base::app);
    // extract event data
    
    G4int eventId = evt -> GetEventID();
    // trajectory container
    G4TrajectoryContainer * trajCont = evt -> GetTrajectoryContainer();
    size_t trajCont_size = trajCont->size();

    if (fdebug>1) std::cout << "trajCont_size: " << trajCont_size << std::endl;
    G4int NtrajCont = trajCont -> entries();
    if (fdebug>1) std::cout << "NtrajCont: " << NtrajCont << std::endl;

    std::vector< G4VTrajectory * > * trajectories = trajCont -> GetVector ();
    for (auto traj:*trajectories){

        G4String ParticleName = traj->GetParticleName() ;
        G4int PDGcode = traj->GetPDGEncoding();
        G4int trackID = traj->GetTrackID ();
        G4int parentID = traj->GetParentID ();
        G4ThreeVector pInit = traj->GetInitialMomentum ();

        // for each particle, we dedicate a line in the csv file
        csvfile
        << eventId          << ","
        << NtrajCont        << ","
        << trackID          << ","
        << parentID         << ","
        << PDGcode          << ","
        << ParticleName     << ","
        << pInit.x()        << ","
        << pInit.y()        << ","
        << pInit.z()        << ","
        << pInit.mag()      << ","
        << Etot
        << std::endl;
        
    }
    
    // hits collection
    if (fdebug>1) std::cout <<  "hits collection: " << std::endl;
    
    if (fdebug>1) std::cout << "G4HCofThisEvent * hitsCol = evt->GetHCofThisEvent() for event " << eventId << std::endl;
    G4HCofThisEvent * hitsCol = evt->GetHCofThisEvent();
    G4int NhitCols = sizeof(hitsCol)/sizeof(hitsCol[0]);
    
    
    if (fdebug>1) std::cout << "we have " << NhitCols << " hit collections in event " << eventId << std::endl;
    if (NhitCols>=1){
        for (G4int hitColIdx=0; hitColIdx<NhitCols; hitColIdx++){
            G4VHitsCollection * HC = hitsCol -> GetHC (hitColIdx);
            G4int NHC = sizeof(HC)/sizeof(HC[0]);
            std::cout << NHC << " hits in hit collection " << hitColIdx << std::endl;
            
            if (NHC>0){
                G4String SDname = HC -> GetSDname ();
                size_t HCsize = HC -> GetSize ();
                if (fdebug>1) std::cout << "hit collection in " << SDname << " of size " << HCsize << std::endl;
                HC -> PrintAllHits ();
//                for (size_t hitIdx=0; hitIdx<HCsize; hitIdx++){
//                    G4VHit * hit = HC -> GetHit(hitIdx);
//                    hit -> Print();
//                }
                
                
            }
            
        }
        
        // close file
        csvfile.close();
        if (fdebug>1) std::cout << "closed file at event " << evt -> GetEventID() << std::endl;
    }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


