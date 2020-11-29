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
{
    // initialize fEdep
    // fEdep is a 2-D array
    // first dimension is volume number
    // (0 = source holder, 1=scintllator 1, 2 = scintillator 2, 4 = world)
    // second dimension is track Id
    
    for (int trackId=0; trackId<NMAXtracks; trackId++){
        ProcessName.push_back("");
        
        for (int iVol=0; iVol<4; iVol++){
        
            fEdep[iVol][trackId] = 0;
            FirstPointInVolume[iVol][trackId] = G4ThreeVector(-999,-999,-999);
            LastPointInVolume[iVol][trackId] = G4ThreeVector(-999,-999,-999);
            FirstPointInVolumeTime[iVol][trackId] = -999;
            LastPointInVolumeTime[iVol][trackId] = -999;
            FirstPointInVolumeEk[iVol][trackId] = -999;
            
        }
    }
    
    
}

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
    
    
    // --------------------------------
    // write a csv file for each particle that was emitted from the source:
    //
    // event ID
    // partice type (PDG code)
    // track ID
    // parent ID
    // initial momentum [MeV/c]
    // final momentum [MeV/c]
    //
    // hit time scintillator 1 [ps] (-999 if did not hit scintillator)
    // start point in scintillator 1 [cm] ((-999,-999,-999) if did not hit scintillator)
    // end point in scintillator 1 [cm] ((-999,-999,-999) if did not hit scintillator)
    // energy deposition in scintilaltor 1 [MeV] (0 if did not hit scintillator)
    // underwent compton in scintillator 1 ?
    //
    // hit time scintillator 2 [ps] (-999 if did not hit scintillator)
    // start point in scintillator 2 [cm] ((-999,-999,-999) if did not hit scintillator)
    // end point in scintillator 2 [cm] ((-999,-999,-999) if did not hit scintillator)
    // energy deposition in scintilaltor 2 [MeV] (0 if did not hit scintillator)
    // underwent compton in scintillator 2 ?
    //
    
    if (fdebug>0) std::cout << "event " << evt -> GetEventID() << ",opening output csv and write data" << std::endl;
    // open output csv and write data
    csvfile.open("particles.csv", std::ios_base::app);
    // extract event data
    
    G4int eventId = evt -> GetEventID();
    // trajectory container
    G4TrajectoryContainer * trajCont = evt -> GetTrajectoryContainer();
    size_t trajCont_size = trajCont->size();
    
    if (fdebug>1) std::cout << "trajCont_size: " << trajCont_size << std::endl;
    std::vector< G4VTrajectory * > * trajectories = trajCont -> GetVector ();
        
    for (auto traj:*trajectories){
        
        G4String ParticleName = traj->GetParticleName() ;
        G4int PDGcode = traj->GetPDGEncoding();
        G4int trackId = traj->GetTrackID ();
        G4int parentId = traj->GetParentID ();
        G4ThreeVector p_init = traj->GetInitialMomentum ();
        
        // for data saving, we do not want to write down neutrinos and heavy isotopes that promptly decay
        if (
            ParticleName=="nu_e" ||
            ParticleName=="Na22" ||
            ParticleName=="Ne22" ||
            ParticleName=="Ne22[1274.577]"
            )
            continue;
        
        //  trajectory points can not help, as they do not provide access for the energy deposition etc.
        // for each particle, we dedicate a line in the csv file
        csvfile
        << eventId          << ","
        << trackId          << ","
        << parentId         << ","
        << ParticleName     << ","
        << p_init.mag()/MeV             << ","
        << ProcessName.at(trackId)      << ",";
        
        
        // track hit-position and energy deposition in scintillators
        for (int iVol=0; iVol<4; iVol++){
            csvfile
            << fEdep[iVol][trackId]/MeV << ",";
        }
        csvfile
        << PDGcode                      << ","
        << p_init.x()/MeV               << ","
        << p_init.y()/MeV               << ","
        << p_init.z()/MeV               << ",";

        for (int iVol=0; iVol<4; iVol++){
            csvfile
            << FirstPointInVolumeEk[iVol][trackId]/MeV << ","
            << FirstPointInVolumeTime[iVol][trackId]/ns << ","
            << FirstPointInVolume[iVol][trackId].x()/mm << ","
            << FirstPointInVolume[iVol][trackId].y()/mm << ","
            << FirstPointInVolume[iVol][trackId].z()/mm << ","
            << LastPointInVolume[iVol][trackId].x()/mm << ","
            << LastPointInVolume[iVol][trackId].y()/mm << ","
            << LastPointInVolume[iVol][trackId].z()/mm << ",";            
        }
        
        // end line
        csvfile << std::endl;
        
    }
    csvfile.close();
    
    
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// AddEdep modified by Erez, Nov-25, 2020
void EventAction::AddEdep(G4int iVol,
                          G4double edep,
                          G4double time,
                          G4int trackId,
                          int fdebug ){
    
    if (fdebug>1) std::cout << "EventAction::AddEdep()" << std::endl;
    
    // omit tracks of id > NMAXtracks
    if (trackId >= NMAXtracks){
        std::cout
        << "(trackId = " << trackId << ") >= (NMAXtracks = " << NMAXtracks << ")"
        <<
        "returning from EventAction::AddEdep() " << std::endl;
        return;
    }
    
    // fEdep is a 2-D array
    // first dimension is volume number
    // (0 = source holder, 1=scintllator 1, 2 = scintillator 2, 4 = world)
    // second dimension is track Id
    fEdep[iVol][trackId] += edep;
    
    
    if (fdebug>1) std::cout << "done EventAction::AddEdep()" << std::endl;
    
    // original
    //    // initialize t0
    //    if (fTime0 < 0.) fTime0 = time;
    //
    //    // out of time window ?
    //    const G4double TimeWindow (1*microsecond);
    //    if (std::fabs(time - fTime0) > TimeWindow) return;
    //    if (iVol == 1) { fEdep1 += edep; fWeight1 += edep*weight;}
    //    if (iVol == 2) { fEdep2 += edep; fWeight2 += edep*weight;}
}


