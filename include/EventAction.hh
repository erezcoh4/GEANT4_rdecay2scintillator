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
/// \file EventAction.hh
/// \brief Definition of the EventAction class
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1
#define NMAXtracks 1000
#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include <fstream>
#include <vector>
#include <string>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
  public:
    EventAction();
   ~EventAction();

  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void   EndOfEventAction(const G4Event*);
    
    void AddEdep (G4int iVol,
                  G4double edep,
                  G4double time,
                  G4int trackId,
                  int fdebug = 0);
           
    
    void SetFirstPointInVolume( G4int iVol, G4int trackId,
                               G4ThreeVector pos,
                               G4double time,
                               G4double Ek,
                               G4String fProcessName){
        int fdebug=0;
        if (fdebug>1){
            std::cout
            << "pos: (" << pos.x() << "," << pos.y() << ","<< pos.z() << ")"
            << ", time:" << time
            << "ns, Ek:" << Ek <<  "MeV, fProcessName:" << fProcessName
            << std::endl;
        }
        
        FirstPointInVolume[iVol][trackId]       = pos;
        FirstPointInVolumeTime[iVol][trackId]   = time;
        FirstPointInVolumeEk[iVol][trackId]     = Ek;
        
        ProcessName.at(trackId) = fProcessName;
        
    };
    void SetLastPointInVolume( G4int iVol, G4int trackId, G4ThreeVector pos, G4double time ){
        LastPointInVolume[iVol][trackId] = pos;
        LastPointInVolumeTime[iVol][trackId] = time;
    };
    
    std::ofstream csvfile;
    
  private:
    G4double fEdep1,   fEdep2;
    G4double fWeight1, fWeight2;
    G4double fTime0;
    
    // fEdep is a 2-D array
    // first dimension is volume number
    // (0 = source holder, 1=scintllator 1, 2 = scintillator 2, 3 = world)
    // second dimension is track Id
    G4double fEdep[4][NMAXtracks];
    
    // first and last point of each track in each volume
    G4ThreeVector FirstPointInVolume[4][NMAXtracks];
    G4ThreeVector LastPointInVolume[4][NMAXtracks];
    G4double FirstPointInVolumeTime[4][NMAXtracks];
    G4double LastPointInVolumeTime[4][NMAXtracks];
    G4double FirstPointInVolumeEk[4][NMAXtracks];
    std::vector<G4String> ProcessName;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
