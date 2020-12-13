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
/// \file rdecay02.cc
/// \brief Main program of the radioactivedecay/rdecay02 example
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Types.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "QBBC.hh"
#include "Randomize.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"
#include "SteppingVerbose.hh"

#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"
#include <time.h>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) {
    int fdebug = 0 ;
    //detect interactive mode (if no arguments) and define UI session
    G4UIExecutive* ui = 0;
    if (argc == 1) ui = new G4UIExecutive(argc,argv);
    
    if (fdebug>1) std::cout << "main(): choose the Random engine" << std::endl;
    //choose the Random engine
    // set a new seed every time
    time_t timer;
    struct tm y2k = {0};
    y2k.tm_hour = 0;   y2k.tm_min = 0; y2k.tm_sec = 0;
    y2k.tm_year = 100; y2k.tm_mon = 0; y2k.tm_mday = 1;
    timer = time(NULL);
    int seconds = int(difftime(timer,mktime(&y2k)));
    G4Random::setTheEngine(new CLHEP::RanecuEngine(seconds));
    
    if (fdebug>1) std::cout << "main(): Construct the default run manager" << std::endl;
    // Construct the default run manager
#ifdef G4MULTITHREADED
    G4MTRunManager* runManager = new G4MTRunManager;
    G4int nThreads = G4Threading::G4GetNumberOfCores();
    if (argc==3) nThreads = G4UIcommand::ConvertToInt(argv[2]);
    runManager->SetNumberOfThreads(nThreads);
#else
    if (fdebug>1) std::cout << "G4VSteppingVerbose::SetInstance(new SteppingVerbose);" << std::endl;
    //my Verbose output class
    G4VSteppingVerbose::SetInstance(new SteppingVerbose);
    if (fdebug>1) std::cout << "before G4RunManager* runManager = new G4RunManager;" << std::endl;
    G4RunManager* runManager = new G4RunManager;
    if (fdebug>1) std::cout << "done G4RunManager* runManager = new G4RunManager;" << std::endl;
#endif
    
    if (fdebug>1) std::cout << "main(): set mandatory initialization classes" << std::endl;
    //set mandatory initialization classes
    DetectorConstruction* det= new DetectorConstruction;
    runManager->SetUserInitialization(det);
    
    if (fdebug>1) std::cout << "main(): construct the physics list" << std::endl;
    // construct the physics list       // from example B1
    G4VModularPhysicsList* physicsList = new QBBC;
    if (fdebug>1) std::cout << "G4VModularPhysicsList* physicsList = new QBBC;" << std::endl;
    physicsList->SetVerboseLevel(1);
    runManager->SetUserInitialization(physicsList);
    
    if (fdebug>1) std::cout << "main(): runManager->SetUserInitialization(physicsList);" << std::endl;
    // from rdecay01
    //    runManager->SetUserInitialization(new PhysicsList);
    //    // original
    PhysicsList* phys = new PhysicsList;
    runManager->SetUserInitialization(phys);
    
    if (fdebug>1) std::cout << "main(): runManager->SetUserInitialization(new ActionInitialization(det));" << std::endl;
    runManager->SetUserInitialization(new ActionInitialization(det));
    
    //initialize visualization
    //  G4VisManager* visManager = nullptr;
    G4VisManager* visManager = new G4VisExecutive;
    // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
    // G4VisManager* visManager = new G4VisExecutive("Quiet");
    visManager->Initialize();
    
    // get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    
    if (fdebug>1) std::cout << "UImanager" << std::endl;
    if (ui)  {
        //   //interactive mode
        //   visManager = new G4VisExecutive;
        //   visManager->Initialize();
        // interactive mode
        if (fdebug>1) std::cout << "UImanager->ApplyCommand(/control/execute erez_vis.mac);" << std::endl;
        //        UImanager->ApplyCommand("/control/execute vis.mac");
                UImanager->ApplyCommand("/control/execute erez_vis.mac");
//                UImanager->ApplyCommand("/control/execute no_vis.mac"); THIS DOESNT WORK FOR SOME REASON
        if (fdebug>1) std::cout << "ui->SessionStart();" << std::endl;
        ui->SessionStart();
        if (fdebug>1) std::cout << "delete ui;" << std::endl;
        delete ui;
    }
    else  {
        //batch mode
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command+fileName);
    }
    
    //job termination
    delete visManager;
    delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
