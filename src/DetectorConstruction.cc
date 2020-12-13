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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),
fTargetMater(0), fLogicTarget(0),
fDetectorMater(0), fLogicDetector(0),
fWorldMater(0), fPhysiWorld(0),
fDetectorMessenger(0)
{
    //    fTargetLength      = 1*cm;
    //    fTargetRadius      = 0.5*cm;
    //    fDetectorLength    = 5*cm;
    //    fDetectorThickness = 2*cm;
    //
    //    fWorldLength = std::max(fTargetLength,fDetectorLength);
    // fWorldRadius = fTargetRadius + fDetectorThickness;
    
//    fWorldLength = 10.*cm;
//    fWorldRadius = 10.*cm;
    DefineMaterials();
    
    fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
    int fdebug = 0;
    // build materials
    G4Element* N  = new G4Element("Nitrogen", "N", 7, 14.01*g/mole);
    G4Element* O  = new G4Element("Oxygen",   "O", 8, 16.00*g/mole);
    G4int ncomponents; G4double fractionmass;
    G4Material* Air20 = new G4Material("Air", 1.205*mg/cm3, ncomponents=2, kStateGas, 293.*kelvin, 1.*atmosphere);
    Air20->AddElement(N, fractionmass=0.7);
    Air20->AddElement(O, fractionmass=0.3);
    fWorldMater = Air20;
    
    

    
    
    G4NistManager* man = G4NistManager::Instance();
    // source holder material - plastic
    fPlasticMater = man->FindOrBuildMaterial("G4_POLYETHYLENE");
    
    // scintillation materials
    // scintillation materials from [http://www.sixiangguo.net/code/geant4/AppDevelop/apas06.html]
    //    fDetectorMater = man->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    //    fDetectorMater = man->FindOrBuildMaterial("G4_SODIUM_IODIDE");
    //    fDetectorMater = man->FindOrBuildMaterial("G4_STILBENE");
    
    // LYSO [http://www.iem.cfmac.csic.es/departamentos/nuclear/fnexp/master-fusion/lab_master_MCsim.pdf]
    // LYSO: d = 7.10 g/cm3 [Lu]=71.45% [Y]= 4.03 % [Si]= 6.37% [O]=18.15%
    // Lu from [https://nexus6.us.es/Tutorial_ITA-IEAv_2015/task1/task1a.html]
    G4Element* Lu = new G4Element("Lutetium", "Lu", 71., 174.97*g/mole);
    // Yttrium from [https://gitlab.physik.uni-kiel.de/geant4/geant4/-/blob/37fff30d2e782e512ccf0193e18ea3a4bbfc5037/examples/advanced/xray_fluorescence/src/XrayFluoDetectorConstruction.cc]
    G4Element* Y  = new G4Element("Yttrium"  ,"Y" , 39., 88.905*g/mole);
    // Si from [https://nexus6.us.es/Tutorial_ITA-IEAv_2015/task1/task1a.html]
    G4Element*  Si = new G4Element("Silicon", "Si", 14., 28.09*g/mole);
    // Oxygen from above world material "Air20"...
    
    G4Material* LYSO = new G4Material("LYSO", 7.1*g/cm3, ncomponents=4,kStateSolid);
    LYSO->AddElement(Lu, fractionmass=0.7145);
    LYSO->AddElement(Y,  fractionmass=0.0403);
    LYSO->AddElement(Si, fractionmass=0.0637);
    LYSO->AddElement(O,  fractionmass=0.1815);
    fDetectorMater = LYSO;
    
    
    // SiPMs
    fSiPMMater = man->FindOrBuildMaterial("G4_Si");;
    
    
    // electronics: PCB / Si
    // glass fiber reinforced
    fElectronicsMaterial = man->FindOrBuildMaterial("G4_Si");;
    
    
    if (fdebug>0){
        std::cout << "detector material: " << fDetectorMater -> GetName() << std::endl;
        std::cout << "material table:" << std::endl;
        std::cout << *(G4Material::GetMaterialTable()) << std::endl;
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
    G4bool DoSourceHolder = TRUE;
    G4bool DoSourcePlaceHolder = FALSE;
    G4bool DoScintillators = TRUE;
    G4bool DoSiPMs = TRUE;
    G4bool DoElectronicPlates = TRUE;
    
    // Cleanup old geometry
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4PhysicalVolumeStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4SolidStore::GetInstance()->Clean();
    
    // definition
    // scintillators
    G4double Scint_dx = 3.*mm;
    G4double Scint_dy = 3.*mm;
    G4double Scint_dz = 5.*mm;
    G4double thickness_wrapping = 0.1*mm;
    // source holder
    G4double sourceHolder_dz = 0.5*mm;
    // SiPMs
    G4double SiPM_dx = 3.*mm;
    G4double SiPM_dy = 3.*mm;
    G4double SiPM_dz = 0.6*mm; // [https://www.ketek.net/wp-content/uploads/2018/12/KETEK-PM3325-WB-D0-Datasheet.pdf]
    // PCB and electronics - rough estimation
    G4double Electronics_dx = 20.*mm;
    G4double Electronics_dy = 20.*mm;
    G4double Electronics_dz = 2*mm;

    
    // World
    fWorldLength =  10*std::max(Scint_dx + Electronics_dx,std::max(Scint_dy + Electronics_dy,Scint_dz));
    G4Box * sWorld = new G4Box("World",                                 //name
                               0.5*fWorldLength, 0.5*fWorldLength, 0.5*fWorldLength); //dimensions
        
    logicWorld = new G4LogicalVolume(sWorld,                    //shape
                                 fWorldMater,                   //material
                                 "World");                      //name
    
    fPhysiWorld = new G4PVPlacement(0,                          //no rotation
                                    G4ThreeVector(),            //at (0,0,0)
                                    logicWorld,                 //logical volume
                                    "World",                    //name
                                    0,                          //mother volume
                                    false,                      //no boolean operation
                                    0);                         //copy number
    
    
    
    
    
    // (plastic) source holder
    if (DoSourceHolder) {
        posSourceHolder = G4ThreeVector(0 , 0 ,0 );
        solidSourceHolder = new G4Box("solidSourceHolder",
                                      1.*cm/2., 2.*cm/2. , sourceHolder_dz/2.);
        logicSourceHolder = new G4LogicalVolume(solidSourceHolder,
                                                                  fPlasticMater ,
                                                                  "logicSourceHolder",
                                                                  0,0,0);
        SourceHolderColor = G4Colour(0.8, 0.6, 0.7);
        SourceHolderVisAttributes = new G4VisAttributes(SourceHolderColor);
        logicSourceHolder->SetVisAttributes(SourceHolderVisAttributes);
        new G4PVPlacement(0,                                    // no rotation
                          G4ThreeVector(0,0,0),                 // at (x,y,z)
                          logicSourceHolder,                    // its logical volume
                          "placementSourceHolder",              // its name
                          logicWorld,                           // its mother  volume
                          false,                                // no boolean operations
                          0);                                   // copy number
        // <--- end source holder
    }


    // place-holder for the source position
    if (DoSourcePlaceHolder) {
        sourcePlaceHolder = new G4Sphere("sourcePlaceHolder",
                                                    0, 0.1 * mm,
                                                    0.*deg,360.*deg,0.*deg,180.*deg);
        logicSourcePlaceHolder = new G4LogicalVolume(sourcePlaceHolder,
                                                                       fWorldMater ,
                                                                       "logicSourcePlaceHolder",
                                                                       0,0,0);
        SourceRed = G4Colour(1.0, 0.1, 0.2);
        sourceVisAttributes = new G4VisAttributes(SourceRed);
        logicSourcePlaceHolder->SetVisAttributes(sourceVisAttributes);
        new G4PVPlacement(0,                                  // no rotation
                          G4ThreeVector(0,0,0),               // at (x,y,z)
                          logicSourcePlaceHolder,             // its logical volume
                          "placementSourcePlaceHolder",       // its name
                          logicSourceHolder,                  // its mother  volume
                          false,                              // no boolean operations
                          0);                                 // copy number
        // <--- end place-holder for the source position
    }
    
    
    // Scintillators from two sides of the source
    if (DoScintillators){
    ScintBlue = G4Colour(0.1, 0.2, 1.0);
    ScintVisAttributes = new G4VisAttributes(ScintBlue);
        
    positionScint_1 = G4ThreeVector(0 , 0 ,
                                    sourceHolder_dz/2. + Scint_dz/2. + thickness_wrapping/2. + Scint_source_dz/2.);
    solidScint_1 = new G4Box("solidScint_1", Scint_dx/2., Scint_dy/2., Scint_dz/2.);
    logicScint_1 = new G4LogicalVolume(solidScint_1,
                                       fDetectorMater ,
                                       "logicScint_1",
                                       0,0,0);
    
    logicScint_1->SetVisAttributes(ScintVisAttributes);
    new G4PVPlacement(0,              // no rotation
                      positionScint_1, // at (x,y,z)
                      logicScint_1,    // its logical volume
                      "Scintillator_1",       // its name
                      logicWorld,      // its mother  volume
                      false,           // no boolean operations
                      0);              // copy number
    G4ThreeVector positionScint_2 = G4ThreeVector(0 , 0 ,
                                                  -sourceHolder_dz/2. - Scint_dz/2. - thickness_wrapping/2. - Scint_source_dz/2. );
    solidScint_2 = new G4Box("solidScint_2", Scint_dx/2., Scint_dy/2., Scint_dz/2.);
    logicScint_2 = new G4LogicalVolume(solidScint_2,
                                       fDetectorMater ,
                                       "logicScint_2",
                                       0,0,0);
    
    logicScint_2->SetVisAttributes(ScintVisAttributes);
    new G4PVPlacement(0,              // no rotation
                      positionScint_2, // at (x,y,z)
                      logicScint_2,    // its logical volume
                      "Scintillator_2",       // its name
                      logicWorld,      // its mother  volume
                      false,           // no boolean operations
                      0);
    // <--- end scintillators
    }
            
    // SiPMs
    if (DoSiPMs){
        SiPMGreen = G4Colour(0.1, 0.9, 0.2);
        SiPMVisAttributes = new G4VisAttributes(SiPMGreen);
        
        positionSiPM_1 = G4ThreeVector(0 , 0 ,
                                        (sourceHolder_dz/2. + Scint_dz
                                         + thickness_wrapping/2.
                                         + Scint_source_dz/2.
                                         + SiPM_dz/2.) );
        solidSiPM_1 = new G4Box("solidSiPM_1", SiPM_dx/2., SiPM_dy/2., SiPM_dz/2.);
        logicSiPM_1 = new G4LogicalVolume(solidSiPM_1,
                                           fSiPMMater ,
                                           "logicSiPM_1",
                                           0,0,0);
        logicSiPM_1->SetVisAttributes(SiPMVisAttributes);
        
        new G4PVPlacement(0,              // no rotation
                          positionSiPM_1, // at (x,y,z)
                          logicSiPM_1,    // its logical volume
                          "SiPM_1",       // its name
                          logicWorld,      // its mother  volume
                          false,           // no boolean operations
                          0);              // copy number
        
        positionSiPM_2 = G4ThreeVector(0 , 0 ,
                                        -(sourceHolder_dz/2. + Scint_dz
                                         + thickness_wrapping/2.
                                         + Scint_source_dz/2.
                                         + SiPM_dz/2.) );
        solidSiPM_2 = new G4Box("solidSiPM_2", SiPM_dx/2., SiPM_dy/2., SiPM_dz/2.);
        logicSiPM_2 = new G4LogicalVolume(solidSiPM_2,
                                           fSiPMMater ,
                                           "logicSiPM_2",
                                           0,0,0);
        
        logicSiPM_2->SetVisAttributes(SiPMVisAttributes);
        new G4PVPlacement(0,              // no rotation
                          positionSiPM_2, // at (x,y,z)
                          logicSiPM_2,    // its logical volume
                          "SiPM_2",       // its name
                          logicWorld,      // its mother  volume
                          false,           // no boolean operations
                          0);              // copy number

    }
        
    
    if (DoElectronicPlates){
        ElectronicsColor = G4Colour(0.6, 0.8, 0.6);
        ElectronicsVisAttributes = new G4VisAttributes(ElectronicsColor);
        
        positionElectronics_1 = G4ThreeVector(0 , 0 ,
                                        (sourceHolder_dz/2. + Scint_dz
                                         + thickness_wrapping/2.
                                         + Scint_source_dz/2.
                                         + SiPM_dz/2.
                                         + Electronics_dz/2.) );
        solidElectronics_1 = new G4Box("solidSiPM_1", Electronics_dx/2., Electronics_dy/2., Electronics_dz/2.);
        logicElectronics_1 = new G4LogicalVolume(solidElectronics_1,
                                           fElectronicsMaterial ,
                                           "logicElectronics_1",
                                           0,0,0);
        logicElectronics_1->SetVisAttributes(ElectronicsVisAttributes);
        
        new G4PVPlacement(0,              // no rotation
                          positionElectronics_1, // at (x,y,z)
                          logicElectronics_1,    // its logical volume
                          "Electronics_1",       // its name
                          logicWorld,      // its mother  volume
                          false,           // no boolean operations
                          0);              // copy number
        
        positionElectronics_2 = G4ThreeVector(0 , 0 ,
                                        -(sourceHolder_dz/2. + Scint_dz
                                          + thickness_wrapping/2.
                                          + Scint_source_dz/2.
                                          + SiPM_dz/2.
                                          + Electronics_dz/2.) );
        solidElectronics_2 = new G4Box("solidElectronics_2", Electronics_dx/2., Electronics_dy/2., Electronics_dz/2.);
        logicElectronics_2 = new G4LogicalVolume(solidElectronics_2,
                                                 fElectronicsMaterial ,
                                           "logicElectronics_2",
                                           0,0,0);
        
        logicElectronics_2->SetVisAttributes(ElectronicsVisAttributes);
        new G4PVPlacement(0,              // no rotation
                          positionElectronics_2, // at (x,y,z)
                          logicElectronics_2,    // its logical volume
                          "Electronics_2",       // its name
                          logicWorld,      // its mother  volume
                          false,           // no boolean operations
                          0);              // copy number
    }

    
    
    
    G4cout << "done G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()." << G4endl;
    // return the root volume
    return fPhysiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
    //    G4cout << "\n Target : Length = " << G4BestUnit(fTargetLength,"Length")
    //    << " Radius = " << G4BestUnit(fTargetRadius,"Length")
    //    << " Material = " << fTargetMater->GetName();
    //    G4cout << "\n Detector : Length = " << G4BestUnit(fDetectorLength,"Length")
    //    << " Tickness = " << G4BestUnit(fDetectorThickness,"Length")
    //    << " Material = " << fDetectorMater->GetName() << G4endl;
    //    G4cout << "\n" << fTargetMater << "\n" << fDetectorMater << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetMaterial(G4String materialChoice)
{
    // search the material by its name
    G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);
    
    if (pttoMaterial) {
        fTargetMater = pttoMaterial;
        if(fLogicTarget) { fLogicTarget->SetMaterial(fTargetMater); }
        G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    } else {
        G4cout << "\n--> warning from DetectorConstruction::SetTargetMaterial : "
        << materialChoice << " not found" << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetectorMaterial(G4String materialChoice)
{
    // search the material by its name
    G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);
    
    if (pttoMaterial) {
        fDetectorMater = pttoMaterial;
        if(fLogicDetector) { fLogicDetector->SetMaterial(fDetectorMater); }
        G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    } else {
        G4cout << "\n--> warning from DetectorConstruction::SetDetectorMaterial : "
        << materialChoice << " not found" << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetRadius(G4double value)
{
    fTargetRadius = value;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTargetLength(G4double value)
{
    fTargetLength = value;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetectorThickness(G4double value)
{
    fDetectorThickness = value;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetDetectorLength(G4double value)
{
    fDetectorLength = value;
    G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetTargetLength()
{
    return fTargetLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetTargetRadius()
{
    return fTargetRadius;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::GetTargetMaterial()
{
    return fTargetMater;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DetectorConstruction::GetLogicTarget()
{
    return fLogicTarget;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetDetectorLength()
{
    return fDetectorLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetDetectorThickness()
{
    return fDetectorThickness;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::GetDetectorMaterial()
{
    return fDetectorMater;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DetectorConstruction::GetLogicDetector()
{
    return fLogicDetector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
