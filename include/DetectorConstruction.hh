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
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"

class G4LogicalVolume;
class G4Material;
class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    DetectorConstruction();
   ~DetectorConstruction();

  public:
  
    virtual G4VPhysicalVolume* Construct();
    
    void SetTargetLength (G4double value);
    void SetTargetRadius (G4double value);
    void SetTargetMaterial (G4String);
    
    void SetDetectorLength(G4double value);           
    void SetDetectorThickness(G4double value);  
    void SetDetectorMaterial(G4String);               
                   
    void PrintParameters();
    
    std::string VolumeName(int iVol=0){
        std::string name = "unknown-volume";
        switch (iVol) {
            case 0:
                name = "source-holder";
                break;
            case 1:
                name = "scint-1";
                break;
            case 2:
                name = "scint-2";
                break;
            case 3:
                name = "world";
                break;

            default:
                name = "unknown-volume";
                break;
        }
        return name;
    }
    
    
  public:
      
    G4double GetTargetLength();
    G4double GetTargetRadius();
    G4Material* GetTargetMaterial();       
    G4LogicalVolume* GetLogicTarget();
    
    G4double GetDetectorLength();
    G4double GetDetectorThickness();
    G4Material* GetDetectorMaterial();                 
    G4LogicalVolume* GetLogicDetector();      
         
    

    
    G4LogicalVolume * logicScint_1;
    G4LogicalVolume * logicScint_2;
    
    
    G4LogicalVolume *       GetLogicScint_1(){ return logicScint_1; };
    G4LogicalVolume *       GetLogicScint_2(){ return logicScint_2; };
    G4LogicalVolume *  GetLogicSourceHolder(){ return logicSourceHolder; };
    G4LogicalVolume *         GetLogicWorld(){ return logicWorld; };
    
    
    
    
  private:
    G4LogicalVolume*   logicWorld;
    G4double           fTargetLength; 
    G4double           fTargetRadius;
    G4Material*        fTargetMater;
    G4LogicalVolume*   fLogicTarget;
                 
    G4double           fDetectorLength;
    G4double           fDetectorThickness;
    G4Material*        fDetectorMater;
    G4LogicalVolume*   fLogicDetector;
               
    G4double           fWorldLength;
    G4double           fWorldRadius;
    G4Material*        fWorldMater;     
    G4VPhysicalVolume* fPhysiWorld;
                
    DetectorMessenger* fDetectorMessenger;
    
    G4Material*         fPlasticMater;

    
    
    G4ThreeVector posSourceHolder;
    G4ThreeVector positionScint_1, positionScint_2;
    
    G4Box * solidSourceHolder;
    G4Sphere * sourcePlaceHolder;
    G4Box * solidScint_1, * solidScint_2;
    
    G4LogicalVolume * logicSourceHolder;
    G4LogicalVolume * logicSourcePlaceHolder;
    
    G4Colour SourceHolderColor;
    G4Colour SourceRed;
    G4Colour ScintBlue;
    
    G4VisAttributes* SourceHolderVisAttributes;
    G4VisAttributes* sourceVisAttributes;
    G4VisAttributes* ScintVisAttributes;

    
  private:
    
    void               DefineMaterials();
    G4VPhysicalVolume* ConstructVolumes();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

