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
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: alexander.howard@cern.ch
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// DetectorConstruction program
// --------------------------------------------------------------

#include "DMXDetectorConstruction.hh"
#include "DMXDetectorMessenger.hh"

#include "DMXScintSD.hh"
#include "DMXPmtSD.hh"


#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4UnitsTable.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4FieldManager.hh"
#include "G4UniformElectricField.hh"
#include "G4TransportationManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4EqMagElectricField.hh"
#include "G4ClassicalRK4.hh"
#include "G4ChordFinder.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UserLimits.hh"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
//Logan
#include "G4NCrystal/G4NCrystal.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
DMXDetectorConstruction::DMXDetectorConstruction()  
{
  // create commands for interactive definition of time cuts:
  detectorMessenger = new DMXDetectorMessenger(this);

  theUserLimitsForRoom     = 0; 
  theUserLimitsForDetector = 0; 
  // default time cut = infinite
  //  - note also number of steps cut in stepping action = MaxNoSteps
  theMaxTimeCuts      = DBL_MAX;
  theMaxStepSize      = DBL_MAX;
  theDetectorStepSize = DBL_MAX;
  theRoomTimeCut      = 1000. * nanosecond;
  theMinEkine         = 250.0*eV; // minimum kinetic energy required in volume
  theRoomMinEkine     = 250.0*eV; // minimum kinetic energy required in volume
  
  //Zero the G4Cache objects to contain logical volumes
  LXeSD.Put(0);
  pmtSD.Put(0);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
DMXDetectorConstruction::~DMXDetectorConstruction() 
{
  delete theUserLimitsForRoom;
  delete theUserLimitsForDetector;
  delete detectorMessenger;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DMXDetectorConstruction::DefineMaterials() 
{

#include "DMXDetectorMaterial.icc"

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4VPhysicalVolume* DMXDetectorConstruction::Construct() {

  DefineMaterials();
  // Get nist material manager
  //G4NistManager* nist = G4NistManager::Instance();
  
  //Logan
  G4Material * mat_aluminium = G4NCrystal::createMaterial("Al_sg225.ncmat");
  G4Material* sapphireNCrystal_mat = G4NCrystal::createMaterial("Al2O3_sg167_Corundum.ncmat;bragg=0"); 
  G4Material* HDPENCrystal_mat = G4NCrystal::createMaterial("Polyethylene_CH2.ncmat;density=0.96gcm3");
  G4Material* vacuumNCrystal_mat = G4NCrystal::createMaterial("void.ncmat");

  // Envelope parameters
  //
  G4double S_l = 7.5*cm, S_w=5*cm, V_l = 3*cm, P_w = 25*cm, P_l = 3*cm, P_p = (S_l+P_l)+10*cm, A_l = 30.48*cm;
  //G4Material* env_mat = nist->FindOrBuildMaterial("G4_WATER");
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeXY = 2.5*m;
  G4double world_sizeZ  = 2.5*m;
  //G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       world_sizeXY, world_sizeXY, world_sizeZ);     //its size
      
    logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        vacuum_mat,           //its material
                        "logicWorld");            //its name
                                   
    physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "physWorld",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

                     

  // Filter
  G4Box* solidPoly = new G4Box("solidPoly", 0.5*(2*P_w+V_l), 0.5*(2*P_w+V_l), 0.5*(P_w+V_l+P_p)); //0.5*(P_w+V_l+P_l+S_l)
  G4Box* vacuumNotch = new G4Box("vacuumNotch", 0.5*S_w, 0.5*S_w, 0.5*(P_p-P_l-S_l)); //0.5*(P_w+V_l+P_l+S_l) //0.5*(P_p-P_l-S_l)
  G4ThreeVector trans = G4ThreeVector(0.,0.,0.5*(V_l+P_l+S_l+P_w));
  G4RotationMatrix rot = G4RotationMatrix(0.,0.,0.);
  G4Transform3D transform = G4Transform3D(rot,trans);
  G4VSolid* solidS = new G4SubtractionSolid("solidS" ,solidPoly, vacuumNotch, transform);
  logicS = new G4LogicalVolume(solidS, HDPENCrystal_mat, "logicS"); //HDPENCrystal_mat                  
  physS = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.5*(P_p-P_w)), logicS, "physS", logicWorld, false, 0);  
 
  // Sapphire Window
  G4Box* solidSap = new G4Box("solidSap", 0.5*S_w, 0.5*S_w, 0.5*S_l); 
  logicSap = new G4LogicalVolume(solidSap, sapphireNCrystal_mat, "logicSap");  //sapphireNCrystal_mat                
  physSap = new G4PVPlacement(0, G4ThreeVector(0.,0.,(0.5*(V_l+S_l)+P_l)-0.5*(P_p-P_w)), logicSap, "physSap", logicS, false, 0); 

  // Test Argon
  G4Box* solidAr = new G4Box("solidAr", 0.5*A_l, 0.5*A_l, 0.5*A_l); 
  logicAr = new G4LogicalVolume(solidAr, LAr_mat, "logicAr");  //sapphire_mat                  
  physAr = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.5*(P_w+V_l+P_p)+0.5*(P_p-P_w)+1*m+0.5*A_l), logicAr, "physAr", logicWorld, false, 0);

  // SD before
 // G4Box* solidSD1 = new G4Box("solidSD1", S_l, S_l, 1*mm); 
  //logicSD1 = new G4LogicalVolume(solidSD1, vacuum_mat, "logicSD1");                    
  //physSD1 = new G4PVPlacement(0, G4ThreeVector(0.,0.,5*cm), logicSD1, "physSD1", logicWorld, false, 0);   
  
  // Empty Inside
  G4Box* solidSD2 = new G4Box("solidSD2", 0.5*V_l, 0.5*V_l, 0.5*V_l); 
  logicSD2 = new G4LogicalVolume(solidSD2, vacuumNCrystal_mat, "logicSD2");                    
  physSD2 = new G4PVPlacement(0, G4ThreeVector(0.,0.,-0.5*(P_p-P_w)), logicSD2, "physSD2", logicS, false, 0);  //0.5*(P_l+S_l-P_w)

 /* // Sapphire Window
  G4Box* solidWindow = new G4Box("solidWindow", 0.25*(S_l-V_l), 0.5*V_l, 0.5*V_l); 
  logicWindow = new G4LogicalVolume(solidWindow, vacuum_mat, "logicWindow");                    
  physWindow = new G4PVPlacement(0, G4ThreeVector(-0.25*(S_l+V_l),0.,0.), logicWindow, "physWindow", logicWorld, false, 0); */

  //
  //always return the physical World
  //
  return physWorld;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXDetectorConstruction::ConstructSDandField()
{
  // ......................................................................
  // sensitive detectors ..................................................
  // ......................................................................

  if (LXeSD.Get() == 0) 
    {    
      G4String name="/DMXDet/LXeSD";
      DMXScintSD* aSD = new DMXScintSD(name);
      LXeSD.Put(aSD);
    }
  G4SDManager::GetSDMpointer()->AddNewDetector(LXeSD.Get()); 
  if(logicS){
      SetSensitiveDetector(logicS,LXeSD.Get());
      SetSensitiveDetector(logicSD2,LXeSD.Get());
      SetSensitiveDetector(logicWorld,LXeSD.Get());
      SetSensitiveDetector(logicSap,LXeSD.Get());
      SetSensitiveDetector(logicAr,LXeSD.Get());}
  /*if (LXe_log)    
    SetSensitiveDetector(LXe_log,LXeSD.Get());

  if (pmtSD.Get() == 0)
    {
      G4String name="/DMXDet/pmtSD";
      DMXPmtSD* aSD = new DMXPmtSD(name);
      pmtSD.Put(aSD);
    }
  G4SDManager::GetSDMpointer()->AddNewDetector(pmtSD.Get()); 
  if (phcath_log)
    SetSensitiveDetector(phcath_log,pmtSD.Get());*/

  return;
}
 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// specific method to G4UserLimits:= SetUserMinEkine
void DMXDetectorConstruction::SetRoomEnergyCut(G4double val)
{
  // set minimum charged particle energy cut - NB: for ROOM
  theRoomMinEkine = val;
  if (theUserLimitsForRoom != 0) 
    {
      theUserLimitsForRoom->SetUserMinEkine(val); 
      G4cout << " Changing Room energy cut to: " << G4BestUnit(val,"Energy")
	     << G4endl;
    }
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// specific method to G4UserLimits:= SetUserMinEkine
void DMXDetectorConstruction::SetEnergyCut(G4double val)
{
  // set minimum charged particle energy cut - NB: for Xenon Detector
  theMinEkine = val;
  if (theUserLimitsForDetector != 0) 
    {
      theUserLimitsForDetector->SetUserMinEkine(val);
      G4cout << "Changing Detector energy cut to: " << G4BestUnit(val,"Energy")
	     << G4endl;
    }
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// specific method to G4UserLimits:= SetUserMaxTime
void DMXDetectorConstruction::SetRoomTimeCut(G4double val)
{
  // set room time cut:
  theRoomTimeCut = val;
  if (theUserLimitsForRoom != 0) 
    {
      theUserLimitsForRoom->SetUserMaxTime(val);
      G4cout << " Changing Room Time cut to: " << G4BestUnit(val,"Time")
	     << G4endl;
    }
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// specific method to G4UserLimits:= SetUserMaxTime
void DMXDetectorConstruction::SetTimeCut(G4double val)
{
  // set detector time cut:
  theMaxTimeCuts = val;
  if (theUserLimitsForDetector != 0) 
    {
      theUserLimitsForDetector->SetUserMaxTime(val);
      G4cout << " Changing Detector Time cut to: " << G4BestUnit(val,"Time")
	     << G4endl;
    }
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



