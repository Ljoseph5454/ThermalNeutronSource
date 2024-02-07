//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The./  Geant4 software  is  copyright of the Copyright Holders  of *
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
#include "G4Polycone.hh"
#include "G4UnionSolid.hh"
#include "G4MultiUnion.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"

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

#include <math.h>

#define PI 3.14159265
using namespace std;

//Ryan added package
#include "G4OpticalSurface.hh"


DMXDetectorConstruction::DMXDetectorConstruction()
{


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

}

DMXDetectorConstruction::~DMXDetectorConstruction()
{
  delete theUserLimitsForRoom;
  delete theUserLimitsForDetector;
  delete detectorMessenger;
}

void DMXDetectorConstruction::DefineMaterials()
{

 #include "DMXDetectorMaterial.icc"

}

G4VPhysicalVolume* DMXDetectorConstruction::Construct()
{  
  DefineMaterials();
  // Get nist material manager
  //G4NistManager* nist = G4NistManager::Instance();
  
  // Envelope parameters
  //
  G4double S_l = 7.5*cm, S_w=5*cm, V_l = 3*cm, P_w = 25*cm, P_l = 10*cm, P_p = (S_l+P_l)+5*cm, A_l = 30.48*cm;
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
  logicS = new G4LogicalVolume(solidS, vacuum_mat, "logicS"); //HDPE_mat                   
  physS = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.5*(P_p-P_w)), logicS, "physS", logicWorld, false, 0);  
 
  // Sapphire Window
  G4Box* solidSap = new G4Box("solidSap", 0.5*S_w, 0.5*S_w, 0.5*S_l); 
  logicSap = new G4LogicalVolume(solidSap, vacuum_mat, "logicSap");  //sapphire_mat                  
  physSap = new G4PVPlacement(0, G4ThreeVector(0.,0.,(0.5*(V_l+S_l)+P_l)-0.5*(P_p-P_w)), logicSap, "physSap", logicS, false, 0);  

  // Test Argon
  G4Box* solidAr = new G4Box("solidAr", 0.5*A_l, 0.5*A_l, 0.5*A_l); 
  logicAr = new G4LogicalVolume(solidAr, vacuum_mat, "logicAr");  //sapphire_mat                  
  physAr = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.5*(P_w+V_l+P_p)+0.5*(P_p-P_w)+1*m+0.5*A_l), logicAr, "physAr", logicWorld, false, 0);  

  // SD before
 // G4Box* solidSD1 = new G4Box("solidSD1", S_l, S_l, 1*mm); 
  //logicSD1 = new G4LogicalVolume(solidSD1, vacuum_mat, "logicSD1");                    
  //physSD1 = new G4PVPlacement(0, G4ThreeVector(0.,0.,5*cm), logicSD1, "physSD1", logicWorld, false, 0);   
  
  // Empty Inside
  G4Box* solidSD2 = new G4Box("solidSD2", 0.5*V_l, 0.5*V_l, 0.5*V_l); 
  logicSD2 = new G4LogicalVolume(solidSD2, vacuum_mat, "logicSD2");                    
  physSD2 = new G4PVPlacement(0, G4ThreeVector(0.,0.,-0.5*(P_p-P_w)), logicSD2, "physSD2", logicS, false, 0);  //0.5*(P_l+S_l-P_w)

  // Test Sphere
  /*G4Sphere* solidSphere = new G4Sphere("solidSphere", 0.0*cm, 10.0*cm, 0.0 * deg, 360.0 *deg, 0.0 * deg, 180 *deg); 
  logicSphere = new G4LogicalVolume(solidSphere, vacuum_mat, "logicSphere");  //sapphire_mat                  
  physSphere = new G4PVPlacement(0, G4ThreeVector(50*cm,0.,0.), logicSphere, "physSphere", logicWorld, false, 0); */ 
 
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

   if (LXeSD.Get() == 0)                                           // Aquí detecto los eventos del argón
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
 //if(logicSD2)
   //   SetSensitiveDetector(logicSD2,LXeSD.Get());
    /*if (pmtSD.Get() == 0)                                        //Aquí detecto los eventos en el SiPM
    {
      G4String name="/DMXDet/pmtSD";
      DMXPmtSD* aSD = new DMXPmtSD(name);
      pmtSD.Put(aSD);
    }
  G4SDManager::GetSDMpointer()->AddNewDetector(pmtSD.Get());
  if (SIPM1_Si_log)
    SetSensitiveDetector(SIPM1_Si_log,pmtSD.Get());
  if (SIPM2_Si_log)
    SetSensitiveDetector(SIPM2_Si_log,pmtSD.Get());
  if (SIPM3_Si_log)
    SetSensitiveDetector(SIPM3_Si_log,pmtSD.Get());
  if (SIPM4_Si_log)
    SetSensitiveDetector(SIPM4_Si_log,pmtSD.Get());
  if (SIPM5_Si_log)
    SetSensitiveDetector(SIPM5_Si_log,pmtSD.Get());


*/


    return;

}


/*
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
DMXDetectorConstruction::DMXDetectorConstruction()
{


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

  // DefineField();

  // make colours
  G4Colour  white   (1.0, 1.0, 1.0) ;
  G4Colour  grey    (0.5, 0.5, 0.5) ;
  G4Colour  lgrey   (.85, .85, .85) ;
  G4Colour  red     (1.0, 0.0, 0.0) ;
  G4Colour  blue    (0.0, 0.0, 1.0) ;
  G4Colour  cyan    (0.0, 1.0, 1.0) ;
  G4Colour  magenta (1.0, 0.0, 1.0) ;
  G4Colour  yellow  (1.0, 1.0, 0.0) ;
  G4Colour  orange  (.75, .55, 0.0) ;
  G4Colour  lblue   (0.0, 0.0, .75) ;
  G4Colour  lgreen  (0.0, .75, 0.0) ;
  G4Colour  green   (0.0, 1.0, 0.0) ;
  G4Colour  brown   (0.7, 0.4, 0.1) ;


  //  un-used colours:
  //  G4Colour  black   (0.0, 0.0, 0.0) ;



  // Universe


  G4double worldWidth  = 2*m ;
  G4double worldLength = 2*m ;
  G4double worldHeight = 1*m ;

  G4Box* world_box = new G4Box
     ("world_box", 0.5*worldWidth, 0.5*worldLength, 0.5*worldHeight );
  world_log  = new G4LogicalVolume(world_box, vacuum_mat, "world_log");
  world_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),"world_phys", world_log, NULL, false,0);

G4double ArcontainerWidth = 1*m ;
G4double ArcontainerLength = 1*m ;
G4double ArcontainerHeight = 3*cm ;

G4Box* Ar_box = new G4Box("Ar_box",0.5*ArcontainerWidth, 0.5*ArcontainerLength, 0.5*ArcontainerHeight);
Arbox_log = new G4LogicalVolume(Ar_box, LAr_mat,"Arbox_log");
Arbox_phys =  new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), Arbox_log, "Arbox_phys", world_log, false, 0);
G4VisAttributes* argon_blue = new G4VisAttributes(blue);
argon_blue->SetVisibility(true);
Arbox_log -> SetVisAttributes(argon_blue);

//add opticalphoton cut
//G4double minEkin = 10*MeV;
//Arbox_log->SetUserLimits(new G4UserLimits(DBL_MAX, DBL_MAX, DBL_MAX,
 //                                         minEkin));


return world_phys;



}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXDetectorConstruction::ConstructSDandField()
{
    // ......................................................................
    // sensitive detectors ..................................................
    // ......................................................................

   if (LXeSD.Get() == 0)                                           // Aquí detecto los eventos del argón
      {
        G4String name="/DMXDet/LXeSD";
        DMXScintSD* aSD = new DMXScintSD(name);
        LXeSD.Put(aSD);
      }
    G4SDManager::GetSDMpointer()->AddNewDetector(LXeSD.Get());
 if(Arbox_log)
      SetSensitiveDetector(Arbox_log,LXeSD.Get());


    /*if (pmtSD.Get() == 0)                                        //Aquí detecto los eventos en el SiPM
    {
      G4String name="/DMXDet/pmtSD";
      DMXPmtSD* aSD = new DMXPmtSD(name);
      pmtSD.Put(aSD);
    }
  G4SDManager::GetSDMpointer()->AddNewDetector(pmtSD.Get());
  if (SIPM1_Si_log)
    SetSensitiveDetector(SIPM1_Si_log,pmtSD.Get());
  if (SIPM2_Si_log)
    SetSensitiveDetector(SIPM2_Si_log,pmtSD.Get());
  if (SIPM3_Si_log)
    SetSensitiveDetector(SIPM3_Si_log,pmtSD.Get());
  if (SIPM4_Si_log)
    SetSensitiveDetector(SIPM4_Si_log,pmtSD.Get());
  if (SIPM5_Si_log)
    SetSensitiveDetector(SIPM5_Si_log,pmtSD.Get());


*/


//   return;

//}