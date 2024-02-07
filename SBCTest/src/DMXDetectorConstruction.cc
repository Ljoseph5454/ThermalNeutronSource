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
  G4double worldHeight = 4*m ;

  G4Box* world_box = new G4Box
     ("world_box", 0.5*worldWidth, 0.5*worldLength, 0.5*worldHeight );
  world_log  = new G4LogicalVolume(world_box, lab_mat, "world_log");
  world_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.),
     "world_phys", world_log, NULL, false,0);


 // Vacuum Jacket 

  //G4double vacuum_vessel_hz = 163.4744*cm ;
  G4double vacuum_vessel_r = 30.48*cm;
  G4RotationMatrix* yRot = new G4RotationMatrix;               
  yRot->rotateY(22.5*deg);
  
  G4RotationMatrix* y180Rot = new G4RotationMatrix;               
  y180Rot->rotateY(180*deg);
  
  G4RotationMatrix* xRot = new G4RotationMatrix;               
  xRot->rotateX(22.5*deg);

  G4RotationMatrix* zzRot = new G4RotationMatrix;  
  zzRot->rotateZ(120.0*deg); 

  G4RotationMatrix* zz2Rot = new G4RotationMatrix;  
  zz2Rot->rotateZ(240.0*deg); 
  
  G4RotationMatrix* yzRot = new G4RotationMatrix;               
  yzRot->rotateY(22.5*deg);
  yzRot->rotateZ(120.0*deg); 
  
  
  //New axis for camera rotation
  G4ThreeVector AxisOfRotation = G4ThreeVector(0.866025,-0.5,0).unit();
  
  G4RotationMatrix * xyRot = new G4RotationMatrix();
  xyRot -> rotate(22.5*deg, AxisOfRotation);
  xyRot -> rotateZ(30.0*deg);
  
  
  G4ThreeVector AxisOfRotation1 = G4ThreeVector(-0.866025,-0.5,0).unit();
  
  G4RotationMatrix * xy1Rot = new G4RotationMatrix();
  xy1Rot -> rotate(22.5*deg, AxisOfRotation1);
  xy1Rot -> rotateZ(60.0*deg);
  
  
  
  
  G4double* z00 = new G4double[8];
  z00[0]= -72.465*cm;                     //Original 52.465
  z00[1]= -72.465*cm + 2*cm;
  z00[2]= -72.465*cm + 4*cm;
  z00[3]= -72.465*cm + 6*cm;
  z00[4]= -72.465*cm + 8*cm;
  z00[5]= -72.465*cm + 10*cm;
  z00[6]= 0;
  z00[7]= 111.0094*cm - 11*cm + 2*cm;



  G4double* rInn00 = new G4double[8];
  rInn00[0]=0.;
  rInn00[1]=0.;
  rInn00[2]=0.;
  rInn00[3]=0.;
  rInn00[4]=0.;
  rInn00[5]=0.;
  rInn00[6]=0.;
  rInn00[7]=0.;


  G4double* rOut00 = new G4double[8];
  rOut00[0]=0.;
  rOut00[1]=sqrt(2/0.01076)*cm;
  rOut00[2]=sqrt(4/0.01076)*cm;
  rOut00[3]=sqrt(6/0.01076)*cm;
  rOut00[4]=sqrt(8/0.01076)*cm;
  rOut00[5]=vacuum_vessel_r;
  rOut00[6]=vacuum_vessel_r;
  rOut00[7]=vacuum_vessel_r;




  G4Polycone* Cil_1 = new G4Polycone("Cil_1", 0.*deg, 360.*deg,
  8, z00, rInn00, rOut00);


 // G4Tubs* Cil_2 = new G4Tubs("Cil_2",0, 5*cm, 4*cm, 0.*deg, 360.*deg);


 //Old hole in VJ for camera
 // G4VSolid* Vacuum_vessel_Cil= new G4SubtractionSolid("Vacuum_vessel_Cil", Cil_1, Cil_2, yRot, G4ThreeVector(-25*cm,0.,111.0094*cm)); 



  Vacuum_vessel_log = new G4LogicalVolume(Cil_1, vessel_mat, "Vacuum_vessel_log");
  Vacuum_vessel_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), Vacuum_vessel_log, "Vacuum_vessel_phys", world_log, false, 0) ;

  G4VisAttributes* ssteal_blue= new G4VisAttributes(blue);
  ssteal_blue->SetVisibility(true);
  Vacuum_vessel_log->SetVisAttributes(ssteal_blue);


  // Holder 1 Vacuum Vessel

G4Box* leg_box = new G4Box
     ("leg_box", 0.5*5*cm, 0.5*10*cm, 0.5*111*cm );

G4Box* leg_end_box = new G4Box
     ("leg_end_box", 0.5*13*cm, 0.5*13*cm, 0.5*2*cm );

G4VSolid* leg= new G4UnionSolid("leg", leg_box, leg_end_box, 0, G4ThreeVector(3*cm,0.,-56.5*cm)); 


  Holder1_VV_log = new G4LogicalVolume(leg, vessel_mat, "Holder1_VV_log");
  Holder1_VV_phys = new G4PVPlacement(0, G4ThreeVector(30.48*cm + 2.5*cm,0.,-68*cm), Holder1_VV_log, "Holder1_VV_phys", world_log, false, 0) ;

 Holder1_VV_log->SetVisAttributes(ssteal_blue);

  // Holder 2 Vacuum Vessel
  G4RotationMatrix* z3Rot = new G4RotationMatrix;
  z3Rot->rotateZ(120*deg);

  Holder2_VV_log = new G4LogicalVolume(leg, vessel_mat, "Holder2_VV_log");
  Holder2_VV_phys = new G4PVPlacement(z3Rot, G4ThreeVector(32.98*cos(2*3.141592/3)*cm,-32.98*sin(2*3.141592/3)*cm,-68*cm), Holder2_VV_log, "Holder2_VV_phys", world_log, false, 0) ;

 Holder2_VV_log->SetVisAttributes(ssteal_blue);


  // Holder 3 Vacuum Vessel
  G4RotationMatrix* z4Rot = new G4RotationMatrix;
  z4Rot->rotateZ(240*deg);

  Holder3_VV_log = new G4LogicalVolume(leg, vessel_mat, "Holder3_VV_log");
  Holder3_VV_phys = new G4PVPlacement(z4Rot, G4ThreeVector(32.98*cos(4*3.141592/3)*cm,-32.98*sin(4*3.141592/3)*cm,-68*cm), Holder3_VV_log, "Holder3_VV_phys", world_log, false, 0) ;

 Holder3_VV_log->SetVisAttributes(ssteal_blue);



 // Inside Vacuum Vessel

  G4double* z01 = new G4double[8];
  z01[0]= -72.465*cm + 1.04*cm;
  z01[1]= -72.465*cm + 1.04*cm + 2*cm;
  z01[2]= -72.465*cm + 1.04*cm + 4*cm;
  z01[3]= -72.465*cm + 1.04*cm + 6*cm;
  z01[4]= -72.465*cm + 1.04*cm + 8*cm;
  z01[5]= -72.465*cm + 1.04*cm + 10*cm;
  z01[6]= 0;
  z01[7]= 111.0094*cm - 2.54*cm - 11*cm + 2*cm;           //Se debe revisar el tamaño total del VJ, hacia abajo


  G4double* rInn01 = new G4double[8];
  rInn01[0]=0.;
  rInn01[1]=0.;
  rInn01[2]=0.;
  rInn01[3]=0.;
  rInn01[4]=0.;
  rInn01[5]=0.;
  rInn01[6]=0.;
  rInn01[7]=0.;


  G4double* rOut01 = new G4double[8];
  rOut01[0]=0.;
  rOut01[1]=sqrt(2/0.0128099)*cm;
  rOut01[2]=sqrt(4/0.0128099)*cm;
  rOut01[3]=sqrt(6/0.0128099)*cm;
  rOut01[4]=sqrt(8/0.0128099)*cm;
  rOut01[5]=vacuum_vessel_r - 2.54*cm;
  rOut01[6]=vacuum_vessel_r - 2.54*cm;
  rOut01[7]=vacuum_vessel_r - 2.54*cm;



  G4Polycone* Cil_3 = new G4Polycone("Cil_3", 0.*deg, 360.*deg,
  8, z01, rInn01, rOut01);

  //Old hole for camera
 // G4VSolid* Inside_vacuum_vessel_Cil= new G4SubtractionSolid("Inside_vacuum_vessel_Cil", Cil_3, Cil_2, yRot, G4ThreeVector(-25*cm,0.,111.0094*cm)); 

 
  Inside_vacuum_vessel_log = new G4LogicalVolume(Cil_3, vacuum_mat, "Inside_vacuum_vessel_log");
 Inside_vacuum_vessel_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), Inside_vacuum_vessel_log, "Inside_vacuum_vessel_phys", Vacuum_vessel_log, false, 0) ;


 //Camera Port Definition of geometry

   G4double* z_cp = new G4double[4];
  z_cp[0]= 0*cm;
  z_cp[1]= 7.9743*cm;
  z_cp[2]= 7.9743*cm + 1.27*cm;
  z_cp[3]= 7.9743*cm + 2.54*cm;


  G4double* rOut_cp = new G4double[4];
  rOut_cp[0]=6.985*cm;
  rOut_cp[1]=6.985*cm;
  rOut_cp[2]=6.0325*cm;
  rOut_cp[3]=5.08*cm;


  G4double* rInn_cp = new G4double[4];
  rInn_cp[0]=6.0325*cm;
  rInn_cp[1]=3.2385*cm;
  rInn_cp[2]=2.794*cm;
  rInn_cp[3]=2.794*cm;

  G4Polycone* Camera_Port = new G4Polycone("Camera_Port", 0.*deg, 360.*deg,
  4, z_cp, rInn_cp, rOut_cp); 

  



  //Pressure vessel                                          

  G4double cil1_hz = (71.85-5.43)*cm ;
  G4double r1 = 20.3*cm ;
  G4double r8 = 55.88*0.5*cm;

  G4double* z1 = new G4double[11];
  
  z1[0]=-19*cm;
  z1[1]=-13.286*cm;
  z1[2]=-13.285*cm; 
  z1[3]= 0.;
  z1[4]= cil1_hz;
  z1[5]= z1[4]+2*r1*cos(25 * PI / 180.0)-sqrt(3)*r1;
  z1[6]= z1[4]+2*r1*cos(20 * PI / 180.0)-sqrt(3)*r1;
  z1[7]= z1[4]+2*r1*cos(15 * PI / 180.0)-sqrt(3)*r1;
  z1[8]= z1[4]+2*r1*cos(10 * PI / 180.0)-sqrt(3)*r1;
  z1[9]= z1[4]+2*r1*cos(5 * PI / 180.0)-sqrt(3)*r1;
  z1[10]= z1[4]+2*r1-sqrt(3)*r1;

  //G4cout<<"Altura P_Vessel\n"<<z1[10]<<G4endl;    // 71.8594 cm


  G4double* rInn = new G4double[12];  
  
  rInn[0]=0.;
  rInn[1]=0.;
  rInn[2]=0.;
  rInn[3]=0.;
  rInn[4]=0.;
  rInn[5]=0.;
  rInn[6]=0.;
  rInn[7]=0.;
  rInn[8]=0.;
  rInn[9]=0.;
  rInn[10]=0.;
  rInn[11]=0.;

  G4double* rOut1 = new G4double[11];
  
  rOut1[0]=r8;
  rOut1[1]=r8;
  rOut1[2]=r1;
  rOut1[3]=r1;
  rOut1[4]=r1;
  rOut1[5]=2*r1*sin(25 * PI / 180.0);
  rOut1[6]=2*r1*sin(20 * PI / 180.0);
  rOut1[7]=2*r1*sin(15 * PI / 180.0);
  rOut1[8]=2*r1*sin(10 * PI / 180.0);
  rOut1[9]=2*r1*sin(5 * PI / 180.0);
  rOut1[10]=0.;

  G4Polycone* contenedor1 = new G4Polycone("contenedor1", 0.*deg, 360.*deg,
  11, z1, rInn, rOut1);
  
  
 // G4VSolid* PV_Geometry = new G4SubtractionSolid("Inter_cp_PV", contenedor1, Camera_Port, yRot, G4ThreeVector(-3.0*2.54*cm,0.,72*cm-2*2.54*cm));  //View port hole in PV
  

  pressure_vessel_log = new G4LogicalVolume(contenedor1, vessel_mat, "pressure_vessel_log");
  pressure_vessel_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), pressure_vessel_log, "pressure_vessel_phys", Inside_vacuum_vessel_log, false, 0) ;

  pressure_vessel_log->SetVisAttributes(ssteal_blue);






  //PV spool Piece

  G4double d1 = 1.3*cm;
  G4double r2 = r1-d1 ;
  G4double r9 = 18.415*0.5*cm ;


  G4double* z02 = new G4double[6];
  z02[0]= -19*cm;
  z02[1]= -25.667*cm;
  z02[2]= -25.668*cm;
  z02[3]= -54.475*cm;
  z02[4]= -54.476*cm;
  z02[5]= -58.769*cm;



  G4double* rInn02 = new G4double[6];
  rInn02[0]=r9;
  rInn02[1]=r9;
  rInn02[2]=r9;
  rInn02[3]=r9;
  rInn02[4]=r9;
  rInn02[5]=r9;
  


  G4double* rOut02 = new G4double[6];
  rOut02[0]=r8;
  rOut02[1]=r8;
  rOut02[2]=20.32*0.5*cm;
  rOut02[3]=20.32*0.5*cm;
  rOut02[4]=41.91*0.5*cm;
  rOut02[5]=41.91*0.5*cm;


  G4Polycone* PV_spool = new G4Polycone("PV_spool", 0.*deg, 360.*deg,
  6, z02, rInn02, rOut02);


  PV_spool_log = new G4LogicalVolume(PV_spool, vessel_mat, "PV_spool_log");
  PV_spool_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), PV_spool_log, "PV_spool_phys", Inside_vacuum_vessel_log, false, 0) ;

  PV_spool_log->SetVisAttributes(ssteal_blue);



 //Bellows Weldment
 

  G4double* z03 = new G4double[8];
  z03[0]= 0.;
  z03[1]=4.292*cm;
  z03[2]=4.293*cm;
  z03[3]=6.038*cm;
  z03[4]=6.039*cm;
  z03[5]=64.468*cm;
  z03[6]=64.467*cm;
  z03[7]=67.078*cm;

  
  G4double* rInnn03 = new G4double[8];
  rInnn03[0]= 13.970*0.5*cm;
  rInnn03[1]= 13.970*0.5*cm;
  rInnn03[2]= 13.970*0.5*cm;
  rInnn03[3]= 13.970*0.5*cm;
  rInnn03[4]= 13.970*0.5*cm;
  rInnn03[5]= 13.970*0.5*cm;
  rInnn03[6]= 0;
  rInnn03[7]= 0;
  
   

  G4double* rOut03 = new G4double[8];
  rOut03[0]=41.910*0.5*cm;
  rOut03[1]=41.910*0.5*cm;
  rOut03[2]=26.975*0.5*cm;
  rOut03[3]=26.975*0.5*cm; 
  rOut03[4]=16.086*0.5*cm; 
  rOut03[5]=16.086*0.5*cm;
  rOut03[6]=16.086*0.5*cm;
  rOut03[7]=16.086*0.5*cm;
  
  


  G4Polycone* Bellows_weldment = new G4Polycone("Bellows_weldment", 0.*deg, 360.*deg,
  8, z03, rInnn03, rOut03);

  Bellows_weldment_log = new G4LogicalVolume(Bellows_weldment, vessel_mat , "Bellows_weldment_log");
  Bellows_weldment_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,-58.769*cm - 4.293*cm - 1.746*cm), Bellows_weldment_log,"Bellows_weldment_phys", Inside_vacuum_vessel_log, false, 0) ;

  Bellows_weldment_log->SetVisAttributes(ssteal_blue);
 
 
 
 // bellows_cont_log
 








  // Hydraulic fluid                   


  G4double cil2_hz = cil1_hz - 1*mm ;               // 1.3*cm es el espesor de la pared mas externa


  G4double* z2 = new G4double[10];
  z2[0]= -19*cm;
  z2[1]= 2.3*cm;
  z2[2]= 2.301*cm;
  z2[3]= cil2_hz;
  z2[4]= z2[3]+2*r2*cos(25 * PI / 180.0)-sqrt(3)*r2;
  z2[5]= z2[3]+2*r2*cos(20 * PI / 180.0)-sqrt(3)*r2;
  z2[6]= z2[3]+2*r2*cos(15 * PI / 180.0)-sqrt(3)*r2;
  z2[7]= z2[3]+2*r2*cos(10 * PI / 180.0)-sqrt(3)*r2;
  z2[8]= z2[3]+2*r2*cos(5 * PI / 180.0)-sqrt(3)*r2;
  z2[9]= z2[3]+2*r2-sqrt(3)*r2;

  //G4cout<<"Altura P_Vessel Interna\n"<<z2[7]<<G4endl;    // 71.41 cm


  G4double* rOut2 = new G4double[10];
  rOut2[0]=r2;
  rOut2[1]=r2;
  rOut2[2]=r2;
  rOut2[3]=r2;
  rOut2[4]=2*r2*sin(25 * PI / 180.0);
  rOut2[5]=2*r2*sin(20 * PI / 180.0);
  rOut2[6]=2*r2*sin(15 * PI / 180.0);
  rOut2[7]=2*r2*sin(10 * PI / 180.0);
  rOut2[8]=2*r2*sin(5 * PI / 180.0);
  rOut2[9]=0.;
  
  
  
  G4double* rInn2 = new G4double[10];
  rInn2[0]=16.1*0.5*cm;
  rInn2[1]=16.1*0.5*cm;
  rInn2[2]=0;
  rInn2[3]=0;
  rInn2[4]=0;
  rInn2[5]=0;
  rInn2[6]=0;
  rInn2[7]=0;
  rInn2[8]=0;
  rInn2[9]=0;
  
  

  G4Polycone* contenedor2 = new G4Polycone("contenedor2", 0.*deg, 360.*deg,
  10, z2, rInn2, rOut2);

  hydraulic_fluid_log = new G4LogicalVolume(contenedor2, the_CF4, "hydraulic_fluid_log");
  hydraulic_fluid_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), hydraulic_fluid_log, "hydraulic_fluid_phys",pressure_vessel_log, false, 0) ;


  G4VisAttributes* CF4_cyan= new G4VisAttributes(cyan);
  CF4_cyan->SetVisibility(true);
  hydraulic_fluid_log->SetVisAttributes(CF4_cyan);



/*


   // Place to simulate emanation


  G4double* rOut_emanacion = new G4double[2];
  rOut_emanacion[0]=r2;
  rOut_emanacion[1]=r2;


  G4double* rInn_emanacion = new G4double[2];
  rInn_emanacion[0]=r2-0.5*mm;
  rInn_emanacion[1]=r2-0.5*mm;


  G4Polycone* contenedor_emanacion = new G4Polycone("contenedor_emanacion", 0.*deg, 360.*deg,
  2, z2, rInn_emanacion, rOut_emanacion);

  contenedor_emanacion_log = new G4LogicalVolume(contenedor_emanacion, the_CF4, "contenedor_emanacion_log");
  contenedor_emanacion_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), contenedor_emanacion_log, "contenedor_emanacion_phys",
                                                         hydraulic_fluid_log, false, 0) ;

*/

// Camera Port Virtual(cp)
          

//Pressure vessel  Inverted (Virtual Volume)


  G4double* z1_pvi = new G4double[9];
  z1_pvi[0]= 0.;
  z1_pvi[1]= cil1_hz;
  z1_pvi[2]= z1[4]+2*r1*cos(25 * PI / 180.0)-sqrt(3)*r1;
  z1_pvi[3]= z1[4]+2*r1*cos(20 * PI / 180.0)-sqrt(3)*r1;
  z1_pvi[4]= z1[4]+2*r1*cos(15 * PI / 180.0)-sqrt(3)*r1;
  z1_pvi[5]= z1[4]+2*r1*cos(10 * PI / 180.0)-sqrt(3)*r1;
  z1_pvi[6]= z1[4]+2*r1*cos(5 * PI / 180.0)-sqrt(3)*r1;
  z1_pvi[7]= z1[4]+2*r1-sqrt(3)*r1;
  z1_pvi[8]= z1_pvi[7] + 10*cm;


G4double* rInn_pvi = new G4double[9];
  rInn_pvi[0]=r1;
  rInn_pvi[1]=r1;
  rInn_pvi[2]=2*r1*sin(25 * PI / 180.0);
  rInn_pvi[3]=2*r1*sin(20 * PI / 180.0);
  rInn_pvi[4]=2*r1*sin(15 * PI / 180.0);
  rInn_pvi[5]=2*r1*sin(10 * PI / 180.0);
  rInn_pvi[6]=2*r1*sin(5 * PI / 180.0);
  rInn_pvi[7]=0.;
  rInn_pvi[8]=0.;


  G4double* rOut_pvi = new G4double[9];
  rOut_pvi[0]=30*cm;
  rOut_pvi[1]=30*cm;
  rOut_pvi[2]=30*cm;
  rOut_pvi[3]=30*cm;
  rOut_pvi[4]=30*cm;
  rOut_pvi[5]=30*cm;
  rOut_pvi[6]=30*cm;
  rOut_pvi[7]=30*cm;
  rOut_pvi[8]=30*cm;

  G4Polycone* pvi = new G4Polycone("pvi", 0.*deg, 360.*deg,
  9, z1_pvi, rInn_pvi, rOut_pvi);

 //Camera Port Virtual Volume Rotated


 G4VSolid* Inter_cp_PV = new G4IntersectionSolid("Inter_cp_PV", pvi, Camera_Port, yRot, G4ThreeVector(-3.0*2.54*cm,0.,72*cm-2*2.54*cm));   


 // Camera Port

 Camera_port_log = new G4LogicalVolume(Inter_cp_PV, vessel_mat, "Camera_port_log");
 Camera_port_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), Camera_port_log, "Camera_port_phys", hydraulic_fluid_log, false, 0) ;

Camera_port_log->SetVisAttributes(ssteal_blue);

 Camera_port2_log = new G4LogicalVolume(Inter_cp_PV, vessel_mat, "Camera_port2_log");
 Camera_port2_phys = new G4PVPlacement(zzRot, G4ThreeVector(0.,0.,0.), Camera_port2_log, "Camera_port2_phys", hydraulic_fluid_log, false, 0) ;

Camera_port2_log->SetVisAttributes(ssteal_blue);

 Camera_port3_log = new G4LogicalVolume(Inter_cp_PV, vessel_mat, "Camera_port3_log");
 Camera_port3_phys = new G4PVPlacement(zz2Rot, G4ThreeVector(0.,0.,0.), Camera_port3_log, "Camera_port3_phys", hydraulic_fluid_log, false, 0) ;

Camera_port3_log->SetVisAttributes(ssteal_blue);



//Camera System  Nanoguide

 G4double camera_system_hz = 20;
 G4double camera_system_hx = 10;
 G4double mov = 20.55;

 G4Box* Camera_System_box = new G4Box
     ("Camera_System_box", 0.5*camera_system_hx*cm, 0.5*camera_system_hx*cm, 0.5*camera_system_hz*cm );
     
 Camera_System_log = new G4LogicalVolume(Camera_System_box, vacuum_mat, "Camera_System_log");
 Camera_System_phys = new G4PVPlacement(yRot, G4ThreeVector(-3.0*2.54*cm - 0.3827*mov*cm,0.,72*cm-2*2.54*cm + 0.9239*mov*cm), Camera_System_log, "Camera_System_phys", Inside_vacuum_vessel_log, false, 0) ;    
 
 
 //Sapphire Viewport
 
 
 G4Tubs* Sapphire_Cil = new G4Tubs("Sapphire_Cil",0, 0.5*3.175*cm, 0.5*1.67*cm, 0.*deg, 360.*deg);
 
 Sapphire_log = new G4LogicalVolume(Sapphire_Cil, sapphire_mat, "Sapphire_log");
 Sapphire_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,-camera_system_hz/2*cm + 1.67/2*cm), Sapphire_log, "Sapphire_phys", Camera_System_log, false, 0) ; 
 
  G4VisAttributes* Sapphire_red= new G4VisAttributes(red);
  Sapphire_red->SetVisibility(true);
  Sapphire_log->SetVisAttributes(Sapphire_red); 
 

 G4Tubs* Sapphire_Ssteal_Cil = new G4Tubs("Sapphire_Ssteal_Cil", 0.5*3.175*cm, 0.5*8.56*cm, 0.5*1.67*cm, 0.*deg, 360.*deg);
 
 Sapphire_Ssteal_log = new G4LogicalVolume(Sapphire_Ssteal_Cil, vessel_mat, "Sapphire_Ssteal_log");
 Sapphire_Ssteal_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,-camera_system_hz/2*cm + 1.67/2*cm), Sapphire_Ssteal_log, "Sapphire_Ssteal_phys", Camera_System_log, false, 0) ;
 
 Sapphire_Ssteal_log->SetVisAttributes(ssteal_blue);


 //Camera Lens Holder
 
 
 G4Tubs* Lens_Cil = new G4Tubs("Lens_Cil",0, 0.5*1.2*cm, 0.5*1.55*cm, 0.*deg, 360.*deg);
 G4Tubs* Lens_Holder1_Cil = new G4Tubs("Lens_Holder_Cil",0.5*1.2*cm, 0.5*9.7*cm, 0.5*1.55*cm, 0.*deg, 360.*deg);
 
 G4VSolid* Lens_Holder_Cil= new G4SubtractionSolid("Lens_Holder_Cil", Lens_Holder1_Cil, Lens_Cil, 0, G4ThreeVector(0,0.,0.));
 
 Lens_Holder_log = new G4LogicalVolume(Lens_Holder_Cil, vessel_mat, "Lens_Holder_log");
 Lens_Holder_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,-camera_system_hz/2*cm + 1.67*cm + 1.55/2*cm), Lens_Holder_log, "Lens_Holder_phys", Camera_System_log, false, 0) ;

 Lens_Holder_log->SetVisAttributes(ssteal_blue);
 
 
 
 
 Lens_log = new G4LogicalVolume(Lens_Cil, panel_mat, "Lens_log");
 Lens_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,-camera_system_hz/2*cm + 1.67*cm + 1.55/2*cm), Lens_log, "Lens_phys", Camera_System_log, false, 0) ;
 
 
 //Alignment Ring
 
 
 G4Tubs* Alignment_Ring1_Cil = new G4Tubs("Alignment_Ring1_Cil",0, 0.5*5*cm, 0.5*2.6*cm, 0.*deg, 360.*deg);
 G4Box* Nanoguide_hole_box = new G4Box
     ("Nanoguide_hole_box", 0.5*0.65*cm, 0.5*0.65*cm, 0.5*16*cm );
 
 
G4VSolid* Alignment_Ring_Cil= new G4SubtractionSolid("Alignment_Ring_Cil", Alignment_Ring1_Cil, Nanoguide_hole_box, 0, G4ThreeVector(0,0.,0.));
 
 Alignment_Ring_log = new G4LogicalVolume(Alignment_Ring_Cil, vessel_mat, "Alignment_Ring_log");
 Alignment_Ring_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,-camera_system_hz/2*cm + 1.67*cm + 1.55*cm + 2.6/2*cm), Alignment_Ring_log, "Alignment_Ring_phys", Camera_System_log, false, 0) ;
 
 
 Alignment_Ring_log->SetVisAttributes(ssteal_blue);
 
 
 // Front Nanoguide Holder
 
 
 G4Box* Front_Nanoguide1_Holder_box = new G4Box
     ("Front_Nanoguide1_Holder_box", 0.5*4*cm, 0.5*1.2*cm, 0.5*2.1*cm );
     
 G4VSolid* Front_Nanoguide_Holder_box= new G4SubtractionSolid("Front_Nanoguide_Holder_box", Front_Nanoguide1_Holder_box, Nanoguide_hole_box, 0, G4ThreeVector(0,0.,0.));    
 
 Front_Nanoguide_Holder_log = new G4LogicalVolume(Front_Nanoguide_Holder_box, vessel_mat, "Front_Nanoguide_Holder_log");
 Front_Nanoguide_Holder_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,-camera_system_hz/2*cm + 1.67*cm + 1.55*cm + 2.6*cm + 0.5*2.1*cm), Front_Nanoguide_Holder_log, "Front_Nanoguide_Holder_phys", Camera_System_log, false, 0) ;
 
  Front_Nanoguide_Holder_log->SetVisAttributes(ssteal_blue);
  
  
 //Peaks Rod
 
 
 G4Tubs* Peak_Rod_Cil = new G4Tubs("Peak_Rod_Cil",0, 0.5*0.635*cm, 0.5*8*cm, 0.*deg, 360.*deg);
 
 Peak_Rod1_log = new G4LogicalVolume(Peak_Rod_Cil, vessel_mat, "Peak_Rod1_log");
 Peak_Rod1_phys = new G4PVPlacement(0, G4ThreeVector(2*cm - 0.5*0.635*cm,0.,-camera_system_hz/2*cm + 1.67*cm + 1.55*cm + 2.6*cm + 2.1*cm + 0.5*8*cm), Peak_Rod1_log, "Peak_Rod1_phys", Camera_System_log, false, 0) ;
 
 Peak_Rod1_log->SetVisAttributes(ssteal_blue);
 
 
 Peak_Rod2_log = new G4LogicalVolume(Peak_Rod_Cil, vessel_mat, "Peak_Rod2_log");
 Peak_Rod2_phys = new G4PVPlacement(0, G4ThreeVector(-2*cm + 0.5*0.635*cm,0.,-camera_system_hz/2*cm + 1.67*cm + 1.55*cm + 2.6*cm + 2.1*cm + 0.5*8*cm), Peak_Rod2_log, "Peak_Rod2_phys", Camera_System_log, false, 0) ;
 
 Peak_Rod2_log->SetVisAttributes(ssteal_blue);
 
 
 
 //Camera OV9281 - 1MP
 
 G4Tubs* Camera_Cil = new G4Tubs("Camera_Cil",0.5*1*cm, 0.5*1.2*cm, 0.5*0.47*cm, 0.*deg, 360.*deg);
 
 
 Camera_log = new G4LogicalVolume(Camera_Cil, panel_mat, "Camera_log");        //Al
 Camera_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,-camera_system_hz/2*cm + 1.67*cm + 1.55*cm + 2.6*cm + 2.1*cm + 8*cm + 3.5*cm + 0.12*cm - 0.5*0.47*cm), Camera_log, "Camera_phys", Camera_System_log, false, 0) ;
 
 
 
 //Rear nanoguide holder
 
 
 
 G4Box* Rear_Nanoguide1_Holder_box = new G4Box
     ("Rear_Nanoguide1_Holder_box", 0.5*4*cm, 0.5*2*cm, 0.5*3.5*cm );
     
  G4VSolid* Rear_Nanoguide_Holder_box= new G4SubtractionSolid("Rear_Nanoguide_Holder_box", Rear_Nanoguide1_Holder_box, Nanoguide_hole_box, 0, G4ThreeVector(0,0.,0.));  
  
  G4VSolid* Rear_Nanoguide2_Holder_box= new G4SubtractionSolid("Rear_Nanoguide2_Holder_box", Rear_Nanoguide_Holder_box, Camera_Cil, 0, G4ThreeVector(0,0.,0.5*3.5*cm - 0.5*0.47*cm + 0.12*cm));     
     
 Rear_Nanoguide_Holder_log = new G4LogicalVolume(Rear_Nanoguide2_Holder_box, vessel_mat, "Rear_Nanoguide_Holder_log");
 Rear_Nanoguide_Holder_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,-camera_system_hz/2*cm + 1.67*cm + 1.55*cm + 2.6*cm + 2.1*cm + 8*cm + 0.5*3.5*cm), Rear_Nanoguide_Holder_log, "Rear_Nanoguide_Holder_phys", Camera_System_log, false, 0) ;
 
 Rear_Nanoguide_Holder_log->SetVisAttributes(ssteal_blue);
 
 
 //Nanoguide
 
 
 Nanoguide_log = new G4LogicalVolume(Nanoguide_hole_box, Photopolymer_mat, "Nanoguide_log");
 Nanoguide_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,-camera_system_hz/2*cm + 0.5*16*cm + 1.55*cm + 1.67*cm), Nanoguide_log, "Nanoguide_phys", Camera_System_log, false, 0) ;  
 
 
 
 
  G4VisAttributes* nanoguide_green= new G4VisAttributes(green);
  nanoguide_green->SetVisibility(true);
  Nanoguide_log->SetVisAttributes(nanoguide_green);
 
 //Sensor Plate
 
  
 G4Box* Sensor_Plate1_box = new G4Box
     ("Sensor_Plate1_box", 0.5*4*cm, 0.5*2*cm, 0.5*0.12*cm );
  
  G4VSolid* Sensor_Plate_box= new G4SubtractionSolid("Sensor_Plate_Holder_box", Sensor_Plate1_box, Nanoguide_hole_box, 0, G4ThreeVector(0,0.,0.));  
  
  G4VSolid* Sensor_Plate2_box= new G4SubtractionSolid("Sensor_Plate2_box", Sensor_Plate_box, Camera_Cil, 0, G4ThreeVector(0,0.,0.));
 

 Sensor_Plate_log = new G4LogicalVolume(Sensor_Plate2_box, CuShield_mat, "Sensor_Plate_log");
 Sensor_Plate_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,-camera_system_hz/2*cm + 1.67*cm + 1.55*cm + 2.6*cm + 2.1*cm + 8*cm + 3.5*cm + 0.5*0.12*cm), Sensor_Plate_log, "Sensor_Plate_phys", Camera_System_log, false, 0) ;
 


  G4VisAttributes* Cu_orange= new G4VisAttributes(orange);
  Cu_orange->SetVisibility(true);
  Sensor_Plate_log->SetVisAttributes(Cu_orange);
  
  
  
  
 //Sensor Cover
 
  G4Box* Sensor_Cover_box = new G4Box
     ("Sensor_Cover_box", 0.5*4*cm, 0.5*2*cm, 0.5*0.3*cm );  
  
  
 Sensor_Cover_log = new G4LogicalVolume(Sensor_Cover_box, CuShield_mat, "Sensor_Cover_log");
 Sensor_Cover_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,-camera_system_hz/2*cm + 1.67*cm + 1.55*cm + 2.6*cm + 2.1*cm + 8*cm + 3.5*cm + 0.12*cm + 0.5*0.3*cm), Sensor_Cover_log, "Sensor_Covere_phys", Camera_System_log, false, 0) ; 
  
 Sensor_Cover_log->SetVisAttributes(Cu_orange); 
 
 
 
 //Camera PCB 
 
 G4Box* Camera_PCB_box = new G4Box
     ("Camera_PCB_box", 0.5*1.45*cm, 0.5*1.45*cm, 0.5*0.058*cm );
     
 Camera_PCB_log = new G4LogicalVolume(Camera_PCB_box,  PBC_mat, "Camera_PCB_log");        //PCB
 Camera_PCB_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,-0.5*0.3*cm + 0.5*0.058*cm), Camera_PCB_log, "Camera_PCB_phys", Sensor_Cover_log, false, 0) ;    
 
  
  
  
 //Screws Camera System
 
 
 G4Box* Screw_box = new G4Box
     ("Screw_box", 0.5*1.1*cm, 0.5*0.665*cm, 0.5*0.665*cm );
     
     
 Screw1_log = new G4LogicalVolume(Screw_box, vessel_mat, "Screw1_log");
 Screw1_phys = new G4PVPlacement(0, G4ThreeVector(2*cm + 0.5*1.1*cm,0.,-camera_system_hz/2*cm + 1.67*cm + 1.55*cm + 2.6*cm + 2.1*cm - 0.5*0.665*cm), Screw1_log, "Screw1_phys", Camera_System_log, false, 0) ;
 
 Screw1_log->SetVisAttributes(ssteal_blue);
 
  
 Screw2_log = new G4LogicalVolume(Screw_box, vessel_mat, "Screw2_log");
 Screw2_phys = new G4PVPlacement(0, G4ThreeVector(-2*cm - 0.5*1.1*cm,0.,-camera_system_hz/2*cm + 1.67*cm + 1.55*cm + 2.6*cm + 2.1*cm - 0.5*0.665*cm), Screw2_log, "Screw2_phys", Camera_System_log, false, 0) ;
 
 Screw2_log->SetVisAttributes(ssteal_blue); 
  
  
 
 Screw3_log = new G4LogicalVolume(Screw_box, vessel_mat, "Screw3_log");
 Screw3_phys = new G4PVPlacement(0, G4ThreeVector(2*cm + 0.5*1.1*cm,0.,-camera_system_hz/2*cm + 1.67*cm + 1.55*cm + 2.6*cm + 2.1*cm + 0.5*0.665*cm + 8*cm), Screw3_log, "Screw3_phys", Camera_System_log, false, 0) ;
 
 Screw3_log->SetVisAttributes(ssteal_blue); 
 
 
 Screw4_log = new G4LogicalVolume(Screw_box, vessel_mat, "Screw4_log");
 Screw4_phys = new G4PVPlacement(0, G4ThreeVector(-2*cm - 0.5*1.1*cm,0.,-camera_system_hz/2*cm + 1.67*cm + 1.55*cm + 2.6*cm + 2.1*cm + 0.5*0.665*cm + 8*cm), Screw4_log, "Screw4_phys", Camera_System_log, false, 0) ;
 
 Screw4_log->SetVisAttributes(ssteal_blue); 
 
 
  //Springs Camera System
 
 
 G4Tubs* Camera_springs = new G4Tubs("Camera_springs",0.5*0.535*cm, 0.5*0.635*cm, 0.5*8*cm, 0.*deg, 360.*deg);
 
 Camera_spring1_log = new G4LogicalVolume(Camera_springs, vessel_mat, "Camera_spring1_log");
 Camera_spring1_phys = new G4PVPlacement(0, G4ThreeVector(2*cm + 0.5*0.635*cm + 0.2*cm,0.,-camera_system_hz/2*cm + 1.67*cm + 1.55*cm + 2.6*cm + 2.1*cm + 0.5*8*cm), Camera_spring1_log, "Camera_spring1_phys", Camera_System_log, false, 0) ;
 
 Camera_spring1_log->SetVisAttributes(ssteal_blue);
 
 Camera_spring2_log = new G4LogicalVolume(Camera_springs, vessel_mat, "Camera_spring2_log");
 Camera_spring2_phys = new G4PVPlacement(0, G4ThreeVector(-2*cm - 0.5*0.635*cm - 0.2*cm,0.,-camera_system_hz/2*cm + 1.67*cm + 1.55*cm + 2.6*cm + 2.1*cm + 0.5*8*cm), Camera_spring2_log, "Camera_spring2_phys", Camera_System_log, false, 0) ;
 
 Camera_spring2_log->SetVisAttributes(ssteal_blue);
 
 
 
 //Camera System  Relay Lens
 

 G4double camera_system1_hx = 10.5;
 G4double camera_system1_hz = 22;
 G4double mov1 = 21.55;

 G4Box* Camera_System1_box = new G4Box
     ("Camera_System_box", 0.5*camera_system1_hx*cm, 0.5*camera_system1_hx*cm, 0.5*camera_system1_hz*cm );
 
 G4double cos30 = 0.86625;
 G4double sin30 = 0.5;
     
 Camera_System1_log = new G4LogicalVolume(Camera_System1_box, vacuum_mat, "Camera_System1_log");
 Camera_System1_phys = new G4PVPlacement(xyRot, G4ThreeVector(3.0*2.54*sin30*cm + 0.3827*mov1*sin30*cm, 3.0*2.54*cos30*cm + 0.3827*mov1*cos30*cm,72*cm-2*2.54*cm + 0.9239*mov1*cm), Camera_System1_log, "Camera_System1_phys", Inside_vacuum_vessel_log, false, 0) ;    
 

 
 //Sapphire Viewport 1
 
 
 
 
 G4Tubs* Sapphire1_Cil = new G4Tubs("Sapphire_Cil",0, 0.5*3.175*cm, 0.5*1.6*cm, 0.*deg, 360.*deg);
 
 Sapphire1_log = new G4LogicalVolume(Sapphire1_Cil, sapphire_mat, "Sapphire1_log");
 Sapphire1_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,-camera_system1_hz/2*cm + 1.60/2*cm), Sapphire1_log, "Sapphire1_phys", Camera_System1_log, false, 0) ; 
 
 Sapphire1_log->SetVisAttributes(Sapphire_red); 
  
  
 G4Tubs* Sapphire1_Ssteal_Cil = new G4Tubs("Sapphire1_Ssteal_Cil", 0.5*3.175*cm, 0.5*8.56*cm, 0.5*1.60*cm, 0.*deg, 360.*deg);
 
 Sapphire1_Ssteal_log = new G4LogicalVolume(Sapphire1_Ssteal_Cil, vessel_mat, "Sapphire1_Ssteal_log");
 Sapphire1_Ssteal_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,-camera_system1_hz/2*cm + 1.60/2*cm), Sapphire1_Ssteal_log, "Sapphire1_Ssteal_phys", Camera_System1_log, false, 0) ;
 
 Sapphire1_Ssteal_log->SetVisAttributes(ssteal_blue); 
  
 // Plastic Flange
 
 
 G4double* zp = new G4double[6];
  zp[0]= 0*cm;
  zp[1]=0.399*cm;
  zp[2]= 0.4*cm;
  zp[3]= 2.399*cm;
  zp[4]= 2.4*cm;
  zp[5]= 2.8*cm;


  G4double* rInnp = new G4double[6];
  rInnp[0]=0.5*9.7*cm;
  rInnp[1]=0.5*9.7*cm;
  rInnp[2]=0.5*4*cm;
  rInnp[3]=0.5*4*cm;
  rInnp[4]=0.5*9.7*cm;
  rInnp[5]=0.5*9.7*cm;


  G4double* rOutp = new G4double[6];
  rOutp[0]=0.5*10.5*cm;
  rOutp[1]=0.5*10.5*cm;
  rOutp[2]=0.5*10.5*cm;
  rOutp[3]=0.5*10.5*cm;
  rOutp[4]=0.5*10.5*cm;
  rOutp[5]=0.5*10.5*cm;




  G4Polycone* Plastic_Flange_Cil = new G4Polycone("Plastic_Flange_Cil", 0.*deg, 360.*deg,
  6, zp, rInnp, rOutp); 
  
  Plastic_Flange_log = new G4LogicalVolume(Plastic_Flange_Cil, sapphire_mat, "Plastic_Flange_log"); //Delrin_mat
  Plastic_Flange_phys = new G4PVPlacement(0, G4ThreeVector(0.,0., -camera_system1_hz/2*cm + 1.60*cm), Plastic_Flange_log, "Plastic_Flange_phys", Camera_System1_log, false, 0) ;
  
  G4VisAttributes* Plastic_brown= new G4VisAttributes(brown);
  Plastic_brown->SetVisibility(true);
  Plastic_Flange_log->SetVisAttributes(Plastic_brown);
 
 
  // Camera Lens Holder
 
  G4double* zl = new G4double[6];
  zl[0]= 0*cm;
  zl[1]= 1.99*cm;
  zl[2]= 2.0*cm;
  zl[3]= 2.399*cm; 
  zl[4]= 2.4*cm;
  zl[5]= 4.8*cm;
 


  G4double* rInnl = new G4double[6];
  rInnl[0]=0.5*1.2*cm;
  rInnl[1]=0.5*1.2*cm;
  rInnl[2]=0.5*1.2*cm;
  rInnl[3]=0.5*1.2*cm;
  rInnl[4]=0.5*4.8*cm;
  rInnl[5]=0.5*4.8*cm;
 


  G4double* rOutl = new G4double[6];
  rOutl[0]=0.5*4*cm;
  rOutl[1]=0.5*4*cm;
  rOutl[2]=0.5*9.7*cm;
  rOutl[3]=0.5*9.7*cm;
  rOutl[4]=0.5*5.6*cm;
  rOutl[5]=0.5*5.6*cm;

  
  G4Polycone* Lens_Holder2_Cil = new G4Polycone("Lens_Holder2_Cil", 0.*deg, 360.*deg,
  6, zl, rInnl, rOutl);
  
  Lens_Holder2_log = new G4LogicalVolume(Lens_Holder2_Cil, vessel_mat, "Lens_Holder2_log");
  Lens_Holder2_phys = new G4PVPlacement(0, G4ThreeVector(0.,0., -camera_system1_hz/2*cm + 1.60*cm + 0.4*cm), Lens_Holder2_log, "Lens_Holder2_phys", Camera_System1_log, false, 0) ;
  
 Lens_Holder2_log->SetVisAttributes(ssteal_blue);  
  
 //Lens
 
 Lens1_log = new G4LogicalVolume(Lens_Cil, panel_mat, "Lens1_log");
 Lens1_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,-camera_system1_hz/2*cm + 1.60*cm + 0.4*cm + 1.55*0.5*cm), Lens1_log, "Lens1_phys", Camera_System1_log, false, 0) ;
 
 
 //Adjustment Ring 
 
 
  G4double* za = new G4double[4];
  za[0]= 0*cm;
  za[1]= 3.99*cm;
  za[2]= 4.0*cm;
  za[3]= 4.3*cm; 

  G4double* rInna = new G4double[4];
  rInna[0]=0.5*4*cm;
  rInna[1]=0.5*4*cm;
  rInna[2]=0.5*4*cm;
  rInna[3]=0.5*4*cm;
 
  G4double* rOuta = new G4double[4];
  rOuta[0]=0.5*4.8*cm;
  rOuta[1]=0.5*4.8*cm;
  rOuta[2]=0.5*5.6*cm;
  rOuta[3]=0.5*5.6*cm;
  
  
   G4Polycone* Adjustment_Cil = new G4Polycone("Adjustment_Cil", 0.*deg, 360.*deg,
  4, za, rInna, rOuta);
  
  Adjustment_log = new G4LogicalVolume(Adjustment_Cil, vessel_mat, "Adjustment_log");
  Adjustment_phys = new G4PVPlacement(0, G4ThreeVector(0.,0., -camera_system1_hz/2*cm + 2.8*cm + 1.6*cm), Adjustment_log, "Adjustment_phys", Camera_System1_log, false, 0) ;
  
   Adjustment_log->SetVisAttributes(Cu_orange);
   
   
   
  //Top Plastic Ring 
  
  G4double* zt = new G4double[6];
  zt[0]= 0*cm;
  zt[1]= 0.099*cm;
  zt[2]= 0.1*cm;
  zt[3]= 1.199*cm; 
  zt[4]= 1.2*cm;
  zt[5]= 1.5*cm;


  G4double* rInnt = new G4double[6];
  rInnt[0]=0.5*4*cm;
  rInnt[1]=0.5*4*cm;
  rInnt[2]=0.5*5*cm;
  rInnt[3]=0.5*5*cm;
  rInnt[4]=0.5*5*cm;
  rInnt[5]=0.5*5*cm;

 
  G4double* rOutt = new G4double[6];
  rOutt[0]=0.5*5.6*cm;
  rOutt[1]=0.5*5.6*cm;
  rOutt[2]=0.5*5.6*cm;
  rOutt[3]=0.5*5.6*cm; 
  rOutt[4]=0.5*6*cm; 
  rOutt[5]=0.5*6*cm; 
 
   G4Polycone* Top_Plastic_Cil = new G4Polycone("Top_Plastic_Cil", 0.*deg, 360.*deg,
  6, zt, rInnt, rOutt);
  
  Top_Plastic_log = new G4LogicalVolume(Top_Plastic_Cil, PEEK_mat, "Top_Plastic_log");
  Top_Plastic_phys = new G4PVPlacement(0, G4ThreeVector(0.,0., -camera_system1_hz/2*cm + 2.8*cm + 1.6*cm + 4.3*cm), Top_Plastic_log, "Top_Plastic_phys", Camera_System1_log, false, 0) ;
  
  Top_Plastic_log->SetVisAttributes(Plastic_brown);
  
  
  
  //Aspheric lens
  
  G4double ral = 4.85*cm;
  
  G4double* zal = new G4double[8];
  zal[0]= 0*cm;
  zal[1]= 0.442*cm;
  zal[2]= zal[1]+2*ral*cos(25 * PI / 180.0)-sqrt(3)*ral;
  zal[3]= zal[1]+2*ral*cos(20 * PI / 180.0)-sqrt(3)*ral;
  zal[4]= zal[1]+2*ral*cos(15 * PI / 180.0)-sqrt(3)*ral;
  zal[5]= zal[1]+2*ral*cos(10 * PI / 180.0)-sqrt(3)*ral;
  zal[6]= zal[1]+2*ral*cos(5 * PI / 180.0)-sqrt(3)*ral;
  zal[7]= zal[1]+2*ral-sqrt(3)*ral;
  
  

  G4double* rInnal = new G4double[8];
  rInnal[0]=0;
  rInnal[1]=0; 
  rInnal[2]=0;
  rInnal[3]=0;
  rInnal[4]=0;
  rInnal[5]=0;
  rInnal[6]=0;
  rInnal[7]=0;

 
  G4double* rOutal = new G4double[8];
  rOutal[0]=0.5*5*cm;
  rOutal[1]=0.5*5*cm;
  rOutal[2]=ral*sin(25 * PI / 180.0);
  rOutal[3]=ral*sin(20 * PI / 180.0);
  rOutal[4]=ral*sin(15 * PI / 180.0);
  rOutal[5]=ral*sin(10 * PI / 180.0);
  rOutal[6]=ral*sin(5 * PI / 180.0);
  rOutal[7]=0.;
  
  
  G4Polycone* Aspheric_lens_Cil = new G4Polycone("Aspheric_lens_Cil", 0.*deg, 360.*deg,
  8, zal, rInnal, rOutal);
  
  
  Aspheric_lens_log = new G4LogicalVolume(Aspheric_lens_Cil, glass_mat, "Aspheric_lens_log");
  Aspheric_lens_phys = new G4PVPlacement(0, G4ThreeVector(0.,0., -camera_system1_hz/2*cm + 2.8*cm + 1.6*cm + 4.3*cm + 0.1*cm), Aspheric_lens_log, "Aspheric_lens_phys", Camera_System1_log, false, 0) ;
  
  
 //Lens Holder
 
 
  G4double* zh = new G4double[6];
  zh[0]= 0*cm;
  zh[1]= 0.299*cm;
  zh[2]= 0.3*cm;
  zh[3]= 1.699*cm;
  zh[4]= 1.7*cm;
  zh[5]= 1.9*cm;
  
  G4double* rInnh = new G4double[6];
  rInnh[0]=0.5*4.8*cm;
  rInnh[1]=0.5*4.8*cm; 
  rInnh[2]=0.5*5.6*cm;
  rInnh[3]=0.5*5.6*cm;
  rInnh[4]=0.5*6*cm;
  rInnh[5]=0.5*6*cm;

  G4double* rOuth = new G4double[6];
  rOuth[0]=0.5*6.7*cm;
  rOuth[1]=0.5*6.7*cm;
  rOuth[2]=0.5*6.7*cm;
  rOuth[3]=0.5*6.7*cm;
  rOuth[4]=0.5*6.7*cm;
  rOuth[5]=0.5*6.7*cm;
 
 
  G4Polycone* Lens_Holder3_Cil = new G4Polycone("Lens_Holder3_Cil", 0.*deg, 360.*deg,
  6, zh, rInnh, rOuth);
  
  G4Tubs* Lens_Holder_tab = new G4Tubs("Tube_SS",0.8*0.5*cm, 2.0*0.5*cm, 0.7*0.5*cm, 0.*deg, 360.*deg);
  
  G4VSolid* tab_ring1= new G4UnionSolid("tab_ring1", Lens_Holder3_Cil, Lens_Holder_tab, 0, G4ThreeVector(6.7*0.5*cm + 0.85*cm,0.,1.9*cm-0.7*0.5*cm)); 
  G4VSolid* tab_ring2= new G4UnionSolid("tab_ring2", tab_ring1, Lens_Holder_tab, 0, G4ThreeVector(-6.7*0.5*cm - 0.85*cm,0.,1.9*cm-0.7*0.5*cm));
  G4VSolid* tab_ring3= new G4UnionSolid("tab_ring3", tab_ring2, Lens_Holder_tab, 0, G4ThreeVector(0.,-6.7*0.5*cm - 0.85*cm,1.9*cm-0.7*0.5*cm));
  G4VSolid* tab_ring4= new G4UnionSolid("tab_ring4", tab_ring3, Lens_Holder_tab, 0, G4ThreeVector(0.,6.7*0.5*cm + 0.85*cm,1.9*cm-0.7*0.5*cm));
  
  
  Lens_Holder3_log = new G4LogicalVolume(tab_ring4, vessel_mat, "Lens_Holder3_log");
  Lens_Holder3_phys = new G4PVPlacement(0, G4ThreeVector(0.,0., -camera_system1_hz/2*cm + 2.8*cm + 1.6*cm + 4.3*cm -0.6*cm), Lens_Holder3_log, "Lens_Holder3_phys", Camera_System1_log, false, 0) ;
  
  Lens_Holder3_log->SetVisAttributes(ssteal_blue); 
  
  
  //Iris Holder
  
  G4double* zi = new G4double[4];
  zi[0]= 0*cm;
  zi[1]= 0.199*cm;
  zi[2]= 0.2*cm;  
  zi[3]= 1.6*cm;

  
  G4double* rInni = new G4double[4];
  rInni[0]=0.5*5.1*cm;
  rInni[1]=0.5*5.1*cm; 
  rInni[2]=0.5*6.2*cm;
  rInni[3]=0.5*6.2*cm;

  G4double* rOuti = new G4double[4];
  rOuti[0]=0.5*7*cm;
  rOuti[1]=0.5*7*cm;
  rOuti[2]=0.5*7*cm;
  rOuti[3]=0.5*7*cm;
  
  
  G4Polycone* Iris_Holder_Cil = new G4Polycone("Iris_Holder_Cil", 0.*deg, 360.*deg,
  4, zi, rInni, rOuti);
  
  G4VSolid* tab_iris1= new G4UnionSolid("tab_iris1", Iris_Holder_Cil, Lens_Holder_tab, 0, G4ThreeVector(6.7*0.5*cm + 0.85*cm,0.,0.7*0.5*cm)); 
  G4VSolid* tab_iris2= new G4UnionSolid("tab_iris2", tab_iris1, Lens_Holder_tab, 0, G4ThreeVector(-6.7*0.5*cm - 0.85*cm,0.,0.7*0.5*cm)); 
  
   Iris_Holder_log = new G4LogicalVolume(tab_iris2, vessel_mat, "Iris_Holder_log");
   Iris_Holder_phys = new G4PVPlacement(0, G4ThreeVector(0.,0., -camera_system1_hz/2*cm + 2.8*cm + 1.6*cm + 4.3*cm + 1.5*cm),  Iris_Holder_log, " Iris_Holder_phys", Camera_System1_log, false, 0) ;
  
  G4VisAttributes* ssteal_magenta= new G4VisAttributes(magenta);
  ssteal_magenta->SetVisibility(true);
 Iris_Holder_log->SetVisAttributes(ssteal_magenta);
 
 
 //Iris
 
 
  G4double* zs = new G4double[4];
  zs[0]= 0*cm;
  zs[1]= 1.3699*cm;
  zs[2]= 1.37*cm;   
  zs[3]= 2.71*cm;

  
  G4double* rInns = new G4double[4];
  rInns[0]=0.5*5.4*cm;
  rInns[1]=0.5*5.4*cm; 
  rInns[2]=0.5*3.7*cm; 
  rInns[3]=0.5*3.7*cm;

  G4double* rOuts = new G4double[4];
  rOuts[0]=0.5*6.2*cm;
  rOuts[1]=0.5*6.2*cm;
  rOuts[2]=0.5*6.2*cm;
  rOuts[3]=0.5*6.2*cm;
  
  G4Polycone* Iris_Cil = new G4Polycone("Iris_Cil", 0.*deg, 360.*deg,
  4, zs, rInns, rOuts);
  
  Iris_log = new G4LogicalVolume(Iris_Cil, PEEK_mat, "Iris_log");
  Iris_phys = new G4PVPlacement(0, G4ThreeVector(0.,0., -camera_system1_hz/2*cm + 2.8*cm + 1.6*cm + 4.3*cm + 1.5*cm + 0.86*cm ),  Iris_log, " Iris_phys", Camera_System1_log, false, 0) ;
  
  
  //Iris Brass Housing
  
  
  G4Tubs* Iris_brass = new G4Tubs("Iris_brass",0.5*3.7*cm, 0.5*5.3*cm, 0.5*0.6*cm, 0.*deg, 360.*deg);
  
  Iris_brass_log = new G4LogicalVolume(Iris_brass, brass_mat, "Iris_brass_log");              //Material Cambiar
 Iris_brass_phys = new G4PVPlacement(0, G4ThreeVector(0.,0., -camera_system1_hz/2*cm + 2.8*cm + 1.6*cm + 4.3*cm + 1.5*cm + 0.86*cm + 1.06*cm),  Iris_brass_log, "Iris_brass_phys", Camera_System1_log, false, 0) ;
  
   Iris_brass_log->SetVisAttributes(Sapphire_red);
   
   
   
   //Iris Blade 
  
  G4Tubs* Iris_blades = new G4Tubs("Iris_blades",0, 0.5*3.7*cm, 0.5*0.06*cm, 0.*deg, 360.*deg);
  
  Iris_blades_log = new G4LogicalVolume(Iris_blades, vessel_mat, "Iris_blades_log");              
  Iris_blades_phys = new G4PVPlacement(0, G4ThreeVector(0.,0., -camera_system1_hz/2*cm + 2.8*cm + 1.6*cm + 4.3*cm + 1.5*cm + 0.86*cm + 1.06*cm),  Iris_blades_log, "Iris_blades_phys", Camera_System1_log, false, 0) ;
  
   Iris_blades_log->SetVisAttributes(Sapphire_red);
  
  
  
  //Second Iris Holder
  
  
   Iris_Holder1_log = new G4LogicalVolume(tab_iris2, vessel_mat, "Iris_Holder1_log");
   Iris_Holder1_phys = new G4PVPlacement(y180Rot, G4ThreeVector(0.,0., -camera_system1_hz/2*cm + 2.8*cm + 1.6*cm + 4.3*cm + 1.5*cm + 2*1.6*cm + 1.23*cm),  Iris_Holder1_log, " Iris_Holder1_phys", Camera_System1_log, false, 0) ;
  
 Iris_Holder1_log->SetVisAttributes(ssteal_magenta);
  
 
 //Second Top Plastic Ring
 
  Top_Plastic1_log = new G4LogicalVolume(Top_Plastic_Cil, PEEK_mat, "Top_Plastic1_log");
  Top_Plastic1_phys = new G4PVPlacement(y180Rot, G4ThreeVector(0.,0., -camera_system1_hz/2*cm + 2.8*cm + 1.6*cm + 4.3*cm + 1.5*cm + 2*1.6*cm + 1.23*cm + 1.5*cm), Top_Plastic1_log, "Top_Plastic1_phys", Camera_System1_log, false, 0) ;
  
  Top_Plastic1_log->SetVisAttributes(Plastic_brown); 
  
 //Second Aspheric Lens  
  
  Aspheric_lens1_log = new G4LogicalVolume(Aspheric_lens_Cil, glass_mat, "Aspheric_lens1_log");
  Aspheric_lens1_phys = new G4PVPlacement(y180Rot, G4ThreeVector(0.,0., -camera_system1_hz/2*cm + 2.8*cm + 1.6*cm + 4.3*cm + 1.5*cm + 2*1.6*cm + 1.23*cm + 1.5*cm - 0.1*cm), Aspheric_lens1_log, "Aspheric_lens1_phys", Camera_System1_log, false, 0) ;
  
  //Second Adjustment Ring
  
  G4double* zaa = new G4double[4];
  zaa[0]= 0*cm;
  zaa[1]= 2.99*cm;
  zaa[2]= 3.0*cm;
  zaa[3]= 3.3*cm; 

  G4double* rInnaa = new G4double[4];
  rInnaa[0]=0.5*4*cm;
  rInnaa[1]=0.5*4*cm;
  rInnaa[2]=0.5*4*cm;
  rInnaa[3]=0.5*4*cm;
 
  G4double* rOutaa = new G4double[4];
  rOutaa[0]=0.5*4.8*cm;
  rOutaa[1]=0.5*4.8*cm;
  rOutaa[2]=0.5*5.6*cm;
  rOutaa[3]=0.5*5.6*cm;
  
  
   G4Polycone* Adjustment1_Cil = new G4Polycone("Adjustment_Cil", 0.*deg, 360.*deg,
  4, zaa, rInnaa, rOutaa);
  
  
  Adjustment1_log = new G4LogicalVolume(Adjustment1_Cil, vessel_mat, "Adjustment1_log");
  Adjustment1_phys = new G4PVPlacement(y180Rot, G4ThreeVector(0.,0., -camera_system1_hz/2*cm + 2.8*cm + 1.6*cm + 4.3*cm + 1.5*cm + 2*1.6*cm + 1.23*cm + 1.5*cm + 3.3*cm), Adjustment1_log, "Adjustment1_phys", Camera_System1_log, false, 0) ;
  
   Adjustment1_log->SetVisAttributes(Cu_orange);
  
  
  //Second Lens Holder
  
  Lens_Holder4_log = new G4LogicalVolume(tab_ring4, vessel_mat, "Lens_Holder4_log");
  Lens_Holder4_phys = new G4PVPlacement(y180Rot, G4ThreeVector(0.,0., -camera_system1_hz/2*cm + 2.8*cm + 1.6*cm + 4.3*cm + 1.5*cm + 2*1.6*cm + 1.23*cm + 1.5*cm + 3.3*cm + 1.9*cm - 4.6*cm), Lens_Holder4_log, "Lens_Holder4_phys", Camera_System1_log, false, 0) ;
  
  Lens_Holder4_log->SetVisAttributes(ssteal_blue); 
  
  
  //Sensor Holder
  
  G4double* zo = new G4double[6];
  zo[0]= 0*cm;
  zo[1]= 1.799*cm;
  zo[2]= 1.8*cm;   
  zo[3]= 2.99*cm; 
  zo[4]= 3.0*cm;
  zo[5]= 3.4*cm;
  
  G4double* rInno = new G4double[6];
  rInno[0]=0.5*4.8*cm;
  rInno[1]=0.5*4.8*cm; 
  rInno[2]=0.5*4.8*cm;
  rInno[3]=0.5*4.8*cm;
  rInno[4]=0;
  rInno[5]=0;  


  G4double* rOuto = new G4double[6];
  rOuto[0]=0.5*6.7*cm;
  rOuto[1]=0.5*6.7*cm;
  rOuto[2]=0.5*5.6*cm;
  rOuto[3]=0.5*5.6*cm;
  rOuto[4]=0.5*5.6*cm;
  rOuto[5]=0.5*5.6*cm;
  
   G4Polycone* Sensor_Holder_Cil = new G4Polycone("Sensor_Holder_Cil", 0.*deg, 360.*deg,
  6, zo, rInno, rOuto);
  
  G4VSolid* tab_sensor1= new G4UnionSolid("tab_sensor1", Sensor_Holder_Cil, Lens_Holder_tab, 0, G4ThreeVector(6.7*0.5*cm + 0.85*cm,0.,-0.7*0.5*cm + 1.8*cm)); 
  G4VSolid* tab_sensor2= new G4UnionSolid("tab_sensor2", tab_sensor1, Lens_Holder_tab, 0, G4ThreeVector(-6.7*0.5*cm - 0.85*cm,0.,-0.7*0.5*cm + 1.8*cm)); 
  
  Sensor_Holder_log = new G4LogicalVolume(tab_sensor2, CuShield_mat, " Sensor_Holder_log");
  Sensor_Holder_phys = new G4PVPlacement(0, G4ThreeVector(0.,0., -camera_system1_hz/2*cm + 2.8*cm + 1.6*cm + 4.3*cm + 1.5*cm + 2*1.6*cm + 1.23*cm + 1.5*cm + 3.3*cm + 1.9*cm - 4.6*cm + 1.2*cm),  Sensor_Holder_log, " Sensor_Holder_phys", Camera_System1_log, false, 0) ;
  
   Sensor_Holder_log->SetVisAttributes(Cu_orange);
   
   
   //Copper Rods
   
  G4Tubs* Copper_Rod_Cil = new G4Tubs("Copper_Rod_Cil",0, 0.5*0.794*cm, 0.5*12.7*cm, 0.*deg, 360.*deg);
   
  Copper_Rod1_log = new G4LogicalVolume(Copper_Rod_Cil, CuShield_mat, "Copper_Rod1_log");
  Copper_Rod1_phys = new G4PVPlacement(0, G4ThreeVector(6.7*0.5*cm + 0.85*cm,0., 3.6*cm),  Copper_Rod1_log, "Copper_Rod1_phys", Camera_System1_log, false, 0) ;
  
  Copper_Rod1_log->SetVisAttributes(Cu_orange);
  
  
  Copper_Rod2_log = new G4LogicalVolume(Copper_Rod_Cil, CuShield_mat, "Copper_Rod2_log");
  Copper_Rod2_phys = new G4PVPlacement(0, G4ThreeVector(-6.7*0.5*cm - 0.85*cm,0., 3.6*cm),  Copper_Rod2_log, "Copper_Rod2_phys", Camera_System1_log, false, 0) ;
  
  Copper_Rod2_log->SetVisAttributes(Cu_orange);
  
  
  // SS Rods
  
  G4Tubs* SS_Rod_Cil = new G4Tubs("SS_Rod_Cil",0, 0.5*0.794*cm, 0.5*6.712*cm, 0.*deg, 360.*deg);
  
  SS_Rod1_log = new G4LogicalVolume(SS_Rod_Cil, vessel_mat, "SS_Rod1_log");
  SS_Rod1_phys = new G4PVPlacement(0, G4ThreeVector(0.,6.7*0.5*cm + 0.85*cm, -camera_system1_hz/2*cm + 0.5*6.712*cm + 1.6*cm + 2.8*cm),  SS_Rod1_log, "SS_Rod1_phys", Camera_System1_log, false, 0) ;
  
    SS_Rod1_log->SetVisAttributes(ssteal_blue); 
    
    
    
  SS_Rod2_log = new G4LogicalVolume(SS_Rod_Cil, vessel_mat, "SS_Rod2_log");
  SS_Rod2_phys = new G4PVPlacement(0, G4ThreeVector(0.,-6.7*0.5*cm - 0.85*cm, -camera_system1_hz/2*cm + 0.5*6.712*cm + 1.6*cm + 2.8*cm),  SS_Rod2_log, "SS_Rod2_phys", Camera_System1_log, false, 0) ;
  
    SS_Rod2_log->SetVisAttributes(ssteal_blue);  
    
    
    
    
    //Camera PCB 
 
     
 Camera_PCB1_log = new G4LogicalVolume(Camera_PCB_box,  PBC_mat, "Camera_PCB1_log");        //PCB
 Camera_PCB1_phys = new G4PVPlacement(0, G4ThreeVector(0.,0., -camera_system1_hz/2*cm + 2.8*cm + 1.6*cm + 4.3*cm + 1.5*cm + 2*1.6*cm + 1.23*cm + 1.5*cm + 3.3*cm + 1.9*cm - 4.6*cm + 1.2*cm - 0.5*0.058*cm +2.99*cm), Camera_PCB1_log, "Camera_PCB1_phys", Camera_System1_log, false, 0) ; 
   
   
    
    
    
    //Camera OV9281 - 1MP
 
 G4Tubs* Camera1_Cil = new G4Tubs("Camera_Cil",0., 0.5*1.2*cm, 0.5*0.47*cm, 0.*deg, 360.*deg);
 
 Camera1_log = new G4LogicalVolume(Camera1_Cil, panel_mat, "Camera1_log");        //Al
 Camera1_phys = new G4PVPlacement(0, G4ThreeVector(0.,0., -camera_system1_hz/2*cm + 2.8*cm + 1.6*cm + 4.3*cm + 1.5*cm + 2*1.6*cm + 1.23*cm + 1.5*cm + 3.3*cm + 1.9*cm - 4.6*cm + 1.2*cm - 0.058*cm +2.99*cm - 0.5*0.47*cm), Camera1_log, "Camera1_phys", Camera_System1_log, false, 0) ;
    
     
     
  //Replicate Relay System Camera
  
 Camera_System1_phys = new G4PVPlacement(xy1Rot, G4ThreeVector(3.0*2.54*sin30*cm + 0.3827*mov1*sin30*cm, -3.0*2.54*cos30*cm - 0.3827*mov1*cos30*cm, 72*cm-2*2.54*cm + 0.9239*mov1*cm), Camera_System1_log, "Camera_System2_phys", Inside_vacuum_vessel_log, false, 0) ;    
     
  
 // Camera_System_phys = new G4PVPlacement(yRot, G4ThreeVector(-3.0*2.54*cm - 0.3827*mov*cm, 0. ,72*cm-2*2.54*cm + 0.9239*mov*cm), Camera_System_log, "Camera_System_phys", Inside_vacuum_vessel_log, false, 0) ;    
  
   
  
  
 
 /*
 //Looking for positions 
 
 G4Tubs* looking_Cil = new G4Tubs("Camera_springs",0, 0.5*17*cm, 0.5*24*cm, 0.*deg, 360.*deg);
 
 
 looking_log = new G4LogicalVolume(looking_Cil, vessel_mat, "looking_log");
 looking_phys = new G4PVPlacement(0, G4ThreeVector(7*cm,14*cm,87*cm), looking_log, "looking_phys", world_log, false, 0) ; 
*/

  // Outer Jar              



  G4double cil3_hz = (60.75-3.21+1.26
)*cm ;
  G4double r3 = 12*cm ;

  G4double* z3 = new G4double[11];
  z3[0]= 18.709*cm + 1.26*cm;            //26.45*cm
  z3[1]= z3[0] + 2.999*cm;
  z3[2]= z3[0] + 3.0*cm;
  z3[3]= z3[0] + 4.0*cm;
  z3[4]= cil3_hz;
  z3[5]= z3[4]+2*r3*cos(25 * PI / 180.0)-sqrt(3)*r3;
  z3[6]= z3[4]+2*r3*cos(20 * PI / 180.0)-sqrt(3)*r3;
  z3[7]= z3[4]+2*r3*cos(15 * PI / 180.0)-sqrt(3)*r3;
  z3[8]= z3[4]+2*r3*cos(10 * PI / 180.0)-sqrt(3)*r3;
  z3[9]= z3[4]+2*r3*cos(5 * PI / 180.0)-sqrt(3)*r3;
  z3[10]= z3[4]+2*r3-sqrt(3)*r3;

  //G4cout<<"Tamano Outer Jar\n"<<z3[10]-z3[0]<<G4endl;

  G4double* rOut3 = new G4double[11];
  rOut3[0]=r3 + 2.5*cm;
  rOut3[1]=r3 + 2.5*cm;
  rOut3[2]=r3 ;
  rOut3[3]=r3;
  rOut3[4]=r3;
  rOut3[5]=2*r3*sin(25 * PI / 180.0);
  rOut3[6]=2*r3*sin(20 * PI / 180.0);
  rOut3[7]=2*r3*sin(15 * PI / 180.0);
  rOut3[8]=2*r3*sin(10 * PI / 180.0);
  rOut3[9]=2*r3*sin(5 * PI / 180.0);
  rOut3[10]=0.;

  G4Polycone* contenedor3 = new G4Polycone("contenedor3", 0.*deg, 360.*deg,
  11, z3, rInn, rOut3);

  outer_jar_log = new G4LogicalVolume(contenedor3, glass_mat, "outer_jar_log");
  outer_jar_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), outer_jar_log, "outer_jar_phys",
                                                         hydraulic_fluid_log, false, 0) ;

  G4VisAttributes* quartz_green= new G4VisAttributes(green);
  quartz_green->SetVisibility(true);
 outer_jar_log->SetVisAttributes(quartz_green);
 
 
 
 // Inner Jar                        


 
  G4double r5 = 10.5*cm ;
  G4double r_IJ_Inn = 10.0*cm ;
  

  G4double* z5 = new G4double[10];
  z5[0]= 0.226*cm;
  z5[1]= z5[0] + 1.999*cm;
  z5[2]= z5[0] + 2.0*cm;
  z5[3]= z5[2] + 34.963*cm; 
  z5[4]= z5[3]+2*r5*cos(25 * PI / 180.0)-sqrt(3)*r5;
  z5[5]= z5[3]+2*r5*cos(20 * PI / 180.0)-sqrt(3)*r5;
  z5[6]= z5[3]+2*r5*cos(15 * PI / 180.0)-sqrt(3)*r5;
  z5[7]= z5[3]+2*r5*cos(10 * PI / 180.0)-sqrt(3)*r5;
  z5[8]= z5[3]+2*r5*cos(5 * PI / 180.0)-sqrt(3)*r5;
  z5[9]= z5[3]+2*r5-sqrt(3)*r5;

  //G4cout<<"Tamano cupula Inner Jar\n"<<z5[7]<<G4endl;

  G4double* rOut5 = new G4double[10];
  rOut5[0]=r5 + 2.5*cm;
  rOut5[1]=r5 + 2.5*cm;
  rOut5[2]=r5;
  rOut5[3]=r5; 
  rOut5[4]=2*r5*sin(25 * PI / 180.0);
  rOut5[5]=2*r5*sin(20 * PI / 180.0);
  rOut5[6]=2*r5*sin(15 * PI / 180.0);
  rOut5[7]=2*r5*sin(10 * PI / 180.0);
  rOut5[8]=2*r5*sin(5 * PI / 180.0);
  rOut5[9]=0.;
  
  
    G4double* rInn5 = new G4double[10];
  rInn5[0]=r_IJ_Inn;
  rInn5[1]=r_IJ_Inn;
  rInn5[2]=r_IJ_Inn;
  rInn5[3]=r_IJ_Inn - 0.1*cm; 
  rInn5[4]=2*r_IJ_Inn*sin(25 * PI / 180.0);
  rInn5[5]=2*r_IJ_Inn*sin(20 * PI / 180.0);
  rInn5[6]=2*r_IJ_Inn*sin(15 * PI / 180.0);
  rInn5[7]=2*r_IJ_Inn*sin(10 * PI / 180.0);
  rInn5[8]=2*r_IJ_Inn*sin(5 * PI / 180.0);
  rInn5[9]=0.;
  
  
  
  

  G4Polycone* contenedor5 = new G4Polycone("contenedor5", 0.*deg, 360.*deg,
  10, z5, rInn5, rOut5);

  inner_jar_log = new G4LogicalVolume(contenedor5, glass_mat, "inner_jar_log");
  inner_jar_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), inner_jar_log, "inner_jar_phys",
                                                        hydraulic_fluid_log, false, 0) ;


  inner_jar_log->SetVisAttributes(quartz_green);



 //Reflector                              
 // Hole for SiPMs Array in Reflector
  
  
  G4RotationMatrix* ZHoleRot1 = new G4RotationMatrix;               
   ZHoleRot1->rotateZ(45*deg);
  G4RotationMatrix* ZHoleRot2 = new G4RotationMatrix;               
   ZHoleRot2->rotateZ(90*deg); 
  G4RotationMatrix* ZHoleRot3 = new G4RotationMatrix;               
   ZHoleRot3->rotateZ(135*deg);  
  G4RotationMatrix* ZHoleRot4 = new G4RotationMatrix;               
   ZHoleRot4->rotateZ(180*deg);    
  G4RotationMatrix* ZHoleRot5 = new G4RotationMatrix;               
   ZHoleRot5->rotateZ(225*deg);   
  G4RotationMatrix* ZHoleRot6 = new G4RotationMatrix;               
   ZHoleRot6->rotateZ(270*deg);  
  G4RotationMatrix* ZHoleRot7 = new G4RotationMatrix;               
   ZHoleRot7->rotateZ(315*deg);         
  

  G4double* z_hole = new G4double[2];
  z_hole[0]= 55.28*cm - 1.1*cm;
  z_hole[1]= 55.28*cm + 1.1*cm;

  
  G4double* rInn_hole = new G4double[2];
  rInn_hole[0]=r3;                                  //r3 = 12 cm
  rInn_hole[1]=r3;
  
  
   G4double* rOut_hole = new G4double[2];
  rOut_hole[0]=r3 + 1*cm;                                  
  rOut_hole[1]=r3 + 1*cm;
  
  G4Polycone* hole_1 = new G4Polycone("hole_PTFE", -5*deg, 10*deg,
  2, z_hole, rInn_hole, rOut_hole); 
  
  
  //Hole For Top SiPMs
  
  G4double* z_hole2 = new G4double[2];
  z_hole2[0]= 62.35*cm - 0.8*cm;
  z_hole2[1]= 62.35*cm + 1.2*cm;
  
  
  G4double* rInn_hole2 = new G4double[2];
  rInn_hole2[0]=12.7*cm;                                  
  rInn_hole2[1]=12.7*cm;
  
  G4double* rOut_hole2 = new G4double[2];
  rOut_hole2[0]=15.8*cm;                                  
  rOut_hole2[1]=15.8*cm;
  
  
    G4Polycone* hole_2 = new G4Polycone("hole_PTFE2", -5*deg, 10*deg,
  2, z_hole2, rInn_hole2, rOut_hole2); 
  
  
  
  

  G4VSolid* Array_hole_12 = new G4UnionSolid("Array_hole_12", hole_1, hole_1, ZHoleRot1, G4ThreeVector(0. , 0., 0.));
  G4VSolid* Array_hole_13 = new G4UnionSolid("Array_hole_13", Array_hole_12, hole_1, ZHoleRot2, G4ThreeVector(0. , 0., 0.));
  G4VSolid* Array_hole_14 = new G4UnionSolid("Array_hole_14", Array_hole_13, hole_1, ZHoleRot3, G4ThreeVector(0. , 0., 0.));  
  G4VSolid* Array_hole_15 = new G4UnionSolid("Array_hole_15", Array_hole_14, hole_1, ZHoleRot4, G4ThreeVector(0. , 0., 0.)); 
  G4VSolid* Array_hole_16 = new G4UnionSolid("Array_hole_16", Array_hole_15, hole_1, ZHoleRot5, G4ThreeVector(0. , 0., 0.)); 
  G4VSolid* Array_hole_17 = new G4UnionSolid("Array_hole_17", Array_hole_16, hole_1, ZHoleRot6, G4ThreeVector(0. , 0., 0.));      
  G4VSolid* Array_hole_18 = new G4UnionSolid("Array_hole_18", Array_hole_17, hole_1, ZHoleRot7, G4ThreeVector(0. , 0., 0.)); 
  G4VSolid* Array_hole_2 = new G4UnionSolid("Array_hole_2", Array_hole_18, Array_hole_18, 0, G4ThreeVector(0. , 0., -4.42*cm)); 
  G4VSolid* Array_hole = new G4UnionSolid("Array_hole", Array_hole_2, Array_hole_2, 0, G4ThreeVector(0. , 0., -2*4.42*cm)); 
  G4VSolid* Array_hole_51 = new G4UnionSolid("Array_hole_51", Array_hole, hole_2, 0, G4ThreeVector(0. , 0., 0.)); 
  G4VSolid* Array_hole_52 = new G4UnionSolid("Array_hole_52", Array_hole_51, hole_2, ZHoleRot1, G4ThreeVector(0. , 0., 0.)); 
  G4VSolid* Array_hole_53 = new G4UnionSolid("Array_hole_53", Array_hole_52, hole_2, ZHoleRot2, G4ThreeVector(0. , 0., 0.)); 
  G4VSolid* Array_hole_54 = new G4UnionSolid("Array_hole_54", Array_hole_53, hole_2, ZHoleRot3, G4ThreeVector(0. , 0., 0.));
  G4VSolid* Array_hole_55 = new G4UnionSolid("Array_hole_55", Array_hole_54, hole_2, ZHoleRot4, G4ThreeVector(0. , 0., 0.));
  G4VSolid* Array_hole_56 = new G4UnionSolid("Array_hole_56", Array_hole_55, hole_2, ZHoleRot5, G4ThreeVector(0. , 0., 0.));
  G4VSolid* Array_hole_57 = new G4UnionSolid("Array_hole_57", Array_hole_56, hole_2, ZHoleRot6, G4ThreeVector(0. , 0., 0.));
  G4VSolid* Array_hole_final = new G4UnionSolid("Array_hole_final", Array_hole_57, hole_2, ZHoleRot7, G4ThreeVector(0. , 0., 0.));


 //Reflector Hole for SiPM´s Holder in Cu
 /////////////////////////////////////////////////////////////////////////
  G4double* z_hole_1 = new G4double[2];
  z_hole_1[0]= 55.28*cm - 1.9*cm;
  z_hole_1[1]= 55.28*cm + 1.9*cm;
 

  G4Polycone* hole_1_1 = new G4Polycone("hole_Cu", -9*deg, 18*deg,
  2, z_hole_1, rInn_hole, rOut_hole); 
  
  
  
  //Hole For Top SiPMs in Cu
  
  G4double* z_hole3 = new G4double[2];
  z_hole3[0]= 62.35*cm - 1.9*cm;
  z_hole3[1]= 62.35*cm + 1.8*cm;
  
  
    G4Polycone* hole_3 = new G4Polycone("hole_Cu_top", -8*deg, 16*deg,
  2, z_hole3, rInn_hole2, rOut_hole2); 
  
  
  G4VSolid* Array_hole_1_12 = new G4UnionSolid("Array_hole_1_12", hole_1_1, hole_1_1, ZHoleRot1, G4ThreeVector(0. , 0., 0.));
  G4VSolid* Array_hole_1_13 = new G4UnionSolid("Array_hole_1_13", Array_hole_1_12, hole_1_1, ZHoleRot2, G4ThreeVector(0. , 0., 0.));
  G4VSolid* Array_hole_1_14 = new G4UnionSolid("Array_hole_1_14", Array_hole_1_13, hole_1_1, ZHoleRot3, G4ThreeVector(0. , 0., 0.));  
  G4VSolid* Array_hole_1_15 = new G4UnionSolid("Array_hole_1_15", Array_hole_1_14, hole_1_1, ZHoleRot4, G4ThreeVector(0. , 0., 0.)); 
  G4VSolid* Array_hole_1_16 = new G4UnionSolid("Array_hole_1_16", Array_hole_1_15, hole_1_1, ZHoleRot5, G4ThreeVector(0. , 0., 0.)); 
  G4VSolid* Array_hole_1_17 = new G4UnionSolid("Array_hole_1_17", Array_hole_1_16, hole_1_1, ZHoleRot6, G4ThreeVector(0. , 0., 0.));      
  G4VSolid* Array_hole_1_18 = new G4UnionSolid("Array_hole_1_18", Array_hole_1_17, hole_1_1, ZHoleRot7, G4ThreeVector(0. , 0., 0.)); 
  G4VSolid* Array_hole_2_1 = new G4UnionSolid("Array_hole_2_1", Array_hole_1_18, Array_hole_1_18, 0, G4ThreeVector(0. , 0., -4.42*cm)); 
  G4VSolid* Array_hole_Cu = new G4UnionSolid("Array_hole", Array_hole_2_1, Array_hole_2_1, 0, G4ThreeVector(0. , 0., -2*4.42*cm)); 
  G4VSolid* Array_hole_Cu_51 = new G4UnionSolid("Array_hole_Cu_51",Array_hole_Cu, hole_3, 0, G4ThreeVector(0. , 0., 0.));
  G4VSolid* Array_hole_Cu_52 = new G4UnionSolid("Array_hole_Cu_52",Array_hole_Cu_51, hole_3, ZHoleRot1, G4ThreeVector(0. , 0., 0.));
  G4VSolid* Array_hole_Cu_53 = new G4UnionSolid("Array_hole_Cu_53",Array_hole_Cu_52, hole_3, ZHoleRot2, G4ThreeVector(0. , 0., 0.));
  G4VSolid* Array_hole_Cu_54 = new G4UnionSolid("Array_hole_Cu_54",Array_hole_Cu_53, hole_3, ZHoleRot3, G4ThreeVector(0. , 0., 0.));
  G4VSolid* Array_hole_Cu_55 = new G4UnionSolid("Array_hole_Cu_55",Array_hole_Cu_54, hole_3, ZHoleRot4, G4ThreeVector(0. , 0., 0.));
  G4VSolid* Array_hole_Cu_56 = new G4UnionSolid("Array_hole_Cu_56",Array_hole_Cu_55, hole_3, ZHoleRot5, G4ThreeVector(0. , 0., 0.));
  G4VSolid* Array_hole_Cu_57 = new G4UnionSolid("Array_hole_Cu_57",Array_hole_Cu_56, hole_3, ZHoleRot6, G4ThreeVector(0. , 0., 0.));
  G4VSolid* Array_hole_Cu_final = new G4UnionSolid("Array_hole_Cu_final",Array_hole_Cu_57, hole_3, ZHoleRot7, G4ThreeVector(0. , 0., 0.));



 /////////////////////////////////////////////////////////////////////////
 // Reflector inside the CF4 surrounding the outer Jar
 // Part of PTFE is a 0.508 mm thick clindrical shell
 
 G4double r_inn_refl = 12.259*cm;
 
   G4double* z_PTFE = new G4double[3];
  z_PTFE[0]= 37.417*cm;
  z_PTFE[1]= z_PTFE[0] + 22.225*cm;
  z_PTFE[2]= z_PTFE[1] + 9*cm - 2.5*cm;

  G4double* rInn_PTFE = new G4double[3];
  rInn_PTFE[0]=r_inn_refl;                                 
  rInn_PTFE[1]=r_inn_refl;
  rInn_PTFE[2]=r_inn_refl + 3.1*cm;

  G4double* rOut_PTFE = new G4double[3];
  rOut_PTFE[0]=r_inn_refl + 0.508*mm;
  rOut_PTFE[1]=r_inn_refl + 0.508*mm;
  rOut_PTFE[2]=r_inn_refl + 0.508*mm + 3.1*cm;

  G4Polycone* contenedor_PTFE = new G4Polycone("contenedor_PTFE", 0.*deg, 360*deg,
  3, z_PTFE, rInn_PTFE, rOut_PTFE);

  G4VSolid* subtract_PTFE_Reflector= new G4SubtractionSolid("Hole-Reflector-SiPM", contenedor_PTFE, Array_hole_final, 0, G4ThreeVector(0,0.,0.));      // Hole for SiPM
  


  reflector_PTFE_log = new G4LogicalVolume(subtract_PTFE_Reflector, PTFE_mat, "reflector_PTFE_log");
  reflector_PTFE_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), reflector_PTFE_log, "reflector_PTFE_phys",
                                                         hydraulic_fluid_log, false, 0) ;


  G4VisAttributes* HDPE_yellow= new G4VisAttributes(yellow);
  HDPE_yellow->SetVisibility(true);
  reflector_PTFE_log->SetVisAttributes(HDPE_yellow);




  // Reflector inside the CF4 surrounding the outer Jar
  // Copper part is a 0.3175 cm thick clindrical shell
  
  
  G4double* z_Cu = new G4double[5];
  z_Cu[0]= 37.417*cm;
  z_Cu[1]= 37.417*cm + 1.904*cm;
  z_Cu[2]= 37.417*cm + 1.905*cm;
  z_Cu[3]= z_Cu[0] + 22.225*cm;
  z_Cu[4]= z_Cu[3] + 9*cm - 2.5*cm;
  
  
  G4double* rInn_Cu = new G4double[5];
  rInn_Cu[0]=r_inn_refl + 0.508*mm;
  rInn_Cu[1]=r_inn_refl + 0.508*mm;
  rInn_Cu[2]=r_inn_refl + 0.508*mm;
  rInn_Cu[3]=r_inn_refl + 0.508*mm;    
  rInn_Cu[4]=r_inn_refl + 0.508*mm + 3.1*cm;
  
  
  
  G4double* rOut_Cu = new G4double[5];
  rOut_Cu[0]=14.605*cm;
  rOut_Cu[1]=14.605*cm;
  rOut_Cu[2]=r_inn_refl + 0.508*mm + 3.175*mm;
  rOut_Cu[3]=r_inn_refl + 0.508*mm + 3.175*mm;  
  rOut_Cu[4]=r_inn_refl + 0.508*mm + 3.175*mm + 3.1*cm;
  

 G4Polycone* contenedor_Cu = new G4Polycone("contenedor_Cu", 0.*deg, 360.*deg,
  5, z_Cu, rInn_Cu, rOut_Cu);


  G4VSolid* subtract_Cu_Reflector= new G4SubtractionSolid("Hole-Reflector-SiPM", contenedor_Cu, Array_hole_Cu_final, 0, G4ThreeVector(0,0.,0.));      // Hole for SiPM


  reflector_Cu_log = new G4LogicalVolume(subtract_Cu_Reflector, CuShield_mat, "reflector_Cu_log");           //Copper
  reflector_Cu_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), reflector_Cu_log, "reflector_Cu_phys",
                                                         hydraulic_fluid_log, false, 0) ;

  G4VisAttributes* cu_cil= new G4VisAttributes(cyan);
  cu_cil->SetVisibility(true);
  reflector_Cu_log->SetVisAttributes(cu_cil);


 //Reflector                            
 // Top Part  
 // PTFE Part 0.0508 cm

  G4double* z_top = new G4double[2];
  z_top[0]= 67.0*cm;
  z_top[1]= 71*cm;
  
 
  G4double* rInn_top_PTFE = new G4double[2];
  rInn_top_PTFE[0]=r3 + 4.32*cm;                                  //r3 + 4.32*cm;
  rInn_top_PTFE[1]=0;

  G4double* rOut_top_PTFE = new G4double[2];
  rOut_top_PTFE[0]=r3 + 4.32*cm + 0.0508*cm;                                
  rOut_top_PTFE[1]=0.0508*cm;                                                  

G4Polycone* reflector_top_PTFE = new G4Polycone("reflector_top_PTFE", 0.*deg, 360.*deg,
  2, z_top, rInn_top_PTFE, rOut_top_PTFE);

  reflector_top_PTFE_log = new G4LogicalVolume(reflector_top_PTFE, PTFE_mat, "reflector_top_PTFE_log");           //Copper
  reflector_top_PTFE_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), reflector_top_PTFE_log, "reflector_top_PTFE_phys",
                                                         hydraulic_fluid_log, false, 0) ;


  reflector_top_PTFE_log->SetVisAttributes(HDPE_yellow);




 //Reflector                            
 // Top Part  
 // Cu Part 0.3175 cm

  G4double* rOut_top_Cu = new G4double[2];
  rOut_top_Cu[0]=r3 + 4.32*cm + 0.0508*cm + 0.3175*cm;
  rOut_top_Cu[1]= 0.0508*cm + 0.3175;   

G4Polycone* reflector_top_Cu = new G4Polycone("reflector_top_Cu", 0.*deg, 360.*deg,
  2, z_top, rOut_top_PTFE, rOut_top_Cu);

  reflector_top_Cu_log = new G4LogicalVolume(reflector_top_Cu, CuShield_mat, "reflector_top_Cul_log");           //Copper
  reflector_top_Cu_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), reflector_top_Cu_log, "reflector_top_Cu_phys",
                                                         hydraulic_fluid_log, false, 0) ;

  reflector_top_Cu_log->SetVisAttributes(cu_cil);

/*

  //Led Ring                          

  G4VSolid* LED_Cilindro =new G4Tubs("LED_Cilindro", 3*cm, 3.5*cm, 0.5*mm, 0.*deg, 360.*deg);
  LED_log  = new G4LogicalVolume(LED_Cilindro, LED_mat, "LED_log");
  LED_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,70*cm), LED_log, "LED_phys", hydraulic_fluid_log, false,0);

*/




  //Liquid argon. Sensitive Detector                  

  G4double cil4_hz = cil3_hz - 1*mm ;
  G4double r4 = (12-0.5)*cm ;                  

  G4double* z4 = new G4double[15];
  z4[0]= z3[0];
  z4[1]= z5[3];
  z4[2]= z5[3]+2*r5*cos(25 * PI / 180.0)-sqrt(3)*r5;
  z4[3]= z5[3]+2*r5*cos(20 * PI / 180.0)-sqrt(3)*r5;
  z4[4]= z5[3]+2*r5*cos(15 * PI / 180.0)-sqrt(3)*r5;
  z4[5]= z5[3]+2*r5*cos(10 * PI / 180.0)-sqrt(3)*r5;
  z4[6]= z5[3]+2*r5*cos(5 * PI / 180.0)-sqrt(3)*r5;
  z4[7]= z5[3]+2*r5-sqrt(3)*r5;
  z4[8]= cil4_hz; 
  z4[9]= z4[8]+2*r4*cos(25 * PI / 180.0)-sqrt(3)*r4;
  z4[10]= z4[8]+2*r4*cos(20 * PI / 180.0)-sqrt(3)*r4;
  z4[11]= z4[8]+2*r4*cos(15 * PI / 180.0)-sqrt(3)*r4;
  z4[12]= z4[8]+2*r4*cos(10 * PI / 180.0)-sqrt(3)*r4;
  z4[13]= z4[8]+2*r4*cos(5 * PI / 180.0)-sqrt(3)*r4;
  z4[14]= z4[8]+2*r4-sqrt(3)*r4;

  //G4cout<<"Tamano contenedor grande\n"<<z4[7]-z4[0]<<G4endl;

  G4double* rOut4 = new G4double[15];
  rOut4[0]=r4;
  rOut4[1]=r4;
  rOut4[2]=r4;
  rOut4[3]=r4;
  rOut4[4]=r4;
  rOut4[5]=r4;
  rOut4[6]=r4;  
  rOut4[7]=r4; 
  rOut4[8]=r4; 
  rOut4[9]=2*r4*sin(25 * PI / 180.0);
  rOut4[10]=2*r4*sin(20 * PI / 180.0);
  rOut4[11]=2*r4*sin(15 * PI / 180.0);
  rOut4[12]=2*r4*sin(10 * PI / 180.0);
  rOut4[13]=2*r4*sin(5 * PI / 180.0);
  rOut4[14]=0.; 
  
  
  G4double* rInn4 = new G4double[15];
  rInn4[0]=r5;
  rInn4[1]=r5;
  rInn4[2]=2*r5*sin(25 * PI / 180.0);
  rInn4[3]=2*r5*sin(20 * PI / 180.0);
  rInn4[4]=2*r5*sin(15 * PI / 180.0);
  rInn4[5]=2*r5*sin(10 * PI / 180.0);
  rInn4[6]=2*r5*sin(5 * PI / 180.0);
  rInn4[7]=0.;
  rInn4[8]=0.;  
  rInn4[9]=0;
  rInn4[10]=0;
  rInn4[11]=0;
  rInn4[12]=0;
  rInn4[13]=0;
  rInn4[14]=0.;
  

  G4Polycone* contenedor4 = new G4Polycone("contenedor4", 0.*deg, 360.*deg,
  15, z4, rInn4, rOut4);
 
  
  

  LAr_log = new G4LogicalVolume(contenedor4, LAr_mat, "LAr_log");
  LAr_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), LAr_log, "LAr_phys",
                                                         outer_jar_log, false, 0) ;

  G4VisAttributes* LAr_red= new G4VisAttributes(red);
  LAr_red->SetVisibility(true);
  LAr_log->SetVisAttributes(LAr_red);




   // Top Flange 
   
   G4double top_flange_outr = 34.290*0.5*cm;
   G4double top_flange_innr = 25.324*0.5*cm;
   G4double top_flange_hz = 1.588*cm;
   
   G4double top_sf_outr = 34.29*0.5*cm;
   G4double top_sf_innr = 27.94*0.5*cm;
   G4double top_sf_hz = 1.588*cm;

   
   G4Tubs* Tube_top_flange = new G4Tubs("Tube_top_flange",top_flange_innr, top_flange_outr, top_flange_hz/2, 0.*deg, 360.*deg);  
   top_flange_log = new G4LogicalVolume(Tube_top_flange, vessel_mat, "top_flange_log");
   top_flange_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,25.312*cm - top_flange_hz*0.5), top_flange_log, "top_flange_phys",
                                                         hydraulic_fluid_log, false, 0) ;
                                                                                                              
                                                         

   //Top Support Flange

   
   G4Tubs* Tube_top_sf = new G4Tubs("Tube_top_sf",top_sf_innr, top_sf_outr, top_sf_hz/2, 0.*deg, 360.*deg);
   G4Tubs* Sf_tab = new G4Tubs("Tube_SS",1.1*0.5*cm, 2.5*0.5*cm, top_sf_hz*0.5, 0.*deg, 360.*deg);
   
   G4VSolid* tab_sf_ring1= new G4UnionSolid("tab_sf_ring1", Tube_top_sf, Sf_tab, 0, G4ThreeVector( top_sf_outr + 2.5*0.5*cm - 0.7*cm ,0. , 0.)); 
   G4VSolid* tab_sf_ring2= new G4UnionSolid("tab_sf_ring2", tab_sf_ring1, Sf_tab, 0, G4ThreeVector( cos(120 * PI / 180.0)*(top_sf_outr + 2.5*0.5*cm - 0.7*cm) , sin(120 * PI / 180.0)*(top_sf_outr + 2.5*0.5*cm - 0.7*cm) , 0.)); 
   G4VSolid* tab_sf_ring3= new G4UnionSolid("tab_sf_ring3", tab_sf_ring2, Sf_tab, 0, G4ThreeVector( cos(240 * PI / 180.0)*(top_sf_outr + 2.5*0.5*cm - 0.7*cm) , sin(240 * PI / 180.0)*(top_sf_outr + 2.5*0.5*cm - 0.7*cm) , 0.));
   
   Top_sf_log = new G4LogicalVolume(tab_sf_ring3, vessel_mat, "Top_sf_log");
   Top_sf_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,25.312*cm - top_flange_hz*0.5 - 5.936*cm - top_sf_hz ), Top_sf_log, "Top_sf_phys", hydraulic_fluid_log, false, 0) ;
   
 
 
   // OJ Spacer
  

   G4Tubs* Spacer_Tub = new G4Tubs("Spacer_Tub",0, 2.223*0.5*cm, 5.936*0.5*cm, 0.*deg, 360.*deg);  
   
   OJ_Spacer_log = new G4LogicalVolume(Spacer_Tub, vessel_mat, "Spacer_Tub_log");
   OJ_Spacer_phys = new G4PVPlacement(0, G4ThreeVector(top_sf_outr - 2.223*0.5*cm - 0.2*cm, 0. ,25.312*cm - top_flange_hz - 5.936*0.5*cm), OJ_Spacer_log, "OJ_Spacer1_phys",  hydraulic_fluid_log, false, 0) ;
   
   //Replicate Spacers
   int i;
   char Spacer_Name[30];
   char s_number[30];
   char phys[10]="_phys";
   
   for(i=1; i<=17; i++)   
  {
    strcpy(Spacer_Name,"OJ_Spacer");
    sprintf(s_number, "%d", i+1);
    strcat(Spacer_Name, s_number);
    strcat(Spacer_Name, phys);
    
    OJ_Spacer_phys = new G4PVPlacement(0, G4ThreeVector(cos(i*20 * PI / 180.0)*(top_sf_outr - 2.223*0.5*cm - 0.2*cm), sin(i*20 * PI / 180.0)*(top_sf_outr - 2.223*0.5*cm - 0.2*cm) ,25.312*cm - top_flange_hz - 5.936*0.5*cm), OJ_Spacer_log, Spacer_Name,  hydraulic_fluid_log, false, 0) ;
       
   Spacer_Name[0]='\0' ;
   s_number[0]='\0' ;
    
  }
 

  //Guide Rod Flange
  
  
   Guide_rod_flange_log = new G4LogicalVolume(tab_sf_ring3, vessel_mat, "Guide_rod_flange_log");
   Guide_rod_flange_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,25.312*cm - top_flange_hz*0.5 - 5.936*cm - 2*top_sf_hz - 10.16*cm), Guide_rod_flange_log, "Guide_rod_flange_phys", hydraulic_fluid_log, false, 0) ;
   
   
   
  //Bottom Flange  
  
  
   G4Tubs* Bottom_Flange_Tub = new G4Tubs("Bottom_Flange_Tub",19.685*0.5*cm, 31.433*0.5*cm, 1.905*0.5*cm, 0.*deg, 360.*deg);
   Bottom_Flange_log = new G4LogicalVolume(Bottom_Flange_Tub, vessel_mat, "Bottom_Flange_log");
   Bottom_Flange_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,25.312*cm - top_flange_hz*0.5 - 5.936*cm - 2.5*top_sf_hz - 10.16*cm - 4.3*cm - 1.905*0.5*cm), Bottom_Flange_log, "Bottom_Flange_phys", hydraulic_fluid_log, false, 0) ;
   
    
   G4Tubs* Spacer1_Tub = new G4Tubs("Spacer_Tub",0, 1.905*0.5*cm, 4.3*0.5*cm, 0.*deg, 360.*deg);  
   
   Bottom_Spacer_log = new G4LogicalVolume(Spacer1_Tub, vessel_mat, "Bottom_Spacer_log");
   Bottom_Spacer_phys = new G4PVPlacement(0, G4ThreeVector(31.433*0.5*cm - 1.905*0.5*cm , 0. ,25.312*cm - top_flange_hz*0.5 - 5.936*cm - 2.5*top_sf_hz - 10.16*cm - 4.3*0.5*cm ), Bottom_Spacer_log, "Bottom_Spacer1_phys",  hydraulic_fluid_log, false, 0) ;   
  
  //Replicate Bottom Spacers
  
  
   for(i=1; i<=17; i++)   
  {
    strcpy(Spacer_Name,"Bottom_Spacer");
    sprintf(s_number, "%d", i+1);
    strcat(Spacer_Name, s_number);
    strcat(Spacer_Name, phys);
    
   Bottom_Spacer_phys = new G4PVPlacement(0, G4ThreeVector(cos(i*20 * PI / 180.0)*(31.433*0.5*cm - 1.905*0.5*cm), sin(i*20 * PI / 180.0)*(31.433*0.5*cm - 1.905*0.5*cm) ,25.312*cm - top_flange_hz*0.5 - 5.936*cm - 2.5*top_sf_hz - 10.16*cm - 4.3*0.5*cm), Bottom_Spacer_log, Spacer_Name,  hydraulic_fluid_log, false, 0) ;
       
   Spacer_Name[0]='\0' ;
   s_number[0]='\0' ;
    
  }
  
  
  
  // Base Support Flange
  
  G4Tubs* Base_SF_Tub = new G4Tubs("Base_SF_Tub",31.75*0.5*cm, 38.0*0.5*cm, 1.905*0.5*cm, 0.*deg, 360.*deg);
   
   Base_SF_log = new G4LogicalVolume( Base_SF_Tub, vessel_mat, "Base_SF_log");
   Base_SF_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,25.312*cm - top_flange_hz*0.5 - 5.936*cm - 2.5*top_sf_hz - 10.16*cm - 4.3*cm - 1.905*1.5*cm - 9.670*cm), Base_SF_log, "Base_SF_phys", hydraulic_fluid_log, false, 0) ;
    
  
  //Guide Rods
  
  G4Tubs* Guide_Rod_Tub = new G4Tubs("Guide_Rod_Tub",0., 1.1*0.5*cm, 27.623*0.5*cm, 0.*deg, 360.*deg);
   
   Guide_Rod_log = new G4LogicalVolume(Guide_Rod_Tub, vessel_mat, "Guide_Rod_log");
   Guide_Rod_phys = new G4PVPlacement(0, G4ThreeVector(top_sf_outr + 2.5*0.5*cm - 0.7*cm,0.,25.312*cm - top_flange_hz*0.5 - 5.936*cm - 1.5*top_sf_hz -  27.623*0.5*cm), Guide_Rod_log, "Guide_Rod1_phys", hydraulic_fluid_log, false, 0) ;
   Guide_Rod_phys = new G4PVPlacement(0, G4ThreeVector(cos(120 * PI / 180.0)*(top_sf_outr + 2.5*0.5*cm - 0.7*cm), sin(120 * PI / 180.0)*(top_sf_outr + 2.5*0.5*cm - 0.7*cm) ,25.312*cm - top_flange_hz*0.5 - 5.936*cm - 1.5*top_sf_hz -  27.623*0.5*cm), Guide_Rod_log, "Guide_Rod2_phys", hydraulic_fluid_log, false, 0) ;
   Guide_Rod_phys = new G4PVPlacement(0, G4ThreeVector(cos(240 * PI / 180.0)*(top_sf_outr + 2.5*0.5*cm - 0.7*cm), sin(240 * PI / 180.0)*(top_sf_outr + 2.5*0.5*cm - 0.7*cm) ,25.312*cm - top_flange_hz*0.5 - 5.936*cm - 1.5*top_sf_hz -  27.623*0.5*cm), Guide_Rod_log, "Guide_Rod3_phys", hydraulic_fluid_log, false, 0) ; 
   
   
   
   //Side Support
   
   
   G4Box* Side_support_box = new G4Box
     ("Side_support_box", 0.5*0.635*cm, 0.5*5.08*cm, 0.5*29.120*cm );
     
     
     G4RotationMatrix* zSide = new G4RotationMatrix;         
    zSide->rotateZ(-60*deg);
     G4RotationMatrix* zSide2 = new G4RotationMatrix;         
    zSide2->rotateZ(-180*deg);
     G4RotationMatrix* zSide3 = new G4RotationMatrix;         
    zSide3->rotateZ(-300*deg);
     
    
   Side_support_log = new G4LogicalVolume( Side_support_box, vessel_mat, "Side_support_log");
   Side_support_phys = new G4PVPlacement(zSide, G4ThreeVector(cos(60 * PI / 180.0)*(34.29*0.5*cm + 0.5*0.635*cm) , sin(60 * PI / 180.0)*(34.29*0.5*cm + 0.5*0.635*cm) ,25.312*cm - top_flange_hz*0.5 - 5.936*cm - 2.5*top_sf_hz - 10.16*cm - 4.3*cm - 1.905*cm - 9.670*cm + 0.5*29.120*cm), Side_support_log, "Side_support1_phys", hydraulic_fluid_log, false, 0) ;  
  
   Side_support_phys = new G4PVPlacement(zSide2, G4ThreeVector(cos(180 * PI / 180.0)*(34.29*0.5*cm + 0.5*0.635*cm) , sin(180 * PI / 180.0)*(34.29*0.5*cm + 0.5*0.635*cm) ,25.312*cm - top_flange_hz*0.5 - 5.936*cm - 2.5*top_sf_hz - 10.16*cm - 4.3*cm - 1.905*cm - 9.670*cm + 0.5*29.120*cm), Side_support_log, "Side_support2_phys", hydraulic_fluid_log, false, 0) ;  

   Side_support_phys = new G4PVPlacement(zSide3, G4ThreeVector(cos(300 * PI / 180.0)*(34.29*0.5*cm + 0.5*0.635*cm) , sin(300 * PI / 180.0)*(34.29*0.5*cm + 0.5*0.635*cm) ,25.312*cm - top_flange_hz*0.5 - 5.936*cm - 2.5*top_sf_hz - 10.16*cm - 4.3*cm - 1.905*cm - 9.670*cm + 0.5*29.120*cm), Side_support_log, "Side_support3_phys", hydraulic_fluid_log, false, 0) ;  
     
  
   
   
  
  
  
  //HYSPAN Bellows Assembly
  
  
  G4double* z_HB = new G4double[8];
  z_HB[0]= z5[2];
  z_HB[1]= z_HB[0] + 2.539*cm;
  z_HB[2]= z_HB[0] + 2.54*cm;
  z_HB[3]= z_HB[2] + 12.699*cm;
  z_HB[4]= z_HB[2] + 12.700*cm;
  z_HB[5]= z_HB[2] + 12.700*cm + 0.5*cm;
  z_HB[6]= z_HB[2] + 12.700*cm + 0.51*cm; 
  z_HB[7]= z_HB[6] + 2.525*cm - 0.51*cm;
  
  
  G4double* rOut_HB = new G4double[8];
  rOut_HB[0]=27.146*0.5*cm;
  rOut_HB[1]=27.146*0.5*cm;
  rOut_HB[2]=25.21*0.5*cm;
  rOut_HB[3]=25.21*0.5*cm; 
  rOut_HB[4]=27.94*0.5*cm;
  rOut_HB[5]=27.94*0.5*cm;
  rOut_HB[6]=29.21*0.5*cm;
  rOut_HB[7]=29.21*0.5*cm;

  
  
  G4double* rInn_HB = new G4double[8];
  rInn_HB[0]=21.971*0.5*cm;
  rInn_HB[1]=21.971*0.5*cm;
  rInn_HB[2]=21.971*0.5*cm;
  rInn_HB[3]=21.971*0.5*cm;
  rInn_HB[4]=21.971*0.5*cm;
  rInn_HB[5]=21.971*0.5*cm;
  rInn_HB[6]=21.971*0.5*cm;
  rInn_HB[7]=21.971*0.5*cm;
  
  
  G4Polycone* HB = new G4Polycone("HB", 0.*deg, 360.*deg,
  8, z_HB, rInn_HB, rOut_HB);
  Hyspan_bellow_log = new G4LogicalVolume( HB, vessel_mat, "Hyspan_bellow_log");
  Hyspan_bellow_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), Hyspan_bellow_log, "Hyspan_bellow_phys", hydraulic_fluid_log, false, 0) ;
  
  


/*




  // Insulating  HDPE  //          
  // Inside Inner Jar

  G4double* z11 = new G4double[2];
  z11[0]= 27.5*cm;
  z11[1]= 37.5*cm;

  G4double* rOut11 = new G4double[2];
  rOut11[0]=r6-0.8*cm;
  rOut11[1]=r6-0.8*cm;


  G4Polycone* contenedor11 = new G4Polycone("contenedor11", 0.*deg, 360.*deg,
  2, z11, rInn, rOut11);

  HDPE_inner_jar_log = new G4LogicalVolume(contenedor11, HDPE_mat, "HDPE_inner_jar_log");
  HDPE_inner_jar_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), HDPE_inner_jar_log, "HDPE_inner_jar_phys",
                                                         inner_jar_tick_log, false, 0) ;

HDPE_inner_jar_log->SetVisAttributes(HDPE_yellow);

*/




  //  Piezo  //



  // Cu Container                 

    G4double Cu_radio = 1.53*cm ;
    G4double Cu_hz = 2.54*cm ;

    G4RotationMatrix* ypiezoRot = new G4RotationMatrix;               // Rotar el cilindro de Cobre
    ypiezoRot->rotateY(90*deg);


  G4Tubs* Cu_Cil = new G4Tubs("Cu_Cil",0, Cu_radio, Cu_hz/2, 0.*deg, 360.*deg);
  Cu_log = new G4LogicalVolume(Cu_Cil, CuShield_mat, "Cu_log");
  Piezo_Cu_phys = new G4PVPlacement(ypiezoRot, G4ThreeVector(13.3*cm,0.,28.805*cm), Cu_log, "Piezo_Cu_1_phys", hydraulic_fluid_log, false, 0) ;              // Esta dentro del liquido hidraulico
  

  Cu_log->SetVisAttributes(cu_cil);



  // Piezo container        //Está hecho de MAS Epoxy       

  G4double Vacio_radio = Cu_radio - 1*mm ;
  G4double Vacio_hz = Cu_hz - 2*mm ;

  G4Tubs* Vacio_Cil = new G4Tubs("Vacio_Cil",0, Vacio_radio, Vacio_hz/2, 0.*deg, 360.*deg);
  Vacio_log = new G4LogicalVolume(Vacio_Cil, Epoxy_mat, "Vacio_log");
  Vacio_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), Vacio_log, "Vacio_phys",
                                                         Cu_log, false, 0) ;




  // PCB Container     (Printed Circuit Board)

  G4double PBC_radio = Vacio_radio ;
  G4double PBC_hz = 1.7*mm ;

  G4Tubs* PBC_Cil = new G4Tubs("PBC_Cil",0, PBC_radio, PBC_hz/2, 0.*deg, 360.*deg);
  PBC_log = new G4LogicalVolume(PBC_Cil, PBC_mat, "PBC_log");
  PBC_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,-1.085*cm), PBC_log, "PBC_phys",
                                                         Vacio_log, false, 0) ;

  G4VisAttributes* pbc_cil= new G4VisAttributes(red);
  pbc_cil->SetVisibility(true);
  PBC_log->SetVisAttributes(pbc_cil);



  // Piezo Cylinder

  G4double Piezo_radio = (4.75/2)*mm ;
  G4double Piezo_hz = 2.17*cm - PBC_hz;

  G4Tubs* Piezo_Cil = new G4Tubs("Piezo_Cil",0, Piezo_radio, Piezo_hz/2, 0.*deg, 360.*deg);
  Piezo_log = new G4LogicalVolume(Piezo_Cil, Piezo_mat, "Piezo_log");                             // Poner material luego
  Piezo_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), Piezo_log, "Piezo_phys", Vacio_log, false, 0) ;


  G4VisAttributes* piezo_cil= new G4VisAttributes(blue);
  piezo_cil->SetVisibility(true);
  Piezo_log->SetVisAttributes(piezo_cil);


  //Silver Epoxy

G4Tubs* S_Epoxy_Cil = new G4Tubs("Piezo_Cil",0, Piezo_radio, PBC_hz/2, 0.*deg, 360.*deg);

  S_Epoxy_log = new G4LogicalVolume(S_Epoxy_Cil, S_Epoxy_mat, "S_Epoxy_log");                             // Poner material luego
  S_Epoxy_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,1.085*cm), S_Epoxy_log, "S_Epoxy_phys", Vacio_log, false, 0) ;

S_Epoxy_log->SetVisAttributes(quartz_green);



 //Insulating  HDPE Castle inside pressure vessel                   

  G4double* z12 = new G4double[9];
  z12[0]= 25.315*cm;
  z12[1]= z12[0] + 6.985*cm;
  z12[2]= z12[1] + 5.079*cm;
  z12[3]= z12[2] + 0.001*cm;
  z12[4]= z12[3] + 18.414*cm;
  z12[5]= z12[4] +  0.001*cm;  
  z12[6]= z12[5] +  2.539*cm;    
  z12[7]= z12[6] +  0.001*cm; 
  z12[8]= z12[7] +  6.985*cm;   
  
 G4cout<<"Altura HDPE Castle\n"<<z12[8]-z12[0]<<G4endl;    
  
  
  G4double* rInn12 = new G4double[9];
  rInn12[0]=12.637*cm; 
  rInn12[1]=12.637*cm;
  rInn12[2]=12.637*cm;
  rInn12[3]=14.923*cm;
  rInn12[4]=14.923*cm;
  rInn12[5]=14.923*cm; 
  rInn12[6]=14.923*cm;  
  rInn12[7]=15.558*cm; 
  rInn12[8]=15.558*cm;   
  
  
  

  G4double* rOut12 = new G4double[9];
  rOut12[0]=18.733*cm;
  rOut12[1]=18.733*cm;
  rOut12[2]=18.733*cm;
  rOut12[3]=18.733*cm;
  rOut12[4]=18.733*cm;  
  rOut12[5]=18.098*cm;   
  rOut12[6]=18.098*cm;   
  rOut12[7]=18.098*cm;   
  rOut12[8]=18.098*cm; 
  
  
  G4RotationMatrix * PiezoRot = new G4RotationMatrix();
   PiezoRot -> rotateY(90.0*deg);
  
  
  G4Polycone* contenedor12 = new G4Polycone("contenedor12", 0.*deg, 360.*deg,
  9, z12, rInn12, rOut12);
  
    G4VSolid* subtract_HDPE_1 = new G4SubtractionSolid("subtract_HDPE_1", contenedor12, Cu_Cil, PiezoRot, G4ThreeVector(13.3*cm,0.,28.805*cm));   
    PiezoRot -> rotateX(-45*deg);  
    G4VSolid* subtract_HDPE_2 = new G4SubtractionSolid("subtract_HDPE_2", subtract_HDPE_1, Cu_Cil, PiezoRot, G4ThreeVector(cos(45 * PI / 180.0)*13.3*cm, sin(45 * PI / 180.0)*13.3*cm,28.805*cm));
    PiezoRot -> rotateX(-45*deg);  
    G4VSolid* subtract_HDPE_3 = new G4SubtractionSolid("subtract_HDPE_3", subtract_HDPE_2, Cu_Cil, PiezoRot, G4ThreeVector(cos(90 * PI / 180.0)*13.3*cm, sin(90 * PI / 180.0)*13.3*cm,28.805*cm));
    PiezoRot -> rotateX(-45*deg);  
    G4VSolid* subtract_HDPE_4 = new G4SubtractionSolid("subtract_HDPE_4", subtract_HDPE_3, Cu_Cil, PiezoRot, G4ThreeVector(cos(135 * PI / 180.0)*13.3*cm, sin(135 * PI / 180.0)*13.3*cm,28.805*cm));
    PiezoRot -> rotateX(-45*deg);  
    G4VSolid* subtract_HDPE_5 = new G4SubtractionSolid("subtract_HDPE_5", subtract_HDPE_4, Cu_Cil, PiezoRot, G4ThreeVector(cos(180 * PI / 180.0)*13.3*cm, sin(180 * PI / 180.0)*13.3*cm,28.805*cm));
    PiezoRot -> rotateX(-45*deg);  
    G4VSolid* subtract_HDPE_6 = new G4SubtractionSolid("subtract_HDPE_6", subtract_HDPE_5, Cu_Cil, PiezoRot, G4ThreeVector(cos(225 * PI / 180.0)*13.3*cm, sin(225 * PI / 180.0)*13.3*cm,28.805*cm));
  PiezoRot -> rotateX(-45*deg);  
    G4VSolid* subtract_HDPE_7 = new G4SubtractionSolid("subtract_HDPE_7", subtract_HDPE_6, Cu_Cil, PiezoRot, G4ThreeVector(cos(270 * PI / 180.0)*13.3*cm, sin(270 * PI / 180.0)*13.3*cm,28.805*cm));
    PiezoRot -> rotateX(-45*deg);
    G4VSolid* subtract_HDPE_8 = new G4SubtractionSolid("subtract_HDPE_8", subtract_HDPE_7, Cu_Cil, PiezoRot, G4ThreeVector(cos(315 * PI / 180.0)*13.3*cm, sin(315 * PI / 180.0)*13.3*cm,28.805*cm));
  

  HDPE_pressure_vessel_log = new G4LogicalVolume(subtract_HDPE_8, HDPE_mat, "HDPE_pressure_vessel_log");
  HDPE_pressure_vessel_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), HDPE_pressure_vessel_log, "HDPE_pressure_vessel_phys", hydraulic_fluid_log, false, 0) ;
   
  G4VisAttributes* hdpe_blue= new G4VisAttributes(blue);
  hdpe_blue->SetVisibility(true);
  HDPE_pressure_vessel_log->SetVisAttributes(hdpe_blue);                                                     







/*

  // Piezo wire                 //Aqui comienzo a construir el cable que va hasta el piezo
  // Plastic cylinder           //Cilindro de plastico

  G4double Cable_Diameter = 1.1*mm ;
  G4double Cable_hz = 37*cm ;
  G4double Cobre_Diameter = 0.86*mm ;
  G4double Cobre_hz = Cable_hz ;

  G4Tubs* Cable_Cil = new G4Tubs("Cable_Cil",0, Cable_Diameter/2, Cable_hz/2, 0.*deg, 360.*deg);
  Cable_log = new G4LogicalVolume(Cable_Cil, Cable_mat, "Cable_log");                            //Poner material correcto
  Cable_phys = new G4PVPlacement(0, G4ThreeVector(12.25*cm,0.,19.5*cm), Cable_log, "Cable_phys",
                                                         contenedor2_log, false, 0) ;


  // Cu wire                   //Cable interior de Cobre

  G4Tubs* Cobre_Cil = new G4Tubs("Cobre_Cil",0, Cobre_Diameter/2, Cobre_hz/2, 0.*deg, 360.*deg);
  Cobre_log = new G4LogicalVolume(Cobre_Cil, CuShield_mat, "Cobre_log");
  Cobre_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), Cobre_log, "Cobre_phys",
                                                         Cable_log, false, 0) ;

*/





  //Calibration Port
 
  G4double Tube_diameter = 5.08*cm;
  G4double Tube_hz = 24*cm;
  G4double Tube_z = z1[7] + Tube_hz/2 + 1.5*cm + 2*cm;
    
  

  G4Tubs* Tube_SS = new G4Tubs("Tube_SS",0, Tube_diameter/2, Tube_hz/2, 0.*deg, 360.*deg);
  calibration_port_log = new G4LogicalVolume(Tube_SS, vessel_mat, "calibration_port_log");
  calibration_port_phys = new G4PVPlacement(0, G4ThreeVector(9.3*cm,0.,Tube_z), calibration_port_log, "calibration_port_phys", Inside_vacuum_vessel_log, false, 0) ;

 calibration_port_log->SetVisAttributes(ssteal_blue);  

// Inside Calibration Port


  G4double Tube_radio = Tube_diameter/2 - 0.1651*cm;
  G4double Tube_Cap_THK = 0.635*cm;
  G4double Tube_hz1 = Tube_hz - Tube_Cap_THK;

  
  
  G4Tubs* Tube_Air = new G4Tubs("Tube_Air",0, Tube_radio, Tube_hz1/2, 0.*deg, 360.*deg);
  calibration_air_log = new G4LogicalVolume(Tube_Air, lab_mat, "calibration_air_log");
  calibration_air_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,Tube_Cap_THK/2), calibration_air_log, "calibration_air_phys",
                                                         calibration_port_log, false, 0) ;

  
  // Berylium Source

  G4double Tube_Be_z = 7.62*cm;

  G4Tubs* Tube_Be = new G4Tubs("Tube_Be",0, Tube_radio, Tube_Be_z/2, 0.*deg, 360.*deg);
  calibration_Be_log = new G4LogicalVolume(Tube_Be, Berilio_mat, "calibration_Be_log");
  calibration_Be_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,-Tube_hz1/2+Tube_Be_z/2), calibration_Be_log, "calibration_Be_phys",
                                                         calibration_air_log, false, 0) ;
 
  G4VisAttributes* Si_brown= new G4VisAttributes(brown);
  Si_brown->SetVisibility(true);
  calibration_Be_log->SetVisAttributes(Si_brown);
 
/*
 //Camara Container

  G4VSolid* Container_camara_Cil =new G4Tubs("Container_camara_Cil", 0.*cm, 5*cm, 4*cm, 0.*deg, 360.*deg);
  Container_camara_log  = new G4LogicalVolume(Container_camara_Cil, vacuum_mat, "Container_camara_log");                         
  Container_camara_phys = new G4PVPlacement(yRot, G4ThreeVector(-25*cm,0.,111.0094*cm), Container_camara_log, "Container_camara_phys", world_log, false,0);


//Lente Camara



  G4VSolid* Lente_Cilindro =new G4Tubs("Lente_Cilindro", 0.*cm, 1.5*cm, 4.0*0.5*cm, 0.*deg, 360.*deg);
  Lente_log  = new G4LogicalVolume(Lente_Cilindro, panel_mat, "Lente_log");                          // Al
  Lente_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,-2*cm), Lente_log, "Lente_phys", Container_camara_log, false,0);

  G4VisAttributes* Al_orange= new G4VisAttributes(orange);
  Al_orange->SetVisibility(true);
  Lente_log->SetVisAttributes(Al_orange);


// Camara


  G4Box* Camara_box = new G4Box
     ("Camara_box", 0.5*4*cm, 0.5*4*cm, 0.5*4*cm );
  Camara_log  = new G4LogicalVolume(Camara_box, panel_mat, "Camara_log");                // Aluminio
  Camara_phys = new G4PVPlacement(0, G4ThreeVector(0.,0.,2*cm), Camara_log, "Camara_phys", Container_camara_log, false,0);

  Camara_log->SetVisAttributes(Al_orange);


//Light Guide   Corregir posición

G4Box* Guia_Cilindro = new G4Box
     ("Guia_Cilindro", 0.5*1.27*cm, 0.5*1.27*cm, 0.5*38.48*cm );

  Guide_log  = new G4LogicalVolume(Guia_Cilindro, acrylic_mat, "Guide_log"); 
  Guide_phys = new G4PVPlacement(yRot, G4ThreeVector(-25*cm + 23.24*0.3827*cm,0.,111.0094*cm - 23.24*0.9239*cm), Guide_log, "Guide_phys", Inside_vacuum_vessel_log, false,0);   

  Guide_log->SetVisAttributes(quartz_green);
*/

 // First SiPM Holder //

G4RotationMatrix* y2Rot = new G4RotationMatrix;               
y2Rot->rotateY(90.0*deg);

G4RotationMatrix* zRot = new G4RotationMatrix;               
  zRot->rotateZ(90.0*deg);
G4RotationMatrix* z1Rot = new G4RotationMatrix;               
  z1Rot->rotateZ(180.0*deg);
G4RotationMatrix* z2Rot = new G4RotationMatrix;               
  z2Rot->rotateZ(270.0*deg);

G4Box* holder_box = new G4Box
     ("holder_box", 0.5*2.0*cm, 0.5*2.0*cm, 0.5*4.5*mm );


G4Box* holder_arm = new G4Box
     ("holder_arm", 0.5*3.7*mm, 0.5*2.0*mm, 0.5*3.0*mm );
     

G4Tubs* holder_ring = new G4Tubs("holder_ring",1.5*mm, 2.5*mm, 0.5*3*mm, 0.*deg, 360.*deg);

G4VSolid* arm_ring= new G4UnionSolid("arm_ring", holder_arm, holder_ring, 0, G4ThreeVector(4.35*mm - 0.08*mm,0.,0.)); 

G4VSolid* holder_1= new G4UnionSolid("holder_1", holder_box, arm_ring, 0, G4ThreeVector(11.85*mm ,0.,-0.75*mm)); 

G4VSolid* holder_2= new G4UnionSolid("holder_2", holder_1, arm_ring, zRot, G4ThreeVector(0., -11.85*mm, -0.75*mm)); 

G4VSolid* holder_3= new G4UnionSolid("holder_3", holder_2, arm_ring, z1Rot, G4ThreeVector(-11.85*mm, 0., -0.75*mm)); 

G4VSolid* holder_4= new G4UnionSolid("holder_4", holder_3, arm_ring, z2Rot, G4ThreeVector(0., 11.85*mm, -0.75*mm)); 


 SiPM_Holder1_log  = new G4LogicalVolume(holder_4, Photopolymer_mat, "SiPM_Holder1_log");
 SiPM_Holder1_phys = new G4PVPlacement(y2Rot, G4ThreeVector(12.6*cm ,0.,55.28*cm), SiPM_Holder1_log, "SiPM_Holder11_phys", hydraulic_fluid_log, false,0);  
         

SiPM_Holder1_log->SetVisAttributes(HDPE_yellow);





  G4RotationMatrix* y1Rot = new G4RotationMatrix;               
  y1Rot->rotateY(-90.0*deg);
  
  

 G4Box* SiPM_Inside_box = new G4Box                                                                                   //CF4 box containing silicon and quartz
     ("SiPM_Inside_box", 0.5*3.787*mm, 0.5*15*mm, 0.5*15*mm );

 SiPM1_Inside_log  = new G4LogicalVolume(SiPM_Inside_box, the_CF4, "SiPM1_Inside_log"); 
 SiPM1_Inside_phys = new G4PVPlacement(y1Rot, G4ThreeVector(0.,0.,+0.5*0.713*mm), SiPM1_Inside_log, "SiPM1_Inside_phys", SiPM_Holder1_log, false,0);   

 // SiPM1_Inside_log->SetVisAttributes(CF4_cyan);


  G4Box* SIPM_Si_box = new G4Box                                                                                      // Silicon array
     ("SIPM_Si_box", 0.5*3.0*mm, 0.5*10.0*mm, 0.5*10.0*mm );


  SIPM1_Si_log  = new G4LogicalVolume(SIPM_Si_box, SIPM_mat, "SIPM1_Si_log");
  SIPM1_Si_phys = new G4PVPlacement(0, G4ThreeVector(-0.5*0.787*mm,0.,0.), SIPM1_Si_log, "SIPM1_Si_phys", SiPM1_Inside_log, false,0);
  
  
    G4VisAttributes* Si_red= new G4VisAttributes(red);
  Si_red->SetVisibility(true);
  SIPM1_Si_log->SetVisAttributes(Si_red);
  
  
   //First SiPM PCB. Inside
  
  G4Box* SiPM_PCB_C_box = new G4Box                                                                                   //SiPM PCB
     ("SiPM_PCB_box", 0.5*0.787*mm, 0.5*15*mm, 0.5*15*mm );

  SiPM_PCB_Inn1_log  = new G4LogicalVolume(SiPM_PCB_C_box, HDPE_mat, "SiPM_PCB_Inn1_log");
  SiPM_PCB_Inn1_phys = new G4PVPlacement(0, G4ThreeVector(0.5*3.0*mm,0.,0.), SiPM_PCB_Inn1_log, "SiPM_PCB_Inn1_phys", SiPM1_Inside_log, false,0);

  SiPM_PCB_Inn1_log->SetVisAttributes(ssteal_blue);
  
  
  //First SiPM PCB. Outside
  
  G4Box* SiPM_plastic_box = new G4Box                                               //Debe ser eliminado
     ("SiPM__plastic_box", 0.5*0.191*cm, 0.5*1.842*cm, 0.5*1.588*cm );
  
  SiPM_PCB_Out1_log  = new G4LogicalVolume(SiPM_PCB_C_box, HDPE_mat, "SiPM_PCB_Out1_log");
  SiPM_PCB_Out1_phys = new G4PVPlacement(0, G4ThreeVector(12.87*cm ,0.,55.28*cm), SiPM_PCB_Out1_log, "SiPM_PCB_Out11_phys", hydraulic_fluid_log, false,0);
  
   SiPM_PCB_Out1_log->SetVisAttributes(ssteal_blue);
  
  
 


// Second SiPM Holder

  SiPM_Holder2_log  = new G4LogicalVolume(holder_4, Photopolymer_mat, "SiPM_Holder2_log");
  SiPM_Holder2_phys = new G4PVPlacement(y2Rot, G4ThreeVector(12.6*cm ,0., 50.86*cm), SiPM_Holder2_log, "SiPM_Holder21_phys", hydraulic_fluid_log, false,0); 

SiPM_Holder2_log->SetVisAttributes(HDPE_yellow);



// Second SiPM //




  SiPM2_Inside_log  = new G4LogicalVolume(SiPM_Inside_box, the_CF4, "SiPM2_Inside_log"); 
  SiPM2_Inside_phys = new G4PVPlacement(y1Rot, G4ThreeVector(0.,0.,+0.5*0.713*mm), SiPM2_Inside_log, "SiPM2_Inside_phys", SiPM_Holder2_log, false,0);   

  //SiPM2_Inside_log->SetVisAttributes(CF4_cyan);
  
 

  SIPM2_Si_log  = new G4LogicalVolume(SIPM_Si_box, SIPM_mat, "SIPM2_Si_log");
  SIPM2_Si_phys = new G4PVPlacement(0, G4ThreeVector(-0.5*0.787*mm,0.,0.), SIPM2_Si_log, "SIPM2_Si_phys", SiPM2_Inside_log, false,0);

   SIPM2_Si_log->SetVisAttributes(Si_red);
   
   
   
   
   //Second SiPM PCB. Inside
  

  SiPM_PCB_Inn2_log  = new G4LogicalVolume(SiPM_PCB_C_box, HDPE_mat, "SiPM_PCB_Inn2_log");
  SiPM_PCB_Inn2_phys = new G4PVPlacement(0, G4ThreeVector(0.5*3.0*mm,0.,0.), SiPM_PCB_Inn2_log, "SiPM_PCB_Inn2_phys", SiPM2_Inside_log, false,0);

  SiPM_PCB_Inn2_log->SetVisAttributes(ssteal_blue);
   


 //Second SiPM PCB. Outside
  
  
  SiPM_PCB_Out2_log  = new G4LogicalVolume(SiPM_PCB_C_box, HDPE_mat, "SiPM_PCB_Out2_log");
  SiPM_PCB_Out2_phys = new G4PVPlacement(0, G4ThreeVector(12.87*cm ,0.,50.86*cm), SiPM_PCB_Out2_log, "SiPM_PCB_Out21_phys", hydraulic_fluid_log, false,0);
  
   SiPM_PCB_Out2_log->SetVisAttributes(ssteal_blue);

 


 


// Third SiPM Holder

  SiPM_Holder3_log  = new G4LogicalVolume(holder_4, Photopolymer_mat, "SiPM_Holder3_log");
  SiPM_Holder3_phys = new G4PVPlacement(y2Rot, G4ThreeVector(12.6*cm ,0., 46.44*cm), SiPM_Holder3_log, "SiPM_Holder31_phys", hydraulic_fluid_log, false,0); 

SiPM_Holder3_log->SetVisAttributes(HDPE_yellow);


// Third SiPM //




 SiPM3_Inside_log  = new G4LogicalVolume(SiPM_Inside_box, the_CF4, "SiPM3_Inside_log"); 
 SiPM3_Inside_phys = new G4PVPlacement(y1Rot, G4ThreeVector(0.,0.,+0.5*0.713*mm), SiPM3_Inside_log, "SiPM3_Inside_phys", SiPM_Holder3_log, false,0);   

  //SiPM3_Inside_log->SetVisAttributes(CF4_cyan);
  
  

  SIPM3_Si_log  = new G4LogicalVolume(SIPM_Si_box, SIPM_mat, "SIPM3_Si_log");
  SIPM3_Si_phys = new G4PVPlacement(0, G4ThreeVector(-0.5*0.787*mm,0.,0.), SIPM3_Si_log, "SIPM3_Si_phys", SiPM3_Inside_log, false,0);

   SIPM3_Si_log->SetVisAttributes(Si_red);
  
  
  
  //Third SiPM PCB. Inside
  

  SiPM_PCB_Inn3_log  = new G4LogicalVolume(SiPM_PCB_C_box, HDPE_mat, "SiPM_PCB_Inn3_log");
  SiPM_PCB_Inn3_phys = new G4PVPlacement(0, G4ThreeVector(0.5*3.0*mm,0.,0.), SiPM_PCB_Inn3_log, "SiPM_PCB_Inn3_phys", SiPM3_Inside_log, false,0);

  SiPM_PCB_Inn3_log->SetVisAttributes(ssteal_blue);
   


 //Third SiPM PCB. Outside
  
  
  SiPM_PCB_Out3_log  = new G4LogicalVolume(SiPM_PCB_C_box, HDPE_mat, "SiPM_PCB_Out3_log");
  SiPM_PCB_Out3_phys = new G4PVPlacement(0, G4ThreeVector(12.87*cm ,0.,46.44*cm), SiPM_PCB_Out3_log, "SiPM_PCB_Out31_phys", hydraulic_fluid_log, false,0);
  
   SiPM_PCB_Out3_log->SetVisAttributes(ssteal_blue);












// Fourth SiPM Holder

  SiPM_Holder4_log  = new G4LogicalVolume(holder_4, Photopolymer_mat, "SiPM_Holder4_log");
  SiPM_Holder4_phys = new G4PVPlacement(y2Rot, G4ThreeVector(12.6*cm ,0., 42.02*cm), SiPM_Holder4_log, "SiPM_Holder41_phys", hydraulic_fluid_log, false,0); 

SiPM_Holder4_log->SetVisAttributes(HDPE_yellow);


// Fourth SiPM //


 SiPM4_Inside_log  = new G4LogicalVolume(SiPM_Inside_box, the_CF4, "SiPM4_Inside_log"); 
 SiPM4_Inside_phys = new G4PVPlacement(y1Rot, G4ThreeVector(0.,0.,+0.5*0.713*mm), SiPM4_Inside_log, "SiPM4_Inside_phys", SiPM_Holder4_log, false,0);   

 // SiPM4_Inside_log->SetVisAttributes(CF4_cyan);
  
  
  
 

  SIPM4_Si_log  = new G4LogicalVolume(SIPM_Si_box, SIPM_mat, "SIPM4_Si_log");
  SIPM4_Si_phys = new G4PVPlacement(0, G4ThreeVector(-0.5*0.787*mm,0.,0.), SIPM4_Si_log, "SIPM4_Si_phys", SiPM4_Inside_log, false,0);

  SIPM4_Si_log->SetVisAttributes(Si_red);



 //Fourth SiPM PCB. Inside
  

  SiPM_PCB_Inn4_log  = new G4LogicalVolume(SiPM_PCB_C_box, HDPE_mat, "SiPM_PCB_Inn4_log");
  SiPM_PCB_Inn4_phys = new G4PVPlacement(0, G4ThreeVector(0.5*3.0*mm,0.,0.), SiPM_PCB_Inn4_log, "SiPM_PCB_Inn4_phys", SiPM4_Inside_log, false,0);

  SiPM_PCB_Inn4_log->SetVisAttributes(ssteal_blue);
   


 //Fourth SiPM PCB. Outside
  
  
  SiPM_PCB_Out4_log  = new G4LogicalVolume(SiPM_PCB_C_box, HDPE_mat, "SiPM_PCB_Out4_log");
  SiPM_PCB_Out4_phys = new G4PVPlacement(0, G4ThreeVector(12.87*cm ,0.,42.02*cm), SiPM_PCB_Out4_log, "SiPM_PCB_Out41_phys", hydraulic_fluid_log, false,0);
  
   SiPM_PCB_Out4_log->SetVisAttributes(ssteal_blue);






 
  
  
  // Fifth  SiPM Holder (Upper Row) 
  
  G4RotationMatrix* SiPM_5_Rot = new G4RotationMatrix;               
  SiPM_5_Rot->rotateY(64.5*deg);

  SiPM_Holder5_log  = new G4LogicalVolume(holder_4, Photopolymer_mat, "SiPM_Holder5_log");
  SiPM_Holder5_phys = new G4PVPlacement(SiPM_5_Rot, G4ThreeVector(13.99*cm ,0., 62.35*cm), SiPM_Holder5_log, "SiPM_Holder51_phys", hydraulic_fluid_log, false,0); 

 SiPM_Holder5_log->SetVisAttributes(HDPE_yellow);
  
  
  
  
  //Fifth SiPM (Upper Row)  
  


 SiPM5_Inside_log  = new G4LogicalVolume(SiPM_Inside_box, the_CF4, "SiPM5_Inside_log"); 
 SiPM5_Inside_phys = new G4PVPlacement(y1Rot, G4ThreeVector(0.,0.,+0.5*0.713*mm), SiPM5_Inside_log, "SiPM5_Inside_phys",  SiPM_Holder5_log, false,0);   

  //SiPM5_Inside_log->SetVisAttributes(CF4_cyan);
  
 

  SIPM5_Si_log  = new G4LogicalVolume(SIPM_Si_box, SIPM_mat, "SIPM5_Si_log");
  SIPM5_Si_phys = new G4PVPlacement(0, G4ThreeVector(-0.5*0.787*mm,0.,0.), SIPM5_Si_log, "SIPM5_Si_phys", SiPM5_Inside_log, false,0);

  SIPM5_Si_log->SetVisAttributes(Si_red);




  //Fifth SiPM PCB. Inside
  

  SiPM_PCB_Inn5_log  = new G4LogicalVolume(SiPM_PCB_C_box, HDPE_mat, "SiPM_PCB_Inn5_log");
  SiPM_PCB_Inn5_phys = new G4PVPlacement(0, G4ThreeVector(0.5*3.0*mm,0.,0.), SiPM_PCB_Inn5_log, "SiPM_PCB_Inn5_phys", SiPM5_Inside_log, false,0);

  SiPM_PCB_Inn5_log->SetVisAttributes(ssteal_blue);




  //Fifth SiPM PCB. Outside
  
  G4RotationMatrix* SiPM_PCB_Rot = new G4RotationMatrix;               
   SiPM_PCB_Rot->rotateY(-25.5*deg);
  
  
  SiPM_PCB_Out5_log  = new G4LogicalVolume(SiPM_PCB_C_box, HDPE_mat, "SiPM_PCB_Out5_log");
  SiPM_PCB_Out5_phys = new G4PVPlacement(SiPM_PCB_Rot, G4ThreeVector(14.30*cm ,0.,62.35*cm), SiPM_PCB_Out5_log, "SiPM_PCB_Out51_phys", hydraulic_fluid_log, false,0);
  
   SiPM_PCB_Out5_log->SetVisAttributes(ssteal_blue);




















  
  
  
  
  
  
  //Replicate SiPMs and Piezos
  ///////////////////////////   
  
  
  G4double SiPM_r = 12.6*cm ;
  G4double SiPM5_r = 13.99*cm;
  G4double PCB_r = 12.87*cm ;
  G4double PCB5_r = 14.30*cm ;
  G4double Piezo_r = 13.3*cm ;
  
  G4double SiPM_Angle = 2.0*3.1415/8.0;
  int j;

  char SiPM_Name[30];
  char Piezo_Name[30];
  char SiPM_PCB_Name[30];
  char number[30];
  char numberj[30];
 
  for(j=1; j<=5; j++)   
  {
  
  for(i=2;i<=8;i++)
  {
   strcpy(SiPM_Name,"SiPM_Holder");
   strcpy(SiPM_PCB_Name,"SiPM_PCB_Out");
   strcpy(Piezo_Name,"Piezo_Cu_");
   sprintf(number, "%d", i);
   sprintf(numberj, "%d", j);
   strcat(SiPM_Name, numberj);
   strcat(SiPM_Name, number);
   strcat(SiPM_Name, phys);
   strcat(SiPM_PCB_Name, numberj);
   strcat(SiPM_PCB_Name, number);
   strcat(SiPM_PCB_Name, phys);
   strcat(Piezo_Name, number);
   strcat(Piezo_Name, phys);
   
   G4RotationMatrix * SiPM1Rot = new G4RotationMatrix();
   SiPM1Rot -> rotateY(90.0*deg);
   SiPM1Rot -> rotateX(-(i-1)*SiPM_Angle*rad); 
   
   
     G4ThreeVector AxisYRot = G4ThreeVector(0,0,1).unit();
   
    G4RotationMatrix * SiPM5Rot = new G4RotationMatrix();
    SiPM5Rot -> rotateZ(-(i-1)*SiPM_Angle*rad); 
    SiPM5Rot -> rotateY(64.5*deg);
   
   
   G4RotationMatrix * PlasticRot = new G4RotationMatrix();
   PlasticRot -> rotateZ(-(i-1)*SiPM_Angle*rad); 
   
   G4RotationMatrix * Plastic5Rot = new G4RotationMatrix();
   Plastic5Rot -> rotateZ(-(i-1)*SiPM_Angle*rad); 
   Plastic5Rot -> rotateY(-25.5*deg);

   
   
   if(j==1)     //First Row
       {
   SiPM_Holder1_phys = new G4PVPlacement(SiPM1Rot, G4ThreeVector(SiPM_r*cos((i-1)*SiPM_Angle) ,SiPM_r*sin((i-1)*SiPM_Angle),55.28*cm), SiPM_Holder1_log, SiPM_Name, hydraulic_fluid_log, false,0);
   SiPM_PCB_Out1_phys = new G4PVPlacement(PlasticRot, G4ThreeVector(PCB_r*cos((i-1)*SiPM_Angle) ,PCB_r*sin((i-1)*SiPM_Angle), 55.28*cm), SiPM_PCB_Out1_log, SiPM_PCB_Name, hydraulic_fluid_log, false,0);
   
   
   
   //Piezo, only one row
   
  Piezo_Cu_phys = new G4PVPlacement(SiPM1Rot, G4ThreeVector(Piezo_r*cos((i-1)*SiPM_Angle) ,Piezo_r*sin((i-1)*SiPM_Angle),28.805*cm), Cu_log, Piezo_Name, hydraulic_fluid_log, false, 0) ; 
   
      }
      
      
      
   if(j==2)   //Second Row
      {   
      SiPM_Holder2_phys = new G4PVPlacement(SiPM1Rot, G4ThreeVector(SiPM_r*cos((i-1)*SiPM_Angle) ,SiPM_r*sin((i-1)*SiPM_Angle),50.86*cm), SiPM_Holder2_log, SiPM_Name, hydraulic_fluid_log, false,0);
     SiPM_PCB_Out2_phys = new G4PVPlacement(PlasticRot, G4ThreeVector(PCB_r*cos((i-1)*SiPM_Angle) ,PCB_r*sin((i-1)*SiPM_Angle), 50.86*cm), SiPM_PCB_Out2_log, SiPM_PCB_Name, hydraulic_fluid_log, false,0);  
      }
      
      
      
      if(j==3)   //Third Row
      {      
      SiPM_Holder3_phys = new G4PVPlacement(SiPM1Rot, G4ThreeVector(SiPM_r*cos((i-1)*SiPM_Angle) ,SiPM_r*sin((i-1)*SiPM_Angle),46.44*cm), SiPM_Holder3_log, SiPM_Name, hydraulic_fluid_log, false,0);
     SiPM_PCB_Out3_phys = new G4PVPlacement(PlasticRot, G4ThreeVector(PCB_r*cos((i-1)*SiPM_Angle) ,PCB_r*sin((i-1)*SiPM_Angle), 46.44*cm), SiPM_PCB_Out3_log, SiPM_PCB_Name, hydraulic_fluid_log, false,0);
      }
      
      
    if(j==4)  //Fourth Row
      {      
      SiPM_Holder4_phys = new G4PVPlacement(SiPM1Rot, G4ThreeVector(SiPM_r*cos((i-1)*SiPM_Angle) ,SiPM_r*sin((i-1)*SiPM_Angle),42.02*cm), SiPM_Holder4_log, SiPM_Name, hydraulic_fluid_log, false,0);
     SiPM_PCB_Out4_phys = new G4PVPlacement(PlasticRot, G4ThreeVector(PCB_r*cos((i-1)*SiPM_Angle) ,PCB_r*sin((i-1)*SiPM_Angle), 42.02*cm), SiPM_PCB_Out4_log, SiPM_PCB_Name, hydraulic_fluid_log, false,0);
      }
     
    
      
    if(j==5)  //Fifth Row
      {      
     SiPM_Holder5_phys = new G4PVPlacement(SiPM5Rot, G4ThreeVector(SiPM5_r*cos((i-1)*SiPM_Angle) ,SiPM5_r*sin((i-1)*SiPM_Angle),62.35*cm), SiPM_Holder5_log, SiPM_Name, hydraulic_fluid_log, false,0);
     SiPM_PCB_Out5_phys = new G4PVPlacement(Plastic5Rot, G4ThreeVector(PCB5_r*cos((i-1)*SiPM_Angle) ,PCB5_r*sin((i-1)*SiPM_Angle), 62.35*cm), SiPM_PCB_Out5_log, SiPM_PCB_Name, hydraulic_fluid_log, false,0);
      }  
      
     
      
   
  
  SiPM_Name[0]='\0' ; 
  SiPM_PCB_Name[0]='\0' ;
  Piezo_Name[0]='\0' ;
  }
  
  }


 //==============================================================================================================================================
  //============================================================== New RTD Geometry ==============================================================
  //============================================================== Added by Gary S. ==============================================================
  //==============================================================================================================================================
  G4RotationMatrix* RTDRot = new G4RotationMatrix; 
  RTDRot->rotateZ(250.0*deg);
  RTDRot->rotateY(90.0*deg);
 
  G4RotationMatrix* RotY = new G4RotationMatrix;             
  RotY->rotateY(90.0*deg);

  //==================================== Define colours ===================================
  G4VisAttributes *copper= new G4VisAttributes(orange);
  copper->SetVisibility(true);
  copper->SetForceSolid(true);
  G4VisAttributes *plastic= new G4VisAttributes(lgrey);
  plastic->SetVisibility(true);
  plastic->SetForceSolid(true);

  //=======================================================================================
  //================================== RTD mother volume ==================================
  //=======================================================================================
  G4Box* RTD_mother_box = new G4Box("RTD_mother_box", 0.5*2.54*cm, 0.5*0.735*cm, 0.5*0.635*cm);
  RTD_mother_log  = new G4LogicalVolume(
    RTD_mother_box, 
    Epoxy_mat, 
    "RTD_mother_log"
  );

  RTD_mother_phys = new G4PVPlacement(RTDRot, 
    G4ThreeVector(cos(20 * PI / 180.0)*(12.7*cm + 0.5*0.625*cm) ,sin(20 * PI / 180.0)*(12.7*cm + 0.5*0.625*cm) ,55*cm ),
    //G4ThreeVector(cos(20 * PI / 180.0)*(12.7*mm + 0.5*0.625*mm) ,sin(20 * PI / 180.0)*(12.7*mm + 0.5*0.625*mm) ,5.5mm ),
    RTD_mother_log, 
    "RTD_mother1_phys", 
    hydraulic_fluid_log, 
    false,
    0
  );

  //=======================================================================================
  //==============================Create the copper baseplate==============================
  //=======================================================================================
  G4Box* RTD_block_box = new G4Box("RTD_block_box", 0.5*2.54*cm, 0.5*0.3175*cm, 0.5*0.635*cm);
     
  RTD_block_log = new G4LogicalVolume(
    RTD_block_box, 
    CuShield_mat, 
    "RTD_block_log"
  );

  RTD_block_phys = new G4PVPlacement(
    0, 
    G4ThreeVector(0.,0.5*4.175*mm,0.), 
    RTD_block_log, 
    "RTD_block_phys", 
    RTD_mother_log, 
    false,
    0
  ); 

  RTD_block_log->SetVisAttributes(copper);

  //=======================================================================================
  //==============================Create 3D plastic baseplate==============================
  //=======================================================================================
  G4Box *RTD_base = new G4Box("RTD_base", 0.5*2.54*cm, 0.5*1*mm, 0.5*0.635*cm);

  G4LogicalVolume *RTD_base_log = new G4LogicalVolume(
    RTD_base,
    Photopolymer_mat,
    "RTD_base_log"
  );
  
  G4VPhysicalVolume *RTD_base_phys = new G4PVPlacement(
    0,
    G4ThreeVector(0.,0.,0.),
    RTD_base_log,
    "RTD_base_phys",
    RTD_mother_log,
    false,
    0,
    true
  );

  RTD_base_log->SetVisAttributes(plastic);

  //=======================================================================================
  //===============================Create plastic wire clamp===============================
  //=======================================================================================
  G4Box *RTD_wire_clamp = new G4Box("RTD_wire_clamp", 0.5*0.508*cm, 0.5*3.175*mm, 0.5*0.635*cm);

  G4LogicalVolume *RTD_wire_clamp_log = new G4LogicalVolume(
    RTD_wire_clamp,
    Photopolymer_mat,
    "RTD_wire_clamp_log"
  );
  
  G4VPhysicalVolume *RTD_wire_clamp_phys = new G4PVPlacement(
    0,
    G4ThreeVector(0.5*17.7*mm,-0.5*4.175*mm,0.),
    RTD_wire_clamp_log,
    "RTD_wire_clamp_phys",
    RTD_mother_log,
    false,
    0,
    true
  );

  RTD_wire_clamp_log->SetVisAttributes(plastic); 

  //=======================================================================================
  //================================Create connector mount=================================
  //=======================================================================================
  G4Box *RTD_connector_mount = new G4Box("RTD_connector_mount", 0.5*0.3*cm, 0.5*3.175*mm, 0.5*0.635*cm);

  //=======================================================================================
  //==================================Create connectors====================================
  //=======================================================================================
  G4Tubs *RTD_connector = new G4Tubs("RTD_connector", 0., 0.5*0.254*mm, 0.5*10.5*mm, 0.*deg, 360.*deg);
  G4VSolid *RTD_connector_1_2 = new G4UnionSolid("RTD_connector_1_2", RTD_connector, RTD_connector, 0, G4ThreeVector(0.5*3.*mm, 0., 0.));
  G4VSolid *RTD_connector_3_4 = new G4UnionSolid("RTD_connector_3_4", RTD_connector_1_2, RTD_connector_1_2, 0, G4ThreeVector(2*0.5*3.*mm, 0., 0.));

  G4LogicalVolume *RTD_connector_1_4_log = new G4LogicalVolume(
    RTD_connector_3_4, 
    brass_mat, 
    "RTD_connector_1_4_log"
  );

  G4VPhysicalVolume *RTD_connector_1_4_phys = new G4PVPlacement(
    RotY, 
    G4ThreeVector(-0.5*11.2*mm, -0.5*4.175*mm, -0.5*4.5*mm), 
    RTD_connector_1_4_log, 
    "RTD_connector_1_4_phys", 
    RTD_mother_log, 
    false, 
    0, 
    true
  );

  RTD_connector_1_4_log->SetVisAttributes(copper);

  //=======================================================================================
  //=============================Create connector feedthrough==============================
  //=======================================================================================
  G4VSolid *RTD_connector_feedthru = new G4SubtractionSolid(
    "connector_feedthru", 
    RTD_connector_mount, 
    RTD_connector_3_4, 
    RotY, 
    G4ThreeVector(0.,0.,-0.5*4.5*mm)
  );

  G4LogicalVolume *RTD_connector_feedthru_log = new G4LogicalVolume(
    RTD_connector_feedthru, 
    Photopolymer_mat, 
    "RTD_connector_feedthru_log"
  );

  G4VPhysicalVolume *RTD_connector_feedthru_phys = new G4PVPlacement(
    0, 
    G4ThreeVector(-0.5*11.2*mm, -0.5*4.175*mm, 0.), 
    RTD_connector_feedthru_log,
    "RTD_connector_feedthru_phys", 
    RTD_mother_log, 
    false, 
    0, 
    true
  );

  RTD_connector_feedthru_log->SetVisAttributes(plastic);
  
  // Create more RTDs
  RTD_mother_phys = new G4PVPlacement(RTDRot, G4ThreeVector(cos(20 * PI / 180.0)*(12.7*cm + 0.5*0.625*cm) ,sin(20 * PI / 180.0)*(12.7*cm + 0.5*0.625*cm) ,48.5*cm ), RTD_mother_log, "RTD_mother2_phys", hydraulic_fluid_log, false,0);
  RTD_mother_phys = new G4PVPlacement(RTDRot, G4ThreeVector(cos(20 * PI / 180.0)*(12.7*cm + 0.5*0.625*cm) ,sin(20 * PI / 180.0)*(12.7*cm + 0.5*0.625*cm) ,42*cm ), RTD_mother_log, "RTD_mother3_phys", hydraulic_fluid_log, false,0);
  


  //// SiPMs Cables, located netx to the SiPM panel
                   
    G4Tubs* RTD_Cable_Kapton = new G4Tubs("RTD_Cable_Kapton",0, 0.5*2.54*0.05*cm, 0.5*20.0*cm, 0.*deg, 360.*deg);
  
   RTD_Cable_Kapton_log  = new G4LogicalVolume(RTD_Cable_Kapton, sapphire_mat, "RTD_Cable_Kapton_log"); //Kapton_mat
   RTD_Cable_Kapton_phys = new G4PVPlacement(0, G4ThreeVector(13.27*cm ,0.0 ,49.4*cm ), RTD_Cable_Kapton_log, "RTD_Cable_Kapton_phys", hydraulic_fluid_log, false,0);
  
  
   G4Tubs* RTD_Cable_Cu = new G4Tubs("RTD_Cable_Cu",0.5*2.54*0.029*cm, 0.5*2.54*0.041*cm, 0.5*20.0*cm, 0.*deg, 360.*deg);
  
   RTD_Cable_Cu_log  = new G4LogicalVolume(RTD_Cable_Cu, CuShield_mat, "RTD_Cable_Cu_log");
   RTD_Cable_Cu_phys = new G4PVPlacement(0, G4ThreeVector(0.0 ,0.0 ,0.0 ), RTD_Cable_Cu_log, "RTD_Cable_Cu_phys", RTD_Cable_Kapton_log, false,0);
  
  
 
   G4Tubs* RTD_Cable_Cu_1 = new G4Tubs("RTD_Cable_Cu_1",0, 0.5*2.54*0.0105*cm, 0.5*20.0*cm, 0.*deg, 360.*deg);
  
   RTD_Cable_Cu_1_log  = new G4LogicalVolume(RTD_Cable_Cu_1, CuShield_mat, "RTD_Cable_Cu_1_log");
   RTD_Cable_Cu_1_phys = new G4PVPlacement(0, G4ThreeVector(0.0 ,0.0 ,0.0 ), RTD_Cable_Cu_1_log, "RTD_Cable_Cu_1_phys", RTD_Cable_Kapton_log, false,0);



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
 if(LAr_log)                                                               
      SetSensitiveDetector(LAr_log,LXeSD.Get());  

/*
     if (pmtSD.Get() == 0)                                        //Aquí detecto los eventos en el SiPM
    {
      G4String name="/DMXDet/pmtSD";
      DMXPmtSD* aSD = new DMXPmtSD(name);
      pmtSD.Put(aSD);
    }
  G4SDManager::GetSDMpointer()->AddNewDetector(pmtSD.Get()); 
  if (SIPM_Si_log)
    SetSensitiveDetector(SIPM_Si_log,pmtSD.Get());                           

*/

    return;

}

