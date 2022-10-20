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
// ScintSD (scintillator sensitive detector definition) program
// --------------------------------------------------------------

#include "DMXScintSD.hh"

#include "DMXScintHit.hh"
#include "DMXDetectorConstruction.hh"

#include "G4VPhysicalVolume.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Ions.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"//

#include "G4OpticalPhoton.hh"

#include <iostream>
#include <fstream>
#include <cstring>
#include "G4UnitsTable.hh"





int n=1;
int HasHit=0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DMXScintSD::DMXScintSD(G4String name) 
  :G4VSensitiveDetector(name)
{
  G4String HCname="scintillatorCollection";
  collectionName.insert(HCname);

  Info.open("Informacion.csv");
  Info1.open("General.csv");
  Info << "Event,"<<"Energy_Cinetica," << "Particle," << "x,"<<"y,"<<"z,"<<"Volume"<<'\n';
  Info1 << "Evento," << "Hits_Generados\n";




}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DMXScintSD::~DMXScintSD(){ Info.close();  Info1.close();}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXScintSD::Initialize(G4HCofThisEvent*)
{
  scintillatorCollection = new DMXScintHitsCollection
    (SensitiveDetectorName,collectionName[0]);

  HitID = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool DMXScintSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{

  //need to know if this is an optical photon and exclude it:
  if(aStep->GetTrack()->GetDefinition()
    == G4OpticalPhoton::OpticalPhotonDefinition()) return false;
  if(aStep->GetTrack()->GetDefinition()
    == G4Electron::ElectronDefinition()) 

    {aStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
    return false;}
  if(aStep->GetTrack()->GetDefinition()
    == G4Proton::ProtonDefinition()) 

    {aStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
    return false;}

  
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4double ek = aStep->GetPreStepPoint()->GetKineticEnergy();
  G4ParticleDefinition* particleType = aStep->GetTrack()->GetDefinition();
  G4String particleName = particleType->GetParticleName();

  G4double posx = aStep->GetPostStepPoint()->GetPosition().x();
  G4double posy = aStep->GetPostStepPoint()->GetPosition().y();
  G4double posz = aStep->GetPostStepPoint()->GetPosition().z();

  G4String Volume = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();
  G4bool FirstStep = aStep-> IsFirstStepInVolume();
  G4bool LastStep = aStep-> IsLastStepInVolume();
  G4double StepLength = aStep->GetStepLength();
  G4bool hit = false;
  


  char Ar36 [10]="Ar36";
  char Ar38 [10]="Ar38";
  char Ar40 [10]="Ar40";

  G4double stepl = 0.;
  if (particleType->GetPDGCharge() != 0.)
    stepl = aStep->GetStepLength();
  
  //if ((edep==0.)&&(stepl==0.)) return false;      


  // fill in hit
  DMXScintHit* newHit = new DMXScintHit();
  newHit->SetEdep(edep);
  newHit->SetPos(aStep->GetPostStepPoint()->GetPosition());
  newHit->SetTime(aStep->GetPreStepPoint()->GetGlobalTime());
  newHit->SetParticle(particleName);
  newHit->SetParticleEnergy(aStep->GetPreStepPoint()->GetKineticEnergy() );

  HitID = scintillatorCollection->insert(newHit);

 //if(strcmp(particleName,"proton") == 0 || strcmp(particleName,"e-") == 0 || strcmp(particleName,Ar40) == 0)
 //if (FirstStep == true || LastStep ==true)
//if (FirstStep == true)
 


  if(Volume == "physS"){
  HasHit = HasHit+1;}

/*  if(Volume == "physS" && LastStep == true){
  Info << n <<",";}
  if(Volume == "OutOfWorld" && hit == true){
  Info << "1" <<'\n' ;}
  if(Volume =="OutOfWorld" && hit == false){
  Info << "0" <<'\n' ;}*/
 
  if(Volume == "physS"){
  Info << '\n' << n <<"," << 1000000*ek <<","<< HasHit;}
  //if(Volume == "physS" && LastStep == true){
  //Info << '\n' << n <<"," << ek;}
  if(Volume == "physS" && HasHit == 1){
  Info << "1" ;}
   
  //if(Volume == "World"){
  //Info << '\n';//}

  return true;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXScintSD::EndOfEvent(G4HCofThisEvent* HCE)
{

  G4String HCname = collectionName[0];
  static G4int HCID = -1;
  if(HCID<0)
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(HCname);
  HCE->AddHitsCollection(HCID,scintillatorCollection);

  G4int nHits = scintillatorCollection->entries();
  if (verboseLevel>=1){

    G4cout << "     LXe collection: " <<  nHits << " hits" << G4endl;
    Info1 << n<<","<<nHits <<"\n";
    n++;

  }
  if (verboseLevel>=2)
    scintillatorCollection->PrintAllHits();



}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXScintSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXScintSD::DrawAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXScintSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


