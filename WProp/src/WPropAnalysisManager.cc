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


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TRandom.h>
#include <TMath.h>

#include <WPropAnalysisManager.hh>
#include "G4UnitsTable.hh"
#include "G4VTouchable.hh"
#include <G4NavigationHistory.hh>
#include <G4Event.hh>
#include <KM3DetectorConstruction.hh>
#include <G4SystemOfUnits.hh>

#include <fstream>
#include <string>
#include <iostream>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

namespace km3net
{
  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

  void LoadEfficiency (string file, vector<float>& toreturn)
  {
    toreturn.clear ();

    ifstream f(file.data ());

    int angle=0;
    float read;
    float maxAA=0;

    while (angle < 180 && f >> read >> read)
      {
        angle++;
        toreturn.push_back (read);
        if (read>maxAA)
          maxAA=read;
      }

    for (int i=181; i<360; i++)
      {
        toreturn.push_back (toreturn[360-i]);
      }

    for (int i=0; i<360; i++)
      toreturn[i]/=maxAA;
  }
  WPropAnalysisManager* WPropAnalysisManager::getIt()
  {
    if(!fManager) {
      fManager = new WPropAnalysisManager();
    }
    return (WPropAnalysisManager*)fManager;
  }


  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

  WPropAnalysisManager::WPropAnalysisManager()
  {
    fLEDDetected[0] = 0;
    fLEDDetected[1] = 1;
    fRandom=new TRandom (time (0));
  }

  void WPropAnalysisManager::Init (string RFNAME)
  {
    KM3AnalysisManager::Init (RFNAME);
    //fTFastAnalysis = new TTree("FastAnalysis","Data used for very fast analysis.");
    //fTFastAnalysis->Branch("FastAnalysis",&fNbOM,"NbOM/l:NbPhoto:NbDetected[3][100]/i:Coincidences[7][3]");
    //fTEvent->Branch("Hit",&fOMHitX,"OMHitX:OMHitY:OMHitZ:OMHitT:OMHitP:PhotoX:PhotoY:PhotoZ:PhotoT:PhotoP:TotalTime:Energy:OMID/S:Limit/b");
    fTEvent->Branch("SphereDetection",&fEmmitingDirX,"EmmitingDirX/f:EmmitingDirY:EmmitingDirZ:FloorCrossPosX:FloorCrossPosY:FloorCrossPosZ:FloorCrossDirX:FloorCrossDirY:FloorCrossDirZ:TotalTime:Energy:WeigthANTARES:WeigthNEMO:ESWeight/i:KopelWeight:FloorNbCross/b:LEDDetected[2]/l");


    LoadEfficiency (string(PROJECT_SOURCE_DIR)+"/WProp/resources/AANemo.dat", fEffNEMO);

    LoadEfficiency (string(PROJECT_SOURCE_DIR)+"/WProp/resources/AAAntares.dat", fEffANTARES);

    string ledtotake;
    ifstream areader(string(PROJECT_SOURCE_DIR)+"/resources/ledselect");  /////to clean
    areader >> ledtotake;
    LoadEfficiency (string(PROJECT_SOURCE_DIR)+"/resources/" + ledtotake, fLEDAngle);

    fhLED=new TH1F("fhLED", "fhLED",fLEDAngle.size (),0, fLEDAngle.size ());
    fhLED->Adopt (fLEDAngle.size (), fLEDAngle.data ());

  }


  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

  WPropAnalysisManager::~WPropAnalysisManager()
  {

  }

  void WPropAnalysisManager::Write  ()
  {
    cout << endl;
    if (fRFile == NULL)
      return;

    fTRunInfos->Fill ();
    fTFastAnalysis->Fill ();

    fRFile->Write ();
    fRFile->Close ();
    delete (fRFile);
    cout << "Analysis file was normaly closed." << endl;
  }

  void WPropAnalysisManager::SetOriginOfEvent (const G4Event* EV)
  {
    KM3AnalysisManager::SetOriginOfEvent (EV);

    fRealEmmiting = EV->GetPrimaryVertex ()->GetPrimary ()->GetMomentum ();
    fRealEmmiting.setTheta (fRealEmmiting.getTheta ()+5*deg); ///////warning!!!/////////////////

    G4ThreeVector globalDir = EV->GetPrimaryVertex ()->GetPrimary ()->GetMomentum ();;

    /*globalDir.setRThetaPhi(1,
                           fhLED->GetRandom ()*deg,
                           twopi*G4iformRand());
    */
    fEmmitingDirX=globalDir.x ();
    fEmmitingDirY=globalDir.y ();
    fEmmitingDirZ=globalDir.z ();

  }

  void WPropAnalysisManager::IncrementProcessWeight (KM3PetzoldScattering::PROCESS PR)
  {
    switch (PR)
      {
      case KM3PetzoldScattering::EINSTEINSMOLUCHOWSKI :
        fESWeight++;
        break;
      case KM3PetzoldScattering::KOPELEVITCH :
        fKopelWeight++;
        break;
      default :
        break;
      }
  }

  void WPropAnalysisManager::AddDetectedParticule (const KM3TrackInformation* INFO, unsigned char LIM)
  {

    fEnergy=INFO->GetOriginalEnergy ();
    if (fRFile == NULL)
      return;
    if (fFastAnalysisOnly)
      return;

    fFloorNbCross=LIM;

    const G4ThreeVector& globalDir = INFO->GetOriginalMomentum ();
    const G4ThreeVector& globalPoint = INFO->GetOriginalPosition ();
    const G4ThreeVector& posVect=INFO->GetOriginalLocalPosition ();

    fTotalTime = INFO->GetOriginalTime ();


    fFloorCrossPosX = posVect.x (); //\theta_s
    fFloorCrossPosY = posVect.y ();
    fFloorCrossPosZ = posVect.z ();

    fFloorCrossDirX=globalDir.x ();//\theta_i
    fFloorCrossDirY=globalDir.y ();
    fFloorCrossDirZ=globalDir.z ();

    float theta_s = globalPoint.angle(fRealEmmiting)/deg;

    unsigned theta_i = unsigned(globalPoint.angle (globalDir)/deg+0.5);

    fWeigthANTARES = fEffANTARES[theta_i];
    fWeigthNEMO    = fEffNEMO[theta_i];

    //2Floor
    if (fWeigthNEMO < 1e-8 || fWeigthANTARES < 1e-8)
      return;

    bool hasAnEvent=false;

    size_t base=1;

    for (unsigned LEDit=1; LEDit<120; LEDit++)
      {
        float LEDProb = TMath::Gaus(theta_s,0,LEDit/10.);

        if (LEDProb < 1e-8)
          continue;
        float rand=fRandom->Uniform ();

        if (rand < LEDProb*fWeigthNEMO)
          {
            if (LEDit < 64)
              fLEDDetected[0] |= base << LEDit;
            else
              fLEDDetected[1] |= base << (LEDit-64);

            hasAnEvent=true;
          }
      }
    if (hasAnEvent)
      {
        fTEvent->Fill ();
      }

  }

  void WPropAnalysisManager::EndOfEvent ()
  {
    if (fRFile == NULL)
      return;

    if (fFastAnalysisOnly)
      return;

    fLEDDetected[0]=0;
    fLEDDetected[1]=0;

    fEmmitingDirX =0;
    fEmmitingDirY  =0;
    fEmmitingDirZ  =0;
    fFloorCrossPosX=0;
    fFloorCrossPosY=0;
    fFloorCrossPosZ=0;
    fFloorCrossDirX=0;
    fFloorCrossDirY=0;
    fFloorCrossDirZ=0;
    fTotalTime = 0;
    fEnergy=0;
    fFloorNbCross=0;
    fESWeight=0;
    fKopelWeight=0;


  }

}

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
