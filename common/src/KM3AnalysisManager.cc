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

#include "KM3AnalysisManager.hh"
#include "G4UnitsTable.hh"
#include "G4VTouchable.hh"
#include <G4NavigationHistory.hh>
#include <G4TransportationManager.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

namespace km3net
{
  string KM3AnalysisManager::fRootFileName = "default.root";
  bool   KM3AnalysisManager::fKeepAll = false;
  bool   KM3AnalysisManager::fFastAnalysisOnly = true;
  short  KM3AnalysisManager::fVerbose = 0;

  KM3AnalysisManager* KM3AnalysisManager::fManager = 0;

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

  KM3AnalysisManager* KM3AnalysisManager::getIt()
  {
    if(!fManager) {
      throw ("uncorrect use: a daugther of KM3AnalysisManager should be initializated first !!");
    }
    return fManager;
  }


  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

  KM3AnalysisManager::KM3AnalysisManager()
  {

  }

  KM3AnalysisManager::~KM3AnalysisManager()
  {

  }

  void KM3AnalysisManager::Init (string RFNAME)
  {
    if (RFNAME.size ())
      fRootFileName=RFNAME;

    fRFile = new TFile (fRootFileName.data (), "RECREATE");

    fTRunInfos = new TTree("RunInfo","Infos about the simulation.");
    fTRunInfos->Branch("NBTrigg",&fNbTrigg,"NbTrigg/l");
    fTEvent = new TTree("Hit","Hits in a TTree");
    fTEvent->Branch("Vertex",&fVertexX,"VertexX:VertexY:VertexZ:VertexT:VertexP:VertexDirX:VertexDirY:VertexDirZ:VertexEnergy:AbsOMHitX:AbsOMHitY:AbsOMHitZ:OMHitT:OMHitP:PhotoX:PhotoY:PhotoZ:PhotoT:PhotoP:OMID/S:PMID");
    fTFastAnalysis = new TTree("FastAnalysis","Data used for very fast analysis.");
    fTFastAnalysis->Branch("PhotonReach",&fNbOM,"NbOM/l:NbPhoto");

  }


  bool KM3AnalysisManager::HasBeenInTheOMYet ()
  {
    if (fOMHitX==0 &&
        fOMHitY==0 &&
        fOMHitZ==0)
      return false;

    return true;
  }

  void KM3AnalysisManager::AddPhotoTubeEntranceInfo (const G4Step* STEP)
  {
    if (fNbPhoto % 1000 ==0)
      {
        cout << "\r\033[20C" << fNbPhoto << " in pmt ";
        cout.flush ();
      }
    fNbPhoto++;
    if (fRFile == NULL)
      return;
    if (fFastAnalysisOnly)
      return;

    fPMID=STEP->GetPostStepPoint ()->GetTouchable ()->GetCopyNumber (1);

    const G4StepPoint* prestep = STEP->GetPreStepPoint();
    const G4ThreeVector& globalPoint = prestep->GetPosition();
    const G4VTouchable*   touch = prestep->GetTouchable();
    const G4ThreeVector& localpoint=touch->GetHistory ()->GetTopTransform().TransformPoint(globalPoint);
    const G4ThreeVector norm=prestep->GetPhysicalVolume()->GetLogicalVolume ()->GetSolid ()->SurfaceNormal( localpoint );

    fPhotoX=localpoint.x ();
    fPhotoY=localpoint.y ();
    fPhotoZ=localpoint.z ();
    fPhotoT=localpoint.theta ();
    fPhotoP=localpoint.phi ();

  }

  void KM3AnalysisManager::AddOMEntranceInfo (const G4Step* STEP)
  {


    const G4StepPoint* prestep = STEP->GetPreStepPoint();
    G4String prestepName=prestep->GetPhysicalVolume()->GetName();

    if (prestepName != "Target")
      {
        return;
      }

    if (fOMBasePrefix.size () == 0)
      fOMBasePrefix=KM3DetectorConstruction::fOMBasePrefix;

    const G4StepPoint* poststep = STEP->GetPostStepPoint();
    G4String poststepName=poststep->GetPhysicalVolume()->GetName();

    if (fNbOM % 1000 ==0)
      {
        cout << "\r\033[40C" << fNbOM << " in OM ";
        cout.flush ();
      }
    fNbOM++;
    if (fRFile == NULL)
      return;
    if (fFastAnalysisOnly)
      return;

    const G4ThreeVector& globalPoint = poststep->GetPosition();
    const G4VTouchable*   touch = poststep->GetTouchable();
    const G4ThreeVector& localpoint=touch->GetHistory ()->GetTopTransform().TransformPoint(globalPoint);

    fOMId = touch->GetVolume ()->GetCopyNo ();
    fOMHitX=globalPoint.x ();
    fOMHitY=globalPoint.y ();
    fOMHitZ=globalPoint.z ();
    fOMHitT = localpoint.theta ();
    fOMHitP = localpoint.phi ();

  }

  void KM3AnalysisManager::SetOriginOfEvent (const G4Event* EV)
  {
    if (fRFile == NULL)
      return;

    fNbTrigg++;

    if (fNbTrigg % 1000 == 0)
      cout  <<"\r"<<fNbTrigg << " trigg" << flush;

    if (fFastAnalysisOnly)
      return;

    const G4ThreeVector& vertexVect= EV->GetPrimaryVertex ()->GetPosition ();
    const G4ThreeVector& vertexDir = EV->GetPrimaryVertex ()->GetPrimary ()->GetMomentumDirection ();

    fVertexX = vertexVect.x ();
    fVertexY = vertexVect.y ();
    fVertexZ = vertexVect.z ();
    fVertexT = vertexVect.theta ();
    fVertexP = vertexVect.phi ();

    fVertexDirX = vertexDir.x ();
    fVertexDirY = vertexDir.y ();
    fVertexDirZ = vertexDir.z ();

    fVertexEnergy = EV->GetPrimaryVertex ()->GetPrimary ()->GetKineticEnergy ();
  }

  void KM3AnalysisManager::SetEnergy (double E)
  {
    if (fVerbose > 0)
      {
        cout << "KM3AnalysisManager::SetEnergy: nothing to do with the Energy " << E << endl;
      }

  }

  void KM3AnalysisManager::EndOfEvent ()
  {
    fVertexX=0;
    fVertexY=0;
    fVertexZ=0;
    fVertexT=0;
    fVertexP=0;
    fOMHitX=0;
    fOMHitY=0;
    fOMHitZ=0;
    fOMHitT=0;
    fOMHitP=0;
    fPhotoX=0;
    fPhotoY=0;
    fPhotoZ=0;
    fPhotoT=0;
    fPhotoP=0;
    fOMId=-1;
    fPMID=0;
  }

}

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
