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

#include "G4ios.hh"

#include "KM3DetectorConstruction.hh"
#include "KM3SteppingAction.hh"
#include "KM3AnalysisManager.hh"
#include <KM3OpticalModuleMgr.hh>
#include <KM3TrackInformation.hh>
#include "G4Track.hh"
#include "globals.hh"
#include "G4SteppingManager.hh"
#include "G4VTouchable.hh"

#include <KM3PMTOpticalModel.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalBorderSurface.hh>
#include <G4LogicalSkinSurface.hh>

#include <G4TransportationManager.hh>

#include <iostream>
#include <fstream>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace CLHEP;

namespace km3net
{

  KM3SteppingAction::KM3SteppingAction()
  {

  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

  KM3SteppingAction::~KM3SteppingAction()
  { }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

  void KM3SteppingAction::UserSteppingAction(const G4Step* fStep)
  {
    //#ifdef G4ANALYSIS_USE
    G4Track* fTrack = fStep->GetTrack();
    G4StepPoint* prestep = fStep->GetPreStepPoint();
    G4StepPoint* poststep = fStep->GetPostStepPoint();

    G4VPhysicalVolume* thePrePV=fStep->GetPreStepPoint()->GetPhysicalVolume();
    G4VPhysicalVolume* thePostPV=fStep->GetPostStepPoint()->GetPhysicalVolume();


    //moved here temporally by Vladimir since it seems that for the photons "info" is not created
    if (thePostPV != NULL) {
      //if (KM3PMTOpticalModel::IsAPhotocathode (fStep))   //VLA: It never enters here!!!!
      {
        if (fTrack->GetDefinition()->GetParticleType() == "opticalphoton" &&
            fTrack->GetTrackStatus() != fStopAndKill)
          {
            //if (thePostPV->GetName ()=="NOREFLEXVaccumInPmt")
            string vac = thePostPV->GetName ();
            if (vac.find("NOREFLEXVaccumInPmt")!=std::string::npos)
              {
                fTrack->SetTrackStatus(fStopAndKill);
                return;
              }
          }
      }
    }



    KM3TrackInformation* info = (KM3TrackInformation*)(fStep->GetTrack()->GetUserInformation());

    if (info == 0) return;

    //#ifdef G4ANALYSIS_USE
    /*
      G4Track* fTrack = fStep->GetTrack();
      G4StepPoint* prestep = fStep->GetPreStepPoint();
      G4StepPoint* poststep = fStep->GetPostStepPoint();
    */


    if (fTrack->GetGlobalTime ()>1e10) //remove K40 lifetime
      {
        fTrack->SetGlobalTime (fTrack->GetLocalTime ());
        prestep->SetGlobalTime(prestep->GetLocalTime());
        poststep->SetGlobalTime(poststep->GetLocalTime());
      }

    //G4VPhysicalVolume* thePrePV=fStep->GetPreStepPoint()->GetPhysicalVolume();
    //G4VPhysicalVolume* thePostPV=fStep->GetPostStepPoint()->GetPhysicalVolume();

    G4VPhysicalVolume* thePreOriginalPV = info->GetOriginalPrePV ();
    G4VPhysicalVolume* thePostOriginalPV = info->GetOriginalPostPV ();


    if (!thePrePV || !thePostPV || !thePreOriginalPV || !thePostOriginalPV )
      return;

    G4String prestepName=thePreOriginalPV->GetLogicalVolume ()->GetName();
    G4String poststepName=thePostOriginalPV->GetLogicalVolume ()->GetName();

    if (KM3PMTOpticalModel::IsAPhotocathode (thePreOriginalPV, thePostOriginalPV) &&
        fTrack->GetDefinition()->GetParticleName() == "e-" &&
        info->GetOriginalParticleType () == "opticalphoton" &&
        //prestepName.find (KM3PhotoTubeMgr::fPhotoName) != string::npos) 
        /*
         * VK: this does not work for the photons going back from the PMT vacuume (if first they were transmitted inside PMT
         * So, this should be used:
         * prestepName.find (KM3PhotoTubeMgr::fPhotoName) != string::npos) || prestepName.find (KM3PhotoTubeMgr::fPhotoName) != string::npos)  
         * We decided not to do this check anyway (since it is enough that parent photon was in the photocathode.
         * Instead, we add below the check that the electron should be in the PMT sensitive area (PMT vacuum).
         */
        fStep->GetPreStepPoint()->GetMaterial()->
                        GetMaterialPropertiesTable ()->ConstPropertyExists("electronSensitive"))

      {
        KM3AnalysisManager::getIt()->AddPhotoTubeEntranceInfo (fStep);
      }
    else //bravo
      return;

    G4LogicalSurface* Surface =
      G4LogicalBorderSurface::GetSurface(thePreOriginalPV, thePostOriginalPV);
    if (!Surface)
      Surface = G4LogicalBorderSurface::GetSurface(thePostPV, thePrePV);
    if (!Surface)
      Surface = G4LogicalSkinSurface::GetSurface(thePostPV->GetLogicalVolume ());
    if (!Surface)
      Surface = G4LogicalSkinSurface::GetSurface(thePrePV->GetLogicalVolume ());

    if (!Surface ||
        !LoadMaterialProps (dynamic_cast <G4OpticalSurface*> (Surface->GetSurfaceProperty())
                            ->GetMaterialPropertiesTable()) )
      return;

    fTrack->SetTrackStatus(fStopAndKill); // kill e-

    fEnergy = info->GetOriginalEnergy ();

    G4ThreeVector globalPoint = fStep->GetPreStepPoint()->GetPosition ();
    const G4VTouchable*   touch = fStep->GetPostStepPoint()->GetTouchable();
    unsigned char detectedGamma = CheckDetection (touch->GetHistory()->GetTopTransform().TransformPoint(globalPoint));
    if (detectedGamma != 0)
      KM3AnalysisManager::getIt()->AddDetectedParticule (info, detectedGamma);
  }




  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

  bool KM3SteppingAction::LoadMaterialProps (G4MaterialPropertiesTable* aMaterialPropertiesTable)
  {


    if (!aMaterialPropertiesTable) {
      return false;
    }
    //angular efficiency is an added propertyfor some additional calibration in addition to the thickness.
    _angular_efficiency_err= aMaterialPropertiesTable->GetProperty("ANGULAR_EFFICIENCY_ERROR");
    if (_angular_efficiency_err == NULL)
      {
        return false;
      }

    _angular_efficiency= aMaterialPropertiesTable->GetProperty("ANGULAR_EFFICIENCY");
    if (_angular_efficiency == NULL)
      {
        return false;
      }

    _efficiency_photocathode= aMaterialPropertiesTable->GetProperty("EFFICIENCY");
    if (_efficiency_photocathode == NULL)
      {
        return false;
      }

    _thickness_photocathode= aMaterialPropertiesTable->GetProperty("THICKNESS");
    if (_thickness_photocathode == NULL)
      {
        return false;
      }

    G4MaterialPropertyVector* kindex_photocathode= aMaterialPropertiesTable->GetProperty("KINDEX");
    if (kindex_photocathode == NULL)
      {
        return false;
      }

    G4MaterialPropertyVector* rindex_photocathode= aMaterialPropertiesTable->GetProperty("RINDEX");
    if (rindex_photocathode == NULL)
      {
        return false;
      }

    if (aMaterialPropertiesTable->ConstPropertyExists ("ElectronAbsLength")) fElectronAbsLength = aMaterialPropertiesTable->GetConstProperty ("ElectronAbsLength")*nm;

    /*aMaterialPropertiesTable->DumpTable ();
      exit (1);
    */
    return true;

  }

  //adding angular efficiency here
  unsigned char KM3SteppingAction::CheckDetection(const G4ThreeVector& localPoint)
  {
    unsigned char detectedGamma = 4;
    //double eta_E = 0.5; //probability that the electron goes in or out the photocathode
    //from the photonis, can be improved. The abs length is ~30nm, the electron is produced at about an half of the photocathode.
    double teta = localPoint.theta ();//acos(costeta);

    double eta_teta=_angular_efficiency->Value(teta/3*2.1);

    double eta_teta_err=_angular_efficiency_err->Value(teta/3*2.1);

    //cout << "DEBUG: eta_teta " << eta_teta  << " +- " << eta_teta_err << ' ' << teta/3*2.1 <<' ' << endl;

    double rndm = G4UniformRand();
    if(rndm < eta_teta*(1.-eta_teta_err))
      {
        detectedGamma = 3;
      }
    else if(rndm < eta_teta)
      {
        detectedGamma = 2;
      }
    else if (rndm < eta_teta*(1+eta_teta_err))
      {
        detectedGamma = 1;
      }
    //}

    return detectedGamma;
  }

}
