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

#include <KM3PhotoTubeMgr.hh>
#include "KM3DetectorConstruction.hh"
#include "XPPhSteppingAction.hh"
#include "KM3AnalysisManager.hh"
#include "G4Track.hh"
#include "globals.hh"
#include "G4SteppingManager.hh"
#include "G4VTouchable.hh"

#include <iostream>
#include <fstream>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace CLHEP;

namespace km3net
{

  XPPhSteppingAction::XPPhSteppingAction()
  {
    int unsize = sizeof (fAbsorptions[0]);
    memset (fAbsorptions,0,91*unsize);
    memset (fReflections,0,91*unsize);
    memset (fTransmissions,0,91*unsize);
    memset (fCounter,0,91*sizeof(fCounter[0]));

  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

  XPPhSteppingAction::~XPPhSteppingAction()
  {

    ofstream ofs ("XPPh.data");
    ofs << "angle absorption transmission reflection" << endl;
    for (unsigned i=0; i<91; i++)
      {
        unsigned n=fCounter[i];
        if (n == 0) n=1;
        ofs << i << ' ' << fAbsorptions [i]/n << ' ' << fTransmissions[i]/n << ' ' << fReflections[i]/n << endl;
      }
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

  void XPPhSteppingAction::UserSteppingAction(const G4Step* fStep)
  {
    G4StepPoint* poststep = fStep->GetPostStepPoint();
    G4StepPoint* prestep = fStep->GetPreStepPoint();
    G4VPhysicalVolume* poststeVol=poststep->GetPhysicalVolume();
    if (poststeVol == 0)
      return;

    string prestepname  = prestep->GetPhysicalVolume()->GetName ();
    string poststepname = poststeVol->GetName ();

    /*    if (poststepname == OpticalModuleMaker::fPhotoName)
      {
        cout << "The propagation length in " << poststepname << " for " << prestep->GetKineticEnergy()/eV << " eV is "
             << poststeVol->GetLogicalVolume ()->GetMaterial ()->GetMaterialPropertiesTable ()->GetProperty ("KINDEX")->Value (prestep->GetKineticEnergy ())/nm
             << " nm " << endl;
             }
    */

    G4Track* fTrack = fStep->GetTrack();
    if (prestep->GetGlobalTime () == 0)
      {
        fCurrentTrackId = fTrack->GetTrackID ();
        fCurrentAngle=poststep->GetMomentumDirection ().angle (G4ThreeVector(0,0,-1))/deg;
        if (fCurrentAngle > 90)
          {
            fTrack->SetTrackStatus (fStopAndKill);
            return;
          }
        fCounter[unsigned(fCurrentAngle+0.5)]++;

      }

    if (prestepname == KM3PhotoTubeMgr::fPhotoName || poststepname == "VacuumInPmt" )
      {
        if (fTrack->GetTrackStatus () == fStopAndKill)
          {
            fAbsorptions[unsigned(fCurrentAngle+0.5)]++;
            return;
          }
      }

    if (fTrack->GetTrackStatus () == fStopAndKill)
      return;

    if ( poststepname == "VacuumInPmt" && poststep->GetMomentumDirection ().z() <0 )
      {
        fTransmissions[unsigned(fCurrentAngle+0.5)]++;
        fTrack->SetTrackStatus (fStopAndKill);
        return;
      }
    if ( poststepname == "OpticalGel" && poststep->GetMomentumDirection ().z() > 0 )
      {
        fReflections[unsigned(fCurrentAngle+0.5)]++;
        fTrack->SetTrackStatus (fStopAndKill);
        return;
      }

    return;


  }
}

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
