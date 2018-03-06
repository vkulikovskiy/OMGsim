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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include <KM3ExperimentTypeMessenger.hh>
#include <KM3DetectorConstruction.hh>

#include "WPropAnalysisManager.hh"
#include "WPropSteppingAction.hh"

#include "K40AnalysisManager.hh"

#include "AAAnalysisManager.hh"
#include <XPPhSteppingAction.hh>

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "globals.hh"

#include <KM3RunManager.hh>
#include <KM3PrimaryGeneratorAction.hh>
#include <KM3PhysicsList.hh>
#include <KM3RunAction.hh>
#include <KM3EventAction.hh>
#include <KM3TrackingAction.hh>

#include <iostream>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace km3net
{
  KM3ExperimentTypeMessenger* KM3ExperimentTypeMessenger::fInstance=0;
  KM3ExperimentTypeMessenger* KM3ExperimentTypeMessenger::GetIt () {
    if (fInstance == 0)
      fInstance = new KM3ExperimentTypeMessenger;
    return fInstance;
  }

  KM3ExperimentTypeMessenger::KM3ExperimentTypeMessenger()
 {
    fDetector = NULL;

    KM3Dir = new G4UIdirectory("/KM3/");
   KM3Dir->SetGuidance("UI commands specific to this example.");

    XPDir = new G4UIdirectory("/KM3/XP/");
    XPDir->SetGuidance("detector control.");

    XPCommand = new G4UIcmdWithAString("/KM3/XP/Type", this);
    XPCommand->SetGuidance("Select The type of XP, AA K40 or WProp.\nAA can contain Antares(default), NEMO or KM3NeT");
    XPCommand->SetParameterName("choice", false);
    XPCommand->AvailableForStates(G4State_PreInit);

  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  KM3ExperimentTypeMessenger::~KM3ExperimentTypeMessenger() {
    delete XPCommand;
    KM3AnalysisManager::getIt ()->Write ();
    if (frunManager != NULL)
      delete frunManager;
    fInstance = NULL;
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  bool KM3ExperimentTypeMessenger::SetAnAnalysisMgr   (std::string sname, tAnMgrGet AMgr) {
    cout << "bon la c'est bon anamgr" << sname << endl;
    fAnalysisMgrCol[sname] = AMgr;
    return true;
  }
  bool KM3ExperimentTypeMessenger::SetASteppingAction (std::string sname, tStepActGet AMgr) {
    fSteppingActionCol[sname] = AMgr;
    cout << "bon la c'est bon stepping action" << sname << endl;
    return true;
  }

  void KM3ExperimentTypeMessenger::SetNewValue(G4UIcommand* command,
                                               G4String newValue) {
    if (command == XPCommand) {
      cout << newValue << " <- detector name " << fDetector <<endl;

      if (fDetector == NULL)
        {
          fDetector = KM3DetectorConstruction::GetIt ();
          frunManager = new KM3RunManager;

          // set mandatory initialization classes

          frunManager->SetUserInitialization(new KM3PhysicsList);

          // set mandatory user action class

          frunManager->SetUserAction(new KM3RunAction);
          frunManager->SetUserAction(new KM3EventAction);

          //first checking the new method, then for retro-compality the old ones
          if (fAnalysisMgrCol.count (newValue))
            {
              fAnalysisMgrCol[newValue] ();
            }
          if (fSteppingActionCol.count(newValue)) {
            cout << "create the stepping "<< newValue << endl;
            fSteppingAction = fSteppingActionCol[newValue] ();
          }

          ////TODO this part should be removed once the concept is validated
          if (newValue == "WProp")
            {
              WPropAnalysisManager::getIt ();
              fSteppingAction=new WPropSteppingAction ();

            }
          else if (newValue.substr (0,4)  == "XPPh")
            {
              AAAnalysisManager::getIt ();
              fSteppingAction=new XPPhSteppingAction ();
            } else if (fSteppingAction == 0)
            fSteppingAction=new KM3SteppingAction ();

          if (newValue == "K40")
            {
              K40AnalysisManager::getIt ()->SetTargetVolume (pow(fDetector->GetTargetFullLength(),3)*M_PI*4./3);
            }
          else if (newValue.substr (0,2) == "AA")
            {
              AAAnalysisManager::getIt ();
            }
          ////TODO

          if (!fSteppingAction)
            fSteppingAction=new KM3SteppingAction ();
          if (!fTrackingAction)
            fTrackingAction=new KM3TrackingAction ();

          KM3AnalysisManager::getIt ()->Init ();
          KM3PrimaryGeneratorAction* PGA=new KM3PrimaryGeneratorAction;
          PGA->SetProducingVolumeRadius((fDetector)->GetTargetFullLength ());
          frunManager->SetUserAction(PGA);
          frunManager->SetUserInitialization(fDetector);
          frunManager->SetUserAction(fSteppingAction);
          frunManager->SetUserAction(fTrackingAction);


        } else {
        throw("XP is already initiated!! The command /KM3/XP/Type should be done at the very beginning and only once!!\n");
      }
    }
  }

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
