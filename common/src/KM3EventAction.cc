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
// CHANGE HISTORY
// --------------
#include "G4ios.hh"
#include "KM3EventAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"

#include <KM3AnalysisManager.hh>

#include <G4SystemOfUnits.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include <iostream>
using namespace std;

namespace km3net
{
  using namespace CLHEP;

  KM3EventAction::KM3EventAction()
  : drawFlag("all")
  {
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

  KM3EventAction::~KM3EventAction()
  {

  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

  void KM3EventAction::BeginOfEventAction(const G4Event* Ev)
  {
    KM3AnalysisManager::getIt()->SetOriginOfEvent (Ev);//->GetPrimaryVertex ()->GetPosition ());
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

  void KM3EventAction::EndOfEventAction(const G4Event* evt)
  {
    G4int event_id       = evt->GetEventID();
    //analysis
    KM3AnalysisManager::getIt()->EndOfEvent ();
    // visualisation
    if (event_id < 1000000 && G4VVisManager::GetConcreteInstance()) {
      G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
      G4int n_trajectories = 0;
      if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
      for (G4int i=0; i<n_trajectories; i++) {
        G4Trajectory* trj = (G4Trajectory *)
          ((*(evt->GetTrajectoryContainer()))[i]);
        if (drawFlag == "all") {
          trj->DrawTrajectory();
        } else if (drawFlag == "charged" && trj->GetCharge() != 0.) {
          trj->DrawTrajectory();
        }
      }
    }
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
