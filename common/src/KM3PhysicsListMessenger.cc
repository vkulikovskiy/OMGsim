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

#include "KM3PhysicsListMessenger.hh"

#include "KM3PhysicsList.hh"
#include <KM3PetzoldScattering.hh>
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADouble.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <iostream>
using namespace std;

namespace km3net
{
  KM3PhysicsListMessenger::KM3PhysicsListMessenger(KM3PhysicsList* physicslist) :
    myPhysicsList(physicslist) {

    KM3Dir = new G4UIdirectory("/KM3/");
    KM3Dir->SetGuidance("UI commands specific to this program.");

    fPhysicsListDir = new G4UIdirectory("/KM3/PhysicsList/");
    fPhysicsListDir->SetGuidance("Factor to modify the water properties");

    fKopelLengthFactorCmd = new G4UIcmdWithADouble ("/KM3/PhysicsList/KopelLengthFactor", this);
    fKopelLengthFactorCmd->SetGuidance("Set the Kopelevitch scattering Length factor.");
    fKopelLengthFactorCmd->SetParameterName("choice", false);
    fKopelLengthFactorCmd->AvailableForStates(G4State_PreInit);

    fESLengthFactorCmd = new G4UIcmdWithADouble ("/KM3/PhysicsList/ESLengthFactor", this);
    fESLengthFactorCmd->SetGuidance("Set the ES scattering Length factor.");
    fESLengthFactorCmd->SetParameterName("choice", false);
    fESLengthFactorCmd->AvailableForStates(G4State_PreInit);
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  KM3PhysicsListMessenger::~KM3PhysicsListMessenger() {
    delete fKopelLengthFactorCmd;
    delete fESLengthFactorCmd;
    delete fPhysicsListDir;
    delete KM3Dir;
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


  void KM3PhysicsListMessenger::SetNewValue(G4UIcommand* command,
                                         G4String newValue) {
     if (command == fKopelLengthFactorCmd) {
       myPhysicsList->fKopelFactor = fKopelLengthFactorCmd->GetNewDoubleValue(newValue);
     }
    else if (command == fESLengthFactorCmd) {
      myPhysicsList->fESFactor=fESLengthFactorCmd->GetNewDoubleValue(newValue);
    }

  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
