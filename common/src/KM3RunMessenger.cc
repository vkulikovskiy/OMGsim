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

#include "KM3RunMessenger.hh"

#include "KM3RunManager.hh"
#include <KM3PrimaryGeneratorAction.hh>
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace km3net
{

  KM3RunMessenger::KM3RunMessenger(KM3RunManager* RUN) :
  myRun(RUN) {
    KM3Dir = new G4UIdirectory("/KM3/");
    KM3Dir->SetGuidance("UI commands specific to this program.");

    fRunDir = new G4UIdirectory("/KM3/GenType/");
    fRunDir->SetGuidance("GenType control 1 normal (AA), 2 radioactivity (K40).");

    TargLengthCmd = new G4UIcmdWithADoubleAndUnit("/KM3/det/setTargetLength",
                                                  this);
    TargLengthCmd->SetGuidance("Set the Target Length.");
    TargLengthCmd->SetUnitCategory("Length");
    TargLengthCmd->SetParameterName("choice", false);
    TargLengthCmd->AvailableForStates(G4State_PreInit);

    fGenType = new G4UIcmdWithAnInteger
      ("/KM3/GenType", this);
    fGenType->SetGuidance("1 for normal emission, 2 for ion emission.");
    fGenType->SetParameterName("choice", false);
    fGenType->AvailableForStates(G4State_PreInit);
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  KM3RunMessenger::~KM3RunMessenger() {
    delete TargLengthCmd;
    delete fGenType;
    delete fRunDir;
    delete KM3Dir;
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


  void KM3RunMessenger::SetNewValue(G4UIcommand* command,
                                         G4String newValue) {
    KM3PrimaryGeneratorAction* the_UPGA=(KM3PrimaryGeneratorAction*)((void *)(myRun->GetUserPrimaryGeneratorAction ()));

    if (command == TargLengthCmd) {
      the_UPGA->fTargetRadius=TargLengthCmd->GetNewDoubleValue(newValue);
    } else if (command == fGenType) {
      int thetype=fGenType->GetNewIntValue(newValue);
      if (thetype == 1)
        the_UPGA->InitializeNormalType ();
      else if (thetype == 2)
        the_UPGA->InitializeDecayType ();
    }
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
