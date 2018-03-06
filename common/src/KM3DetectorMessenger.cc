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


#include <KM3DetectorMessenger.hh>
#include <KM3Material.hh>
#include "KM3DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "globals.hh"

#include <iostream>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace km3net
{

  KM3DetectorMessenger::KM3DetectorMessenger(KM3DetectorConstruction* DET) :
    myDetector(DET) {
    KM3Dir = new G4UIdirectory("/KM3/");
    KM3Dir->SetGuidance("UI commands specific to this example.");

    detDir = new G4UIdirectory("/KM3/det/");
    detDir->SetGuidance("detector control.");

    TargMatCmd = new G4UIcmdWithAString("/KM3/det/setTargetMaterial", this);
    TargMatCmd->SetGuidance("Select Material of the Target.");
    TargMatCmd->SetParameterName("choice", false);
    TargMatCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    TargScatCmd = new G4UIcmdWithABool("/KM3/det/setTargetScattering", this);
    TargScatCmd->SetGuidance("Select if the Material of the Target has scattering (default true).");
    TargScatCmd->SetParameterName("choice", false);
    TargScatCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    TargLengthCmd = new G4UIcmdWithADoubleAndUnit("/KM3/det/setTargetLength",
                                                  this);
    TargLengthCmd->SetGuidance("Set the Target Length.");
    TargLengthCmd->SetUnitCategory("Length");
    TargLengthCmd->SetParameterName("choice", false);
    TargLengthCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

    WorldLengthCmd = new G4UIcmdWithADoubleAndUnit("/KM3/det/setWorldLength",
                                                   this);
    WorldLengthCmd->SetGuidance("Set the World Length.");
    WorldLengthCmd->SetUnitCategory("Length");
    WorldLengthCmd->SetParameterName("choice", false);
    WorldLengthCmd->AvailableForStates(G4State_PreInit);

    GenerationVolumeLengthCmd = new G4UIcmdWithADoubleAndUnit(
                                                              "/KM3/det/setGenerationVolumeLength", this);
    GenerationVolumeLengthCmd->SetGuidance("Set the Generation Volume Length.");
    GenerationVolumeLengthCmd->SetUnitCategory("Length");
    GenerationVolumeLengthCmd->SetParameterName("choice", false);
    GenerationVolumeLengthCmd->AvailableForStates(G4State_PreInit);

    GenerationVolumeRadiusCmd = new G4UIcmdWithADoubleAndUnit(
                                                              "/KM3/det/setGenerationVolumeRadius", this);
    GenerationVolumeRadiusCmd->SetGuidance("Set the Generation Volume Radius.");
    GenerationVolumeRadiusCmd->SetUnitCategory("Length");
    GenerationVolumeRadiusCmd->SetParameterName("choice", false);
    GenerationVolumeRadiusCmd->AvailableForStates(G4State_PreInit);

    StoreyRadiusCmd = new G4UIcmdWithADoubleAndUnit("/KM3/det/setStoreyRadius",
                                                    this);
    StoreyRadiusCmd->SetGuidance("Set the Storey Radius.");
    StoreyRadiusCmd->SetUnitCategory("Length");
    StoreyRadiusCmd->SetParameterName("choice", false);
    StoreyRadiusCmd->AvailableForStates(G4State_PreInit);

    OpticalModulePhi = new G4UIcmdWithADoubleAndUnit(
                                                     "/KM3/det/setOpticalModulePhi", this);
    OpticalModulePhi->SetGuidance("Set the Optical Module Orientation.");
    OpticalModulePhi->SetUnitCategory("Angle");
    OpticalModulePhi->SetParameterName("choice", false);
    OpticalModulePhi->AvailableForStates(G4State_PreInit);

    OpticalModuleTheta = new G4UIcmdWithADoubleAndUnit(
                                                       "/KM3/det/setOpticalModuleTheta", this);
    OpticalModuleTheta->SetGuidance("Set the Optical Module Orientation.");
    OpticalModuleTheta->SetUnitCategory("Angle");
    OpticalModuleTheta->SetParameterName("choice", false);
    OpticalModuleTheta->AvailableForStates(G4State_PreInit);

    OpticalModulePsi = new G4UIcmdWithADoubleAndUnit(
                                                     "/KM3/det/setOpticalModulePsi", this);
    OpticalModulePsi->SetGuidance("Set the Optical Module Orientation.");
    OpticalModulePsi->SetUnitCategory("Angle");
    OpticalModulePsi->SetParameterName("choice", false);
    OpticalModulePsi->AvailableForStates(G4State_PreInit);

    MeanCapRadius = new G4UIcmdWithADoubleAndUnit("/KM3/det/setMeanCapRadius", this);
    MeanCapRadius->SetGuidance("for nemo struct. 0 is without cap (antares).");
    MeanCapRadius->SetUnitCategory("Length");
    MeanCapRadius->SetParameterName("choice", false);
    MeanCapRadius->AvailableForStates(G4State_PreInit);

    DontDraw = new G4UIcmdWithAString
      ("/KM3/det/DontDraw", this);
    DontDraw->SetGuidance("To avoid to draw some detector elements. Parameters are bento, gel, absorber. Whatever can be a separator (even without it is valid).");
    DontDraw->SetParameterName("choice", false);
    DontDraw->AvailableForStates(G4State_PreInit);

    DetectorTypeName = new G4UIcmdWithAString
      ("/KM3/det/DetectorTypeName", this);
    DetectorTypeName->SetGuidance("new way to create OM dispositions. The name should correspond of one from common/data/KM3Det*.dat.");
    DetectorTypeName->SetParameterName("choice", false);
    DetectorTypeName->AvailableForStates(G4State_PreInit);

    DataFilesSources = new G4UIcmdWithAString
      ("/KM3/det/DataFilesSources", this);
    DataFilesSources->SetGuidance("additionnal folders for data files (as in common/data/KM3Det*.dat). \":\" is the list separator.");
    DataFilesSources->SetParameterName("choice", false);
    DataFilesSources->AvailableForStates(G4State_PreInit);

    ParametersList = new G4UIcmdWithAString
      ("/KM3/det/ParametersList", this);
    ParametersList->SetGuidance("additive parameters for construction, should be in the format \"par1__ = 5;par2__ = 6...");
    ParametersList->SetParameterName("choice", false);
    ParametersList->AvailableForStates(G4State_PreInit);


  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  KM3DetectorMessenger::~KM3DetectorMessenger() {
    delete TargMatCmd;
    delete TargScatCmd;
    delete TargLengthCmd;
    delete GenerationVolumeLengthCmd;
    delete GenerationVolumeRadiusCmd;
    delete StoreyRadiusCmd;
    delete OpticalModuleTheta;
    delete OpticalModulePhi;
    delete OpticalModulePsi;
    delete MeanCapRadius;
    delete detDir;
    delete KM3Dir;
    delete DetectorTypeName;
    delete DataFilesSources;
    delete ParametersList;
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


  void KM3DetectorMessenger::SetNewValue(G4UIcommand* command,
                                         G4String newValue) {
    if (command == TargMatCmd) {
      myDetector->SetTargetMaterial(newValue);
    } else if (command == TargScatCmd) {
      if (!TargScatCmd->GetNewBoolValue(newValue))
        KM3Material::GetIt ()->RemoveScattering ();
    } else if (command == TargLengthCmd) {
      myDetector->SetTargetLength(TargLengthCmd->GetNewDoubleValue(newValue));
    } else if (command == GenerationVolumeLengthCmd) {
      myDetector->SetGenerationVolumeLength
        (GenerationVolumeLengthCmd->GetNewDoubleValue(newValue));
    } else if (command == GenerationVolumeRadiusCmd) {
      myDetector->SetGenerationVolumeRadius
        (GenerationVolumeRadiusCmd->GetNewDoubleValue(newValue));
    } else if (command == WorldLengthCmd) {
      myDetector->SetTargetLength(WorldLengthCmd->GetNewDoubleValue(newValue));
    } else if (command == OpticalModuleTheta) {
      myDetector->SetOrientationTheta
        (OpticalModuleTheta->GetNewDoubleValue(newValue));
    } else if (command == OpticalModulePhi) {
      myDetector->SetOrientationPhi
        (OpticalModulePhi->GetNewDoubleValue(newValue));
    } else if (command == OpticalModulePsi) {
      myDetector->SetOrientationPsi
        (OpticalModulePsi->GetNewDoubleValue(newValue));
    } else if (command == MeanCapRadius) {
      myDetector->SetMeanCapRadius
        (MeanCapRadius->GetNewDoubleValue(newValue));
    } else if (command == DontDraw) {
      myDetector->DontDraw
        (newValue);
    } else if (command == DetectorTypeName) {
      myDetector->SetDetectorTypeName (newValue);
    } else if (command == DataFilesSources) {
      myDetector->SetDataFilesSources (newValue);
    } else if (command == ParametersList) {
      myDetector->SetParametersListFromMac (newValue);
    }


  }

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
