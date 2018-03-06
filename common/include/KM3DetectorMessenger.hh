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

#ifndef KM3DetectorMessenger_h
#define KM3DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"


class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace km3net
{

  class KM3DetectorConstruction;

  class KM3DetectorMessenger: public G4UImessenger
  {
  public:
    KM3DetectorMessenger(KM3DetectorConstruction*);
    ~KM3DetectorMessenger();

    void SetNewValue(G4UIcommand*, G4String);

  private:
    KM3DetectorConstruction* myDetector;

    G4UIdirectory*             KM3Dir;
    G4UIdirectory*             detDir;
    G4UIcmdWithAString*        TargMatCmd;
    G4UIcmdWithABool*          TargScatCmd;
    G4UIcmdWithADoubleAndUnit* TargLengthCmd;
    G4UIcmdWithADoubleAndUnit* GenerationVolumeLengthCmd;
    G4UIcmdWithADoubleAndUnit* GenerationVolumeRadiusCmd;
    G4UIcmdWithADoubleAndUnit* StoreyRadiusCmd;
    G4UIcmdWithADoubleAndUnit* OpticalModuleTheta;
    G4UIcmdWithADoubleAndUnit* OpticalModulePhi;
    G4UIcmdWithADoubleAndUnit* OpticalModulePsi;
    G4UIcmdWithADoubleAndUnit* MeanCapRadius;
    G4UIcmdWithAString*        DontDraw;
    G4UIcmdWithAString*        DetectorTypeName;
    G4UIcmdWithAString*        DataFilesSources;
    G4UIcmdWithAString*        ParametersList;
    G4UIcmdWithADoubleAndUnit* WorldLengthCmd;

  };
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
