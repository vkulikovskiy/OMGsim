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

#ifndef KM3ExperimentTypeMessenger_h
#define KM3ExperimentTypeMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

#include <map>
#include <functional>

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4RunManager;
class G4SteppingAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace km3net
{

  class KM3DetectorConstruction;
  class KM3AnalysisManager;
  class KM3SteppingAction;
  class KM3TrackingAction;


  class KM3ExperimentTypeMessenger: public G4UImessenger
  {
  private:
    typedef std::function<KM3AnalysisManager*()> tAnMgrGet;
    typedef std::function<KM3SteppingAction*()> tStepActGet;
  public:
    static KM3ExperimentTypeMessenger* GetIt ();
    bool SetAnAnalysisMgr   (std::string sname, tAnMgrGet AMgr);
    bool SetASteppingAction (std::string sname, tStepActGet AnAction);

  private:
    KM3ExperimentTypeMessenger( );
  public:
    ~KM3ExperimentTypeMessenger();

    void SetNewValue(G4UIcommand*, G4String);

  private:
    KM3DetectorConstruction*  fDetector;
    KM3SteppingAction*        fSteppingAction = NULL;
    KM3TrackingAction*        fTrackingAction =NULL;
    G4RunManager*             frunManager = NULL;

    G4UIdirectory*            KM3Dir;
    G4UIdirectory*            XPDir;
    G4UIcmdWithAString*       XPCommand;

    std::map <std::string, tStepActGet> fSteppingActionCol;
    std::map <std::string, tAnMgrGet > fAnalysisMgrCol;

    static KM3ExperimentTypeMessenger* fInstance;

  };
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

