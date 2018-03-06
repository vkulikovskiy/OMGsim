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

/*
  Might be good to limit for normal mode
*/

#ifndef KM3PhysicsList_h
#define KM3PhysicsList_h 1

//#include "G4VModularPhysicsList.hh"
#include <G4VUserPhysicsList.hh>
#include "globals.hh"
#include <vector>
#include <G4SystemOfUnits.hh>

class G4VPhysicsConstructor;
class G4ProductionCuts;
class G4Cerenkov;
class G4Scintillation;
class G4OpAbsorption;
//G4OpBoundaryProcess;
class G4OpWLS;
class G4OpRayleigh;
class G4RadioactiveDecayPhysics;

namespace km3net
{
  class KM3PhysicsListMessenger;
  class KM3ParametrisedScattering;
  class KM3PetzoldScattering;
  class KM3PMTOpticalModel;
  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  class KM3PhysicsList: public G4VUserPhysicsList {//G4VModularPhysicsList {
  public:
    friend class KM3PhysicsListMessenger;
    friend class KM3PMTOpticalModel;

    KM3PhysicsList();
    virtual ~KM3PhysicsList();

    void ConstructParticle();

    void SetCuts();
    void SetCutForGamma(G4double);
    void SetCutForElectron(G4double);
    void SetCutForPositron(G4double);

    void SelectPhysicsList(const G4String& name);
    void ConstructProcess();

    void SetTargetCut(G4double val);
    void SetDetectorCut(G4double val);

    void AddParameterisation();

    KM3PMTOpticalModel* GetOpticalModel ()
    {
      return theBoundaryProcess;
    }

  private:

    void AddExtraBuilders(G4bool flagHP);
    void ConstructBosons();
    void ConstructLeptons();
    void ConstructMesons();
    void ConstructBaryons();
    void ConstructIons();
    void ConstructGeneral();
    void ConstructEM();
    void ConstructOp();
    void ConstructNeutronPhysics ();
    void SetVerbose(G4int);

    // hide assignment operator
    KM3PhysicsList & operator=(const KM3PhysicsList &right);
    KM3PhysicsList(const KM3PhysicsList&);

    G4double cutForGamma;
    G4double cutForElectron;
    G4double cutForPositron;

    G4VPhysicsConstructor* emPhysicsList;
    G4VPhysicsConstructor* fRaddecayList;
    G4VPhysicsConstructor* fParticleList;
    G4VPhysicsConstructor* hadPhysicsList;

    std::vector<G4VPhysicsConstructor*> hadronPhys;
    G4int nhadcomp;

    KM3PhysicsListMessenger* pMessenger;
    G4ProductionCuts* DetectorCuts;
    G4ProductionCuts* TargetCuts;

    G4Cerenkov* theCerenkovProcess;
    G4Scintillation* theScintillationProcess;
    G4OpAbsorption* theAbsorptionProcess;
    G4OpRayleigh* theRayleighScatteringProcess;
    KM3PMTOpticalModel* theBoundaryProcess;
    //G4OpBoundaryProcess* theBoundaryProcess;
    G4OpWLS* theWLSProcess;

    KM3PetzoldScattering* thePetzoldScattering;

    double fKopelFactor = 1;
    double fESFactor = 1;

  };

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
}
#endif
