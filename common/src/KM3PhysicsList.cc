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

#include <G4HadronFissionProcess.hh>
#ifdef USE_FISSION_NEW
#include <G4FissLib_new.hh>
#endif
#include <G4IonConstructor.hh>
/////--------------------------------------
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessTable.hh"

// Processes

#include "G4HadronElasticProcess.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4NeutronHPThermalScatteringData.hh"
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPThermalScattering.hh"

#include "G4NeutronInelasticProcess.hh"
#include "G4NeutronHPInelasticData.hh"
#include "G4NeutronHPInelastic.hh"

#include "G4HadronCaptureProcess.hh"
#include "G4NeutronHPCaptureData.hh"
#include "G4NeutronHPCapture.hh"

#include "G4HadronFissionProcess.hh"
#include "G4NeutronHPFissionData.hh"
#include "G4NeutronHPFission.hh"

/////--------------------------------------


#include "KM3PhysicsList.hh"

#include "globals.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4ProcessManager.hh"

#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "G4OpAbsorption.hh"
#include <G4OpWLS.hh>
#include "G4OpRayleigh.hh"
#include "G4OpMieHG.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4Decay.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4DecayPhysics.hh"

#include "G4LossTableManager.hh"
#include "G4EmSaturation.hh"

#include "G4RegionStore.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include <G4ionIonisation.hh>
#include <G4hMultipleScattering.hh>
#include <G4PreCompoundModel.hh>
#include <G4GeneratorPrecompoundInterface.hh>
#include <G4CrossSectionDataSetRegistry.hh>
#include <G4TheoFSGenerator.hh>
#include <G4CascadeInterface.hh>
#include <G4AlphaInelasticProcess.hh>
#include <G4GGNuclNuclCrossSection.hh>
#include <G4FTFModel.hh>
#include <G4ExcitedStringDecay.hh>
#include <G4LundStringFragmentation.hh>
#include <G4NuclearStopping.hh>

#include "G4hIonisation.hh"
#include "G4FastSimulationManagerProcess.hh"

#include <KM3PetzoldScattering.hh>
#include <KM3PMTOpticalModel.hh>
#include <KM3PhysicsListMessenger.hh>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <iostream>
using namespace std;


namespace km3net
{

  KM3PhysicsList::KM3PhysicsList() {
    //G4VModularPhysicsList() {
    //G4LossTableManager::Instance();
    defaultCutValue = 0.1 * mm;
    cutForGamma = defaultCutValue;
    cutForElectron = defaultCutValue;
    cutForPositron = defaultCutValue;

    DetectorCuts = 0;
    TargetCuts = 0;

    //SetVerboseLevel(1);

    theCerenkovProcess = NULL;
    theScintillationProcess = NULL;
    theAbsorptionProcess = NULL;
    theRayleighScatteringProcess = NULL;
    theBoundaryProcess = NULL;
    //default physics
    fParticleList = new G4DecayPhysics();
    //default physics
    fRaddecayList = new G4RadioactiveDecayPhysics();
    //km3 petzold process
    thePetzoldScattering = NULL;
    pMessenger=new KM3PhysicsListMessenger (this);

    theBoundaryProcess = new KM3PMTOpticalModel ();//create here for the messenger of opticalmodel
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  KM3PhysicsList::~KM3PhysicsList() {
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void KM3PhysicsList::ConstructParticle() {
    ConstructBosons();
    ConstructLeptons();
    ConstructMesons();
    ConstructBaryons();
    ConstructIons();
    G4IonConstructor pIonConstructor;
    pIonConstructor.ConstructParticle();
    fParticleList->ConstructParticle();
  }

  void KM3PhysicsList::ConstructIons() {
    G4Alpha::Definition ();
    //other can be H2O, H2...
  }

  void KM3PhysicsList::ConstructBosons() {
    // gamma
    G4Gamma::GammaDefinition();
    // optical photon
    G4OpticalPhoton::OpticalPhotonDefinition();
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void KM3PhysicsList::ConstructLeptons() {
    // leptons
    //  e+/-
    G4Electron::ElectronDefinition();
    G4Positron::PositronDefinition();
    // mu+/-
    G4MuonPlus::MuonPlusDefinition();
    G4MuonMinus::MuonMinusDefinition();
    // nu_e
    G4NeutrinoE::NeutrinoEDefinition();
    G4AntiNeutrinoE::AntiNeutrinoEDefinition();
    // nu_mu
    G4NeutrinoMu::NeutrinoMuDefinition();
    G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void KM3PhysicsList::ConstructMesons() {
    //  mesons
    G4PionPlus::PionPlusDefinition();
    G4PionMinus::PionMinusDefinition();
    G4PionZero::PionZeroDefinition();
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void KM3PhysicsList::ConstructBaryons() {
    //  barions
    G4Proton::ProtonDefinition();
    G4AntiProton::AntiProtonDefinition();

    G4Neutron::NeutronDefinition();
    G4AntiNeutron::AntiNeutronDefinition();
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void KM3PhysicsList::ConstructProcess() {
    AddTransportation();
    //AddParameterisation();
    ConstructGeneral();
    ConstructEM();
    ConstructOp();
    ConstructNeutronPhysics ();


  }

  void KM3PhysicsList::ConstructGeneral() {
    /*// Add Decay Process
      G4Decay* theDecayProcess = new G4Decay();
      theParticleIterator->reset();
      while ((*theParticleIterator)()) {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      if (theDecayProcess->IsApplicable(*particle)) {
      pmanager->AddProcess(theDecayProcess);
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager->SetProcessOrdering(theDecayProcess, idxAtRest);
      }
      }*/
    fParticleList->ConstructProcess();
    fRaddecayList->ConstructProcess();

  }

  double rng4llnlfisslib(void)
  {
    return G4UniformRand();
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  void KM3PhysicsList::ConstructNeutronPhysics () {
    G4ProcessManager* pManager = G4Neutron::Neutron()->GetProcessManager();

    // process: elastic
    //
    G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess();
    pManager->AddDiscreteProcess(theElasticProcess);
    //
    // cross section data set
    G4NeutronHPElasticData* dataSet1a = new G4NeutronHPElasticData();
    G4NeutronHPThermalScatteringData* dataSet1b
      = new G4NeutronHPThermalScatteringData();
    theElasticProcess->AddDataSet(dataSet1a);
    theElasticProcess->AddDataSet(dataSet1b);
    //
    // models
    G4NeutronHPElastic*           model1a = new G4NeutronHPElastic();
    G4NeutronHPThermalScattering* model1b = new G4NeutronHPThermalScattering();
    model1a->SetMinEnergy(4*eV);
    theElasticProcess->RegisterMe(model1a);
    theElasticProcess->RegisterMe(model1b);

    // process: inelastic
    //
    G4NeutronInelasticProcess* theInelasticProcess = new G4NeutronInelasticProcess();
    pManager->AddDiscreteProcess(theInelasticProcess);
    //
    // cross section data set
    G4NeutronHPInelasticData* dataSet2 = new G4NeutronHPInelasticData();
    theInelasticProcess->AddDataSet(dataSet2);
    //
    // models
    G4NeutronHPInelastic* model2 = new G4NeutronHPInelastic();
    theInelasticProcess->RegisterMe(model2);

    // process: nCapture
    //
    G4HadronCaptureProcess* theCaptureProcess = new G4HadronCaptureProcess();
    pManager->AddDiscreteProcess(theCaptureProcess);
    //
    // cross section data set
    G4NeutronHPCaptureData* dataSet3 = new G4NeutronHPCaptureData();
    theCaptureProcess->AddDataSet(dataSet3);
    //
    // models
    G4NeutronHPCapture* model3 = new G4NeutronHPCapture();
    theCaptureProcess->RegisterMe(model3);

    // process: nFission
    //
    G4HadronFissionProcess* theFissionProcess = new G4HadronFissionProcess();
    pManager->AddDiscreteProcess(theFissionProcess);
    //
    // cross section data set
    G4NeutronHPFissionData* dataSet4 = new G4NeutronHPFissionData();
    theFissionProcess->AddDataSet(dataSet4);
    //
    // models
#ifdef USE_FISSION_NEW
#ifdef FISSION_NEW
    // pass the random number generator to class fissionEvent
    fissionEvent::setRNGd(rng4llnlfisslib);
#ifdef USEFREYA
    G4cout << "NeutronHPphysics: using new version of LLNL Fission Library (overriding built-in version) with FREYA turned on\n";
    fissionEvent::setCorrelationOption(3);
#else
    G4cout << "NeutronHPphysics: using new version of LLNL Fission Library (overriding built-in version) without FREYA\n";
#endif
    G4FissLib_new* theFissionModel = new G4FissLib_new;
#else
    G4cout <<"NeutronHPphysics: using built-in version of LLNL Fission Library\n";
    G4FissLib* theFissionModel = new G4FissLib;
#endif
    theFissionProcess->RegisterMe(theFissionModel);
#endif

  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void KM3PhysicsList::ConstructEM() {
    theParticleIterator->reset();
    while ((*theParticleIterator)()) {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();

      if (particleName == "gamma") {
        // gamma
        // Construct processes for gamma
        pmanager->AddDiscreteProcess(new G4GammaConversion());
        pmanager->AddDiscreteProcess(new G4ComptonScattering());
        pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());

      } else if (particleName == "e-") {
        //electron
        // Construct processes for electron
        pmanager->AddProcess(new G4eMultipleScattering(), -1, 1, 1);
        pmanager->AddProcess(new G4eIonisation(), -1, 2, 2);
        pmanager->AddProcess(new G4eBremsstrahlung(), -1, 3, 3);

      } else if (particleName == "e+") {
        //positron
        // Construct processes for positron
        pmanager->AddProcess(new G4eMultipleScattering(), -1, 1, 1);
        pmanager->AddProcess(new G4eIonisation(), -1, 2, 2);
        pmanager->AddProcess(new G4eBremsstrahlung(), -1, 3, 3);
        pmanager->AddProcess(new G4eplusAnnihilation(), 0, -1, 4);

      } else if (particleName == "mu+" || particleName == "mu-") {
        //muon
        // Construct processes for muon
        pmanager->AddProcess(new G4MuMultipleScattering(), -1, 1, 1);
        pmanager->AddProcess(new G4MuIonisation(), -1, 2, 2);
        pmanager->AddProcess(new G4MuBremsstrahlung(), -1, 3, 3);
        pmanager->AddProcess(new G4MuPairProduction(), -1, 4, 4);

      }
      else if (particleName == "alpha" || particleName == "He3") {
        pmanager->AddProcess(new G4hMultipleScattering(), ordInActive,1,1);
        G4ionIonisation* ionIoni = new G4ionIonisation();
        ionIoni->SetStepFunction(0.1, 20*um);
        pmanager->AddProcess(ionIoni,                   -1, 2, 2);


        const G4double theFTFMin0 =    100*GeV;
        //const G4double theFTFMin1 =    4.0*GeV;
        const G4double theFTFMax =   100.0*TeV;
        const G4double theBERTMin0 =   0.0*GeV;
        //const G4double theBERTMin1 =  19.0*MeV;
        const G4double theBERTMax =    5.0*GeV;

        G4FTFModel * theStringModel = new G4FTFModel;
        G4ExcitedStringDecay * theStringDecay = new G4ExcitedStringDecay( new G4LundStringFragmentation );
        theStringModel->SetFragmentationModel( theStringDecay );


        G4PreCompoundModel * thePreEquilib = new G4PreCompoundModel( new G4ExcitationHandler );
        G4GeneratorPrecompoundInterface * theCascade = new G4GeneratorPrecompoundInterface( thePreEquilib );

        G4VCrossSectionDataSet * theGGNuclNuclData = G4CrossSectionDataSetRegistry::Instance()->  GetCrossSectionDataSet(G4GGNuclNuclCrossSection::Default_Name());
        G4TheoFSGenerator * theFTFModel1 = new G4TheoFSGenerator( "FTFP" );
        theFTFModel1->SetHighEnergyGenerator( theStringModel );
        theFTFModel1->SetTransport( theCascade );
        theFTFModel1->SetMinEnergy( theFTFMin0 );
        theFTFModel1->SetMaxEnergy( theFTFMax );

        G4CascadeInterface * theBERTModel0 = new G4CascadeInterface;
        theBERTModel0->SetMinEnergy( theBERTMin0 );
        theBERTModel0->SetMaxEnergy( theBERTMax );


        G4AlphaInelasticProcess* theInelasticProcess =
          new G4AlphaInelasticProcess("inelastic");
        theInelasticProcess->AddDataSet( theGGNuclNuclData );
        theInelasticProcess->RegisterMe( theFTFModel1 );
        theInelasticProcess->RegisterMe( theBERTModel0 );
        pmanager->AddDiscreteProcess( theInelasticProcess );
        pmanager->AddProcess (new G4NuclearStopping, -1, 3, -1);
      }
      else {
        if ((particle->GetPDGCharge() != 0.0)
            && (particle->GetParticleName() != "chargedgeantino")) {
          // all others charged particles except geantino
          pmanager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
          G4ionIonisation* ionIoni = new G4ionIonisation();
          ionIoni->SetStepFunction(0.1, 20*um);
          pmanager->AddProcess(ionIoni,                   -1, 2, 2);



        }
      }
    }

  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void KM3PhysicsList::ConstructOp() {
    theCerenkovProcess = new G4Cerenkov("Cerenkov");
    theAbsorptionProcess = new G4OpAbsorption();
    theRayleighScatteringProcess = new G4OpRayleigh();
    //theBoundaryProcess = new KM3PMTOpticalModel (); done in the constructor for messenger
    theWLSProcess = new G4OpWLS("OpWLS");

    thePetzoldScattering = new KM3PetzoldScattering ();
    thePetzoldScattering->SetKopelFactor (fKopelFactor);
    thePetzoldScattering->SetESFactor (fESFactor);

    theCerenkovProcess->SetMaxNumPhotonsPerStep(20);
    theCerenkovProcess->SetMaxBetaChangePerStep(10.0);
    theCerenkovProcess->SetTrackSecondariesFirst(true);

    // default scintillation process
    G4Scintillation* theScintProcessDef = new G4Scintillation("Scintillation");
    // theScintProcessDef->DumpPhysicsTable();
    theScintProcessDef->SetTrackSecondariesFirst(true);
    theScintProcessDef->SetScintillationYieldFactor(1.0); //
    theScintProcessDef->SetScintillationExcitationRatio(0.0); //

    // scintillation process for alpha:
    G4Scintillation* theScintProcessAlpha = new G4Scintillation("Scintillation");
    // theScintProcessNuc->DumpPhysicsTable();
    theScintProcessAlpha->SetTrackSecondariesFirst(true);
    theScintProcessAlpha->SetScintillationYieldFactor(1.1);
    theScintProcessAlpha->SetScintillationExcitationRatio(1.0);

    // scintillation process for heavy nuclei
    G4Scintillation* theScintProcessNuc = new G4Scintillation("Scintillation");
    // theScintProcessNuc->DumpPhysicsTable();
    theScintProcessNuc->SetTrackSecondariesFirst(true);
    theScintProcessNuc->SetScintillationYieldFactor(0.2);
    theScintProcessNuc->SetScintillationExcitationRatio(1.0);

    theParticleIterator->reset();
    while ((*theParticleIterator)()) {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
      if (theCerenkovProcess->IsApplicable(*particle)) {
        pmanager->AddProcess(theCerenkovProcess);
        pmanager->SetProcessOrdering(theCerenkovProcess, idxPostStep);
      }
      if (theScintProcessDef->IsApplicable(*particle)) {
        if(particle->GetParticleName() == "GenericIon") {
	  pmanager->AddProcess(theScintProcessNuc); // AtRestDiscrete
	  pmanager->SetProcessOrderingToLast(theScintProcessNuc,idxAtRest);
	  pmanager->SetProcessOrderingToLast(theScintProcessNuc,idxPostStep);
	}
	else if(particle->GetParticleName() == "alpha") {
	  pmanager->AddProcess(theScintProcessAlpha);
	  pmanager->SetProcessOrderingToLast(theScintProcessAlpha,idxAtRest);
	  pmanager->SetProcessOrderingToLast(theScintProcessAlpha,idxPostStep);
	}
	else {
	  pmanager->AddProcess(theScintProcessDef);
	  pmanager->SetProcessOrderingToLast(theScintProcessDef,idxAtRest);
	  pmanager->SetProcessOrderingToLast(theScintProcessDef,idxPostStep);
	}
      }
      if (particleName == "opticalphoton") {
        G4cout << " AddDiscreteProcess to OpticalPhoton " << G4endl;
        pmanager->AddDiscreteProcess(theAbsorptionProcess);
        pmanager->AddDiscreteProcess(theRayleighScatteringProcess);
        pmanager->AddDiscreteProcess(thePetzoldScattering);
        pmanager->AddDiscreteProcess(theBoundaryProcess);
        pmanager->AddDiscreteProcess(thePetzoldScattering);
        pmanager->AddDiscreteProcess(theWLSProcess);
      }
    }
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void KM3PhysicsList::SetVerbose(G4int verbose) {
    theCerenkovProcess->SetVerboseLevel(verbose);
    theScintillationProcess->SetVerboseLevel(verbose);
    theAbsorptionProcess->SetVerboseLevel(verbose);
    theRayleighScatteringProcess->SetVerboseLevel(verbose);
    theBoundaryProcess->SetVerboseLevel(verbose);
    theWLSProcess->SetVerboseLevel (verbose);
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void KM3PhysicsList::SetCuts() {

    SetCutValue(cutForGamma, "gamma");
    SetCutValue(cutForElectron, "e-");
    SetCutValue(cutForPositron, "e+");
    G4cout << "world cuts are set" << G4endl;

    /*if (!TargetCuts)
      SetTargetCut(cutForElectron);
      G4Region* region = (G4RegionStore::GetInstance())->GetRegion("Target");
      region->SetProductionCuts(TargetCuts);
      G4cout << "Target cuts are set" << G4endl;

      if (!DetectorCuts)
      SetDetectorCut(cutForElectron);
      region = (G4RegionStore::GetInstance())->GetRegion("Detector");
      region->SetProductionCuts(DetectorCuts);
      G4cout << "Detector cuts are set" << G4endl; changedfortest
    */
    if (verboseLevel > 0)
      DumpCutValuesTable();
  }
  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void KM3PhysicsList::SetCutForGamma(G4double cut) {
    cutForGamma = cut;
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void KM3PhysicsList::SetCutForElectron(G4double cut) {
    cutForElectron = cut;
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void KM3PhysicsList::SetCutForPositron(G4double cut) {
    cutForPositron = cut;
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void KM3PhysicsList::SetTargetCut(G4double cut) {
    if (!TargetCuts)
      TargetCuts = new G4ProductionCuts();

    TargetCuts->SetProductionCut(cut, idxG4GammaCut);
    TargetCuts->SetProductionCut(cut, idxG4ElectronCut);
    TargetCuts->SetProductionCut(cut, idxG4PositronCut);

  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void KM3PhysicsList::SetDetectorCut(G4double cut) {
    if (!DetectorCuts)
      DetectorCuts = new G4ProductionCuts();

    DetectorCuts->SetProductionCut(cut, idxG4GammaCut);
    DetectorCuts->SetProductionCut(cut, idxG4ElectronCut);
    DetectorCuts->SetProductionCut(cut, idxG4PositronCut);
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


  void KM3PhysicsList::AddParameterisation() {
    G4FastSimulationManagerProcess* theFastSimulationManagerProcess=
      new G4FastSimulationManagerProcess("OpticalPhotocathodeProcess");
    theParticleIterator->reset();
    while ((*theParticleIterator)()) {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      // both postStep and alongStep action are required if the detector
      // makes use of ghost volumes. If no ghost, the postStep
      // is sufficient (and faster?).
      //pmanager->AddDiscreteProcess(fastSimProcess_massGeom);
      pmanager->AddProcess(theFastSimulationManagerProcess, -1, 0, 0);
    }

  }
}
