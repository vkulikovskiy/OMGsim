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
#include "KM3RunManager.hh"
#include "G4UImanager.hh"
#ifdef G4UI_USE_GAG
#include "G4UIGAG.hh"
#endif
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif


#include "K40DetectorConstruction.hh"
#include "KM3PhysicsList.hh"
#include "KM3EventAction.hh"
#include "KM3RunAction.hh"
#include "KM3PrimaryGeneratorAction.hh"
#include "K40AnalysisManager.hh"
#include "KM3SteppingAction.hh"
#include "Randomize.hh"
#include <sys/time.h>
#include <unistd.h>

#include <random>
#include <iostream>

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

using namespace std;
using namespace km3net;

int main(int argc_,char** argv_)
{
  G4String srf="default.root";
  string MacFileName;




  bool debug = false;
  bool bsession = false;
  bool keepall = false;
  bool analyse_all = false;
  int iarg = 1;
  while (iarg < argc_)
    {
      string token = argv_[iarg];

      if (token[0] == '-')
        {
          string option = token;
          if ((option == "-d") || (option == "--debug"))
            {
              debug = true;
            }
          else
            if ((option == "-s") || (option == "--session"))
              {
                bsession = true;
              }
            else
              if ((option == "-ka") || (option == "--keep-all"))
                {
                  keepall = true;
                }
              else
                if ((option == "-fa") || (option == "--full-analysis"))
                  {
                    analyse_all=true;
                  }

          /* Here you may add more switches...
           * else if (...)
           *  {
           *    ...
           *  }
           */
                else
                  {
                    clog << "warning: ignoring option '" << option << "'!" << endl;
                  }
        }
      else
        {
          string argument = token;
          /* Here you may add more argument handlers... */
          if (!MacFileName.size ())
            {
              MacFileName = token;
            }
          else if (srf == "default.root")
            {
              srf=token;
            }
          else
            {
              clog << "warning: ignoring argument '" << argument << "'!" << endl;
            }
        }
      iarg++;
    }

  try {

    // Creation of the analysis manager as soon as possible to be sure to get the daugther
    KM3AnalysisManager* analysis = K40AnalysisManager::getIt();
    analysis->fKeepAll=keepall;
    analysis->fFastAnalysisOnly=!analyse_all;
    analysis->Init (srf.data ());
    analysis->fVerbose=(short)debug;

    // random engine
    //CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
    //CLHEP::HepRandom::setTheEngine(new CLHEP::Hurd160Engine);
    CLHEP::HepRandom::setTheEngine(new CLHEP::DualRand);
    //CLHEP::HepRandom::restoreEngineStatus();
    //G4int ttime = time(0);
    //G4cout << "ttime=" << ttime << G4endl;
    timeval tim;
    gettimeofday(&tim, NULL);
    long ttime = tim.tv_sec xor tim.tv_usec;
    G4cout << "ttime=" << ttime << G4endl;

    G4int pid = getpid();
    G4cout << "pid=" << pid << G4endl;

    long Seed = ttime xor (pid << 8);
    // automatic (time-based) random seeds and filenames for each run
    G4cout << "******************" << G4endl;
    G4cout << "*** AUTOSEED ON ***" << G4endl;
    G4cout << "*******************" << G4endl;
    long seeds[2];
    seeds[0] =  Seed;
    seeds[1] =  (long)(Seed*G4UniformRand());
    G4cout << "seed1: " << seeds[0] << "; seed2: " << seeds[1] << G4endl;

    CLHEP::HepRandom::getTheEngine()->setSeeds(seeds,-1);
    //G4cout << "Index: " << CLHEP::HepRandom::getTheSeed() << G4endl;

    //long table[8];
    //CLHEP::HepRandom::getTheTableSeeds(table, 0);
    //for (int i=0; i< sizeof(table); i++)
    //G4cout << table[i] << G4endl;

    //G4cout << "****" << G4endl;

    //CLHEP::HepRandom::getTheTableSeeds(table, -1);
    //for (int i=0; i< sizeof(table); i++)
    // G4cout << table[i] << G4endl;


    //  CLHEP::HepRandom::getTheEngine()->setSeeds(seeds,-1);
    CLHEP::HepRandom::showEngineStatus();
    //return 0;
    //CLHEP::HepRandom::showEngineStatus();

    // table = CLHEP::HepRandom::getTheSeeds();
    // for (int i=0; i< sizeof(table); i++)
    //   G4cout << table[i] << G4endl;

    G4RunManager* runManager = new KM3RunManager;

    // set mandatory initialization classes

    K40DetectorConstruction* Detector = new K40DetectorConstruction;
    Detector->SetOrientationTheta (90* deg );
    runManager->SetUserInitialization(Detector);
    runManager->SetUserInitialization(new KM3PhysicsList);

    KM3PrimaryGeneratorAction* PGA=new KM3PrimaryGeneratorAction;
    PGA->SetProducingVolumeRadius(Detector->GetTargetFullLength ());

    // set mandatory user action class
    runManager->SetUserAction(new KM3PrimaryGeneratorAction);
    runManager->SetUserAction(new KM3RunAction);
    runManager->SetUserAction(new KM3EventAction);
    runManager->SetUserAction(new KM3SteppingAction);


    /*((KM3PrimaryGeneratorAction*)runManager->GetUserPrimaryGeneratorAction())->
      SetDetectorConstructor(Detector);
    */
    //
    //
    G4UIsession* session=0;
    G4VisManager* visManager=NULL;
    if (bsession)
      {
        session = new G4UIXm(0,argv_);//new G4UIterminal();
        // visualization manager
        visManager = new G4VisExecutive;
        visManager->Initialize();
      }

    // Initialize G4 kernel
    // do this at run time so the geometry/physics can be changed
    //  runManager->Initialize();

    // get the pointer to the User Interface manager
    G4UImanager* UI = G4UImanager::GetUIpointer();
    G4String command = "/control/execute ";

    if (session)   // Define UI session for interactive mode.
      {
        UI->ApplyCommand(command+MacFileName);
        session->SessionStart();
        delete session;
      }
    else           // Batch mode
      {
        UI->ApplyCommand(command+MacFileName);
        //session->SessionStart ();
      }
    analysis->Write ();
    //CLHEP::HepRandom::showEngineStatus();
    //CLHEP::HepRandom::getTheTableSeeds(table, 0);
    //for (int i=0; i< sizeof(table); i++)
    //  G4cout << table[i] << G4endl;
    //CLHEP::HepRandom::saveEngineStatus();
    // job termination
    if (visManager)
      {
        delete visManager;
      }
    cout << "baaaaaaaaaaaaah" << endl;

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(-17*2.54/2,17*2.54/2);

    delete runManager;
    cout << "baaaaaaaaaaaaah" << endl;
  }
  catch (const char* err)
    {
      cerr << "Exception occured: " << err << endl;
    }

  cout << "baaaaaaaaaaaaah" << endl;
  exit (0);

}
