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


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
#include <TTree.h>
#include <TFile.h>

#include "AAAnalysisManager.hh"
#include "G4UnitsTable.hh"
#include "G4VTouchable.hh"
#include <G4NavigationHistory.hh>
#include <G4TransportationManager.hh>
#include <KM3DetectorConstruction.hh>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

namespace km3net
{
  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

  AAAnalysisManager* AAAnalysisManager::getIt()
  {
    if(!fManager) {
      fManager = new AAAnalysisManager();
    }
    return (AAAnalysisManager*)fManager;
  }


  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

  AAAnalysisManager::AAAnalysisManager()
  {
    memset (fNbDetected,0 , sizeof(unsigned)*4*100);
  }

  void AAAnalysisManager::Init (string RFNAME)
  {
    KM3AnalysisManager::Init (RFNAME);
    fTFastAnalysis->Branch ("Detection" ,fNbDetected,"NbDetected[4][100]/i");
    fTEvent->Branch        ("Hit"       ,&fTotalTime,"TotalTime:Energy:Polarity/O:Limit/b");
    }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

  AAAnalysisManager::~AAAnalysisManager()
  {

  }

  void AAAnalysisManager::Write  ()
  {
    if (fRFile == NULL)
      return;

    fTRunInfos->Fill ();
    fTFastAnalysis->Fill ();


    fRFile->Write ();
    fRFile->Close ();
    delete (fRFile);
    fRFile=0;
    cout << "Analysis file was normaly closed." << endl;
  }

  void AAAnalysisManager::AddDetectedParticule (const KM3TrackInformation* INFO, unsigned char LIM)
  {

    fEnergy=INFO->GetOriginalEnergy ();
    unsigned Erange=(unsigned)((10e-06 - fEnergy)/(10e-06)*100+0.5);

    if (Erange > 99
        //|| STEP->GetPostStepPoint ()->GetTouchable ()->GetCopyNumber (1) != 1
        )
      return;


    Erange=100-Erange;
    fNbDetected[LIM-1][Erange]++;
    if (fNbDetected[LIM - 1][Erange] % 100 ==1)
      {
        cout << "\r\033[60C" << fNbDetected[LIM - 1][Erange] << " \r\033[65C detected lim-Erange " << (int)LIM << '-' << Erange;
        cout.flush ();
      }


    if (fRFile == NULL)
      return;
    if (fFastAnalysisOnly)
      return;

    fPolarity = INFO->GetOriginalPolarization ().isParallel(INFO->GetOriginalMomentum(),0.1);


    fTotalTime = INFO->GetOriginalTime ();

    fLimit = LIM;

    fTEvent->Fill ();
  }

  void AAAnalysisManager::EndOfEvent ()
  {
    if (fRFile == NULL)
      return;
    if (fFastAnalysisOnly)
      return;
    // if (fKeepAll         ||
    //     (fPhotoX != 0.   ||
    //      fPhotoY != 0.   ||
    //      fPhotoZ != 0. ) ||
    //     (fOMHitX != 0.   ||
    //      fOMHitY != 0.   ||
    //      fOMHitZ != 0. ))
    // {
    //   fTEvent->Fill ();
    // }

    KM3AnalysisManager::EndOfEvent ();
    fTotalTime=.0;
    fEnergy=0;
    fLimit=0;
  }

  void AAAnalysisManager::SetEnergy (double E)
  {
    fEnergy=E;
  }
}

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
