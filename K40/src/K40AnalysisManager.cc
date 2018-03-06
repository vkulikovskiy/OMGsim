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

#include "K40AnalysisManager.hh"
#include "G4UnitsTable.hh"
#include "G4VTouchable.hh"
#include <G4NavigationHistory.hh>
#include <KM3DetectorConstruction.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

namespace km3net
{
  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

  K40AnalysisManager* K40AnalysisManager::getIt()
  {
    if(!fManager) {
      fManager = new K40AnalysisManager();
    }
    return (K40AnalysisManager*)fManager;
  }

  void K40AnalysisManager::SetTargetVolume (float V)
  {
    fTargetVolume=V;
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

  K40AnalysisManager::K40AnalysisManager()
  {

  }

  void K40AnalysisManager::Init (string RFNAME)
  {
    AAAnalysisManager::Init (RFNAME);
    fTFastAnalysis->Branch("FastAnalysis",&fCoincidence,"Coincidences[7][3]/i");
    fTEvent->Branch("HitNumber",&fNbTrigg,"EventNb/l");
    fTRunInfos->Branch("TargetVolume", &fTargetVolume, "TargetVolume/f");
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

  K40AnalysisManager::~K40AnalysisManager()
  {
  }

  void K40AnalysisManager::Write  ()
  {
    float targetlength=KM3DetectorConstruction::GetIt()->GetTargetFullLength ()/m;
    fTargetVolume=pow(targetlength,3)*M_PI*4./3;
    AAAnalysisManager::Write ();
  }

  void K40AnalysisManager::AddDetectedParticule (const KM3TrackInformation* INFO, unsigned char LIM)
  {
    AAAnalysisManager::AddDetectedParticule (INFO, LIM);
    fListOfOM.push_back(fOMId);
    if (fSmallestLimit > LIM - 1)
      fSmallestLimit=LIM - 1;
  }

  void K40AnalysisManager::EndOfEvent ()
  {
    if (fRFile == NULL)
      return;

    fListOfOM.sort ();
    fListOfOM.unique ();


    //1?       <=> 0
    //11 12    <=> 1
    //11 13    <=> 2
    //12 13    <=> 3
    //11 12 13 <=> 4
    //2?       <=> 5
    //21 22    <=> 6


    if (fListOfOM.size () == 1)
      {
        if (fListOfOM.front () < 20)
          {
            fCoincidence[0][(unsigned)(fSmallestLimit)]++;
          }
        else
          {
            fCoincidence[5][(unsigned)(fSmallestLimit)]++;
          }
      }
    else if (fListOfOM.size () == 2)
      {
        if (fListOfOM.front () == 21 && fListOfOM.back () == 22)
          {
            fCoincidence[6][(unsigned)(fSmallestLimit)]++;
          }
        if (fListOfOM.front () == 11)
          {
            if (fListOfOM.back () == 12)
              {
                fCoincidence[1][(unsigned)(fSmallestLimit)]++;
              }
            else if (fListOfOM.back () == 13)
              {
                fCoincidence[2][(unsigned)(fSmallestLimit)]++;
              }
          }
        else if (fListOfOM.front () == 12 && fListOfOM.back () == 13)
          {
            fCoincidence[3][(unsigned)(fSmallestLimit)]++;
          }
      }
    else if (fListOfOM.size () == 3 &&
             fListOfOM.front () == 11 && fListOfOM.back () == 13 && *(++(fListOfOM.begin ())) == 12)
      {
        fCoincidence[4][(unsigned)(fSmallestLimit)]++;
      }

    fSmallestLimit=3;
    fListOfOM.clear ();

    AAAnalysisManager::EndOfEvent ();

  }

}

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
