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
#ifndef K40AnalysisManager_h
#define K40AnalysisManager_h 1

//---------------------------------------------------------------------------
//
// ClassName:   K40AnalysisManager
//
// Description: Singleton class to hold analysis parameters and build histograms.
//              User cannot access to the constructor.
//              The pointer of the only existing object can be got via
//              K40AnalysisManager::GetInstance() static method.
//              The first invokation of this static method makes
//              the singleton object.
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "globals.hh"
#include <vector>
#include <list>
#include <string>
#include "G4Track.hh"
#include <time.h>
#include <AAAnalysisManager.hh>
#include "KM3DetectorConstruction.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

namespace km3net
{

  class K40AnalysisManager : public AAAnalysisManager
  {

  public:
    // With description

    static K40AnalysisManager* getIt();
    virtual void Init (std::string RF = "");
    virtual void Write ();

  protected:

    K40AnalysisManager();
    ~K40AnalysisManager();

  public: // Without description
    virtual void AddDetectedParticule (const KM3TrackInformation* INFO, unsigned char LIM);
    virtual void EndOfEvent ();
    virtual void SetTargetVolume (float V);

  private:
    //for run info
    float         fTargetVolume=0;

    // MEMBERS

    unsigned fCoincidence[7][4]={{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}}; //[Combinaison][LIM]
    //OM Antares single,  12 13 23 123, Nemo single 45 (45 for Nemo)

    std::list<unsigned short> fListOfOM;
    char fSmallestLimit=3;

  };

}
#endif
