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
#ifndef WPropAnalysisManager_h
#define WPropAnalysisManager_h 1

//---------------------------------------------------------------------------
//
// ClassName:   WPropAnalysisManager
//
// Description: Singleton class to hold analysis parameters and build histograms.
//              User cannot access to the constructor.
//              The pointer of the only existing object can be got via
//              WPropAnalysisManager::GetInstance() static method.
//              The first invokation of this static method makes
//              the singleton object.
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "globals.hh"
#include <vector>
#include <list>
#include <vector>
#include <string>
#include "G4Track.hh"
#include <G4ThreeVector.hh>
#include <time.h>
#include <KM3AnalysisManager.hh>
#include "KM3DetectorConstruction.hh"
#include <KM3PetzoldScattering.hh>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
class TH1F;
class TRandom;

namespace km3net
{

  class WPropAnalysisManager : public KM3AnalysisManager
  {

  public:
    // With description

    static WPropAnalysisManager* getIt();
    virtual void Init (std::string RF = "");
    virtual void Write ();

  private:

    WPropAnalysisManager();
    ~WPropAnalysisManager();

  public: // Without description
    virtual void AddDetectedParticule (const KM3TrackInformation* INFO, unsigned char LIM);
    virtual void EndOfEvent ();
    virtual void SetOriginOfEvent (const G4Event* EV);

    virtual void IncrementProcessWeight (KM3PetzoldScattering::PROCESS PR);

  private:
    size_t fcounter=0;
    size_t funcounter=0;
    size_t fomcounter=0;

    // MEMBERS
    //static WPropAnalysisManager* fManager;
    //the trees

    float         fEmmitingDirX  =0;
    float         fEmmitingDirY  =0;
    float         fEmmitingDirZ  =0;
    float         fFloorCrossPosX=0;
    float         fFloorCrossPosY=0;
    float         fFloorCrossPosZ=0;
    float         fFloorCrossDirX=0;
    float         fFloorCrossDirY=0;
    float         fFloorCrossDirZ=0;
    float         fTotalTime     =0;
    float         fEnergy        =0;
    float         fWeigthANTARES =0;
    float         fWeigthNEMO    =0;
    unsigned      fESWeight      =0;
    unsigned      fKopelWeight   =0;
    unsigned char fFloorNbCross  =0;
    size_t        fLEDDetected[2]; //from 0 to 12 deg, 0.1 steps in binary

    std::vector <float>         fLEDAngle;
    std::vector <float>         fEffNEMO;
    std::vector <float>         fEffANTARES;

    TRandom *fRandom;

  private:
#if __GNUC__>=5
    ::TH1F* fhLED;
#else
    TH1F* fhLED;
#endif
    G4ThreeVector fRealEmmiting;


  };

}
#endif
