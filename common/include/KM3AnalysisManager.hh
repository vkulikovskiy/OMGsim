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
#ifndef KM3AnalysisManager_h
#define KM3AnalysisManager_h 1

//---------------------------------------------------------------------------
//
// ClassName:   KM3AnalysisManager
//
// Description: Singleton class to hold analysis parameters and build histograms.
//              User cannot access to the constructor.
//              The pointer of the only existing object can be got via
//              KM3AnalysisManager::GetInstance() static method.
//              The first invokation of this static method makes
//              the singleton object.
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "globals.hh"
#include <vector>
#include <string>
#include "G4Track.hh"
#include <G4Event.hh>
#include <time.h>
#include "KM3DetectorConstruction.hh"
#include <KM3TrackInformation.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
class TTree;
class TFile;
class G4Event;

namespace km3net
{

  class KM3AnalysisManager
  {

  public:
    // With description

    static KM3AnalysisManager* getIt();
    virtual void Init (std::string RF = "");
    virtual void Write () = 0;

  protected:

    KM3AnalysisManager();
    virtual ~KM3AnalysisManager();

  public: // Without description
    virtual void SetOriginOfEvent (const G4Event* EV);
    virtual void AddDetectedParticule (const KM3TrackInformation* info, unsigned char LIM) = 0;
    virtual void AddPhotoTubeEntranceInfo (const G4Step* fStep);
    virtual void AddOMEntranceInfo (const G4Step* Step);
    virtual void EndOfEvent ();
    virtual void SetEnergy (double);

    virtual bool HasBeenInTheOMYet ();

    //virtual void AddAbranch (std::string TreeName, void* address, std::string branchList);

    static bool fKeepAll;
    static bool fFastAnalysisOnly;
    static short fVerbose;
    static std::string fRootFileName;

    std::string fOMBasePrefix;

  protected:
    // MEMBERS
    static KM3AnalysisManager* fManager;

    //the trees
    TFile* fRFile=NULL;
    TTree* fTRunInfos;
    TTree* fTEvent;
    TTree* fTFastAnalysis;


    unsigned long fNbOM=0;
    unsigned long fNbPhoto=0;


    //here some data to put in the trees
    unsigned long fNbTrigg=0;

    float fVertexX=0;
    float fVertexY=0;
    float fVertexZ=0;
    float fVertexT=0;
    float fVertexP=0;
    float fVertexDirX=0;
    float fVertexDirY=0;
    float fVertexDirZ=0;
    float fVertexEnergy=0;
    float fOMHitX=0;
    float fOMHitY=0;
    float fOMHitZ=0;
    float fOMHitT=0;
    float fOMHitP=0;
    float fPhotoX=0;
    float fPhotoY=0;
    float fPhotoZ=0;
    float fPhotoT=0;
    float fPhotoP=0;
    short fOMId=-1;
    short fPMID=0;

  };

}
#endif
