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
#ifndef AAAnalysisManager_h
#define AAAnalysisManager_h 1

//---------------------------------------------------------------------------
//
// ClassName:   AAAnalysisManager
//
// Description: Singleton class to hold analysis parameters and build histograms.
//              User cannot access to the constructor.
//              The pointer of the only existing object can be got via
//              AAAnalysisManager::GetInstance() static method.
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
#include <time.h>
#include <KM3AnalysisManager.hh>
#include "KM3DetectorConstruction.hh"
#include <KM3TrackInformation.hh>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

namespace km3net
{

  class AAAnalysisManager : public KM3AnalysisManager
  {

  public:
    // With description

    static AAAnalysisManager* getIt();
    virtual void Init (std::string RF = "");
    virtual void Write ();

  protected:

    AAAnalysisManager();
    ~AAAnalysisManager();

  public: // Without description
    virtual void AddDetectedParticule (const KM3TrackInformation* INFO, unsigned char LIM);
    virtual void EndOfEvent ();
    virtual void SetEnergy (double E);

  protected:
    size_t fcounter=0;
    size_t funcounter=0;
    size_t fomcounter=0;

    // MEMBERS

    unsigned fNbDetected[4][100]; //[limit][energies] energies in [1-5] eV

    float fTotalTime=0;
    float fEnergy=0;
    bool fPolarity =0;
    unsigned char fLimit=0;

  };

}
#endif
