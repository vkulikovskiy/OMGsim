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


//
// **********************************************************************

#ifndef KM3PrimaryGeneratorAction_h
#define KM3PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "KM3DetectorConstruction.hh"

#include <G4SystemOfUnits.hh>

#include <vector>

class G4GeneralParticleSource;
class G4ParticleGun;
class G4Event;
class G4Navigator;
class G4VPhysicalVolume;
class G4VPrimaryGenerator;
class TH1F;

namespace km3net
{

  class KM3PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
  {
    friend class KM3RunMessenger;
  public:
    KM3PrimaryGeneratorAction();
    ~KM3PrimaryGeneratorAction();


  public:
    void GeneratePrimaries(G4Event* anEvent);

    void SetProducingVolumeRadius (double RAD) {fTargetRadius=RAD;}

  private:
    void InitializeNormalType ();
    void InitializeDecayType ();
    G4int PrevTime;

  private:
    G4VPrimaryGenerator* particleGun;
    //  G4ParticleGun* particleGun;
    G4int randcount;
    G4Navigator* Navigator;
    G4double fTargetRadius = 500.*m;
    unsigned char fGunType=0;//1 for normal, 2 for decay

    std::vector <float>         fLEDAngle;

#if __GNUC__>=5
    ::TH1F* fhLED;
#else
    TH1F* fhLED;
#endif
  };
}
#endif
