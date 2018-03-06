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

// **********************************************************************

#include <G4Neutron.hh>

#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleGun.hh"
#include "Randomize.hh"
#include <G4Navigator.hh>
#include <G4VPhysicalVolume.hh>
#include <G4TransportationManager.hh>

#include "KM3PrimaryGeneratorAction.hh"
#include <KM3DetectorConstruction.hh>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;
using namespace CLHEP;
namespace km3net
{

  KM3PrimaryGeneratorAction::KM3PrimaryGeneratorAction()
  {
    particleGun = 0;//new G4GeneralParticleSource();
  }

  KM3PrimaryGeneratorAction::~KM3PrimaryGeneratorAction()
  {
    if (particleGun != NULL)
      delete particleGun;
  }

  void KM3PrimaryGeneratorAction::InitializeNormalType ()
  {
    if (fGunType)
      throw ("Two guntype init in the mac file !!");

    cout << "Initializing normal type" << endl;

    particleGun = new G4GeneralParticleSource ();
    fGunType = 1;
  }

  void KM3PrimaryGeneratorAction::InitializeDecayType ()
  {
    if (fGunType)
      throw ("Two guntype init in the mac file !!");

    fTargetRadius=KM3DetectorConstruction::GetIt()->GetTargetFullLength();

    particleGun = new G4ParticleGun ();
    randcount = 0;
    Navigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
    fGunType = 2;
  }

  void KM3PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
  {
    string toprint;
    if (fGunType == 1)
      {
        particleGun->GeneratePrimaryVertex(anEvent);
        //cout << "Event generated at: " << anEvent->GetPrimaryVertex()->GetPosition().mag() << endl;
        return;
      }

    if (fGunType == 2)
      {
        G4VPhysicalVolume* initialVolume = NULL;

        G4ThreeVector InitPos;

        do
          {
            InitPos.setRThetaPhi(fTargetRadius * pow(G4UniformRand(),1./3),
                                 acos(2*G4UniformRand()-1),
                                 twopi*G4UniformRand());

            initialVolume = Navigator->LocateGlobalPointAndSetup(InitPos);

            initialVolume?toprint=initialVolume->GetName ():toprint="NULL!";
          }
        while (initialVolume == NULL || initialVolume->GetName () != "Target");



        particleGun->SetParticlePosition (InitPos);


        ((G4GeneralParticleSource*)particleGun)->GeneratePrimaryVertex(anEvent);


        return;

      }

    string tosend("Gun type not initialized in mac file: fGunType");
    tosend+=(char)(fGunType+91);

    throw (tosend.data ());

  }
}
