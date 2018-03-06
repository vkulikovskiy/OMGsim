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
// $Id: KM3ParametrisedScattering.hh,v 1.9 2006/06/29 21:08:40 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
//
////////////////////////////////////////////////////////////////////////
// Optical Photon Rayleigh Scattering Class Definition
////////////////////////////////////////////////////////////////////////
//
// File:        KM3ParametrisedScattering.hh
// Description: Discrete Process -- Rayleigh scattering of optical photons
// Version:     1.0
// Created:     1996-05-31
// Author:      Juliet Armstrong
// Updated:     2005-07-28 add G4ProcessType to constructor
//              1999-10-29 add method and class descriptors
//              1997-04-09 by Peter Gumplinger
//              > new physics/tracking scheme
// mail:        gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#ifndef KM3ParametrisedScattering_h
#define KM3ParametrisedScattering_h 1

/////////////
// Includes
/////////////

#include "globals.hh"
#include "templates.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleMomentum.hh"
#include "G4Step.hh"
#include "G4VDiscreteProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4OpticalPhoton.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsOrderedFreeVector.hh"

// Class Description:
// Discrete Process -- Rayleigh scattering of optical photons.
// Class inherits publicly from G4VDiscreteProcess.
// Class Description - End:

/////////////////////
// Class Definition
/////////////////////
#include <G4SystemOfUnits.hh>

namespace km3net
{

  class KM3ParametrisedScattering : public G4VDiscreteProcess
  {
  private:

    //////////////
    // Operators
    //////////////

    // KM3ParametrisedScattering& operator=(const KM3ParametrisedScattering &right);

  public: // Without description

    ////////////////////////////////
    // Constructors and Destructor
    ////////////////////////////////

    KM3ParametrisedScattering(const G4String& processName = "ParametrisedScattering",
                              G4ProcessType type = fOptical);

    // KM3ParametrisedScattering(const KM3ParametrisedScattering &right);

    ~KM3ParametrisedScattering();

    ////////////
    // Methods
    ////////////

  public: // With description
    void SetSmithBacker (bool B=true);

    G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
    // Returns true -> 'is applicable' only for an optical photon.

    G4double GetMeanFreePath(const G4Track& aTrack,
                             G4double ,
                             G4ForceCondition* );
    // Returns the mean free path for Rayleigh scattering in water.
    // --- Not yet implemented for other materials! ---

    G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                    const G4Step&  aStep);
    // This is the method implementing Rayleigh scattering.

    G4PhysicsTable* GetPhysicsTable() const;
    // Returns the address of the physics table.

    void DumpPhysicsTable() const;
    // Prints the physics table.

  private:
    void SetPartDens (G4double pd)
    {
      PartDens=pd;
    }

    void BuildThePhysicsTable();

    /////////////////////
    // Helper Functions
    /////////////////////

    G4PhysicsOrderedFreeVector* RayleighAttenuationLengthGenerator(
                                                                   G4MaterialPropertiesTable *aMPT);

    ///////////////////////
    // Class Data Members
    ///////////////////////

  protected:

    G4PhysicsTable* thePhysicsTable;
    //  A Physics Table can be either a cross-sections table or
    //  an energy table (or can be used for other specific
    //  purposes).

  private:
    bool fSmithBacker=true;
    double* fprobList=NULL;
    unsigned fnBins=1000;
    G4bool DefaultWater;
    G4double PartDens=1.;

  };

  ////////////////////
  // Inline methods
  ////////////////////

  inline
  G4bool KM3ParametrisedScattering::IsApplicable(const G4ParticleDefinition& aParticleType)
  {
    return ( &aParticleType == G4OpticalPhoton::OpticalPhoton() );
  }

  inline
  void KM3ParametrisedScattering::DumpPhysicsTable() const

  {
    G4int PhysicsTableSize = thePhysicsTable->entries();
    G4PhysicsOrderedFreeVector *v;

    for (G4int i = 0 ; i < PhysicsTableSize ; i++ )
      {
        v = (G4PhysicsOrderedFreeVector*)(*thePhysicsTable)[i];
        v->DumpValues();
      }
  }

  inline G4PhysicsTable* KM3ParametrisedScattering::GetPhysicsTable() const
  {
    return thePhysicsTable;
  }

}
#endif /* KM3ParametrisedScattering_h */
