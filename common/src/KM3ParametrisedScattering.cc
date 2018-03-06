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
// $Id: KM3ParametrisedScattering.cc,v 1.17 2008/10/24 19:51:12 gum Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
//
////////////////////////////////////////////////////////////////////////
// Optical Photon Rayleigh Scattering Class Implementation
////////////////////////////////////////////////////////////////////////
//
// File:        KM3ParametrisedScattering.cc
// Description: Discrete Process -- Rayleigh scattering of optical
//		photons
// Version:     1.0
// Created:     1996-05-31
// Author:      Juliet Armstrong
// Updated:     2005-07-28 - add G4ProcessType to constructor
//              2001-10-18 by Peter Gumplinger
//              eliminate unused variable warning on Linux (gcc-2.95.2)
//              2001-09-18 by mma
//		>numOfMaterials=G4Material::GetNumberOfMaterials() in BuildPhy
//              2001-01-30 by Peter Gumplinger
//              > allow for positiv and negative CosTheta and force the
//              > new momentum direction to be in the same plane as the
//              > new and old polarization vectors
//              2001-01-29 by Peter Gumplinger
//              > fix calculation of SinTheta (from CosTheta)
//              1997-04-09 by Peter Gumplinger
//              > new physics/tracking scheme
// mail:        gum@triumf.ca
//
////////////////////////////////////////////////////////////////////////

#include "G4ios.hh"
#include "G4OpProcessSubType.hh"

#include "KM3ParametrisedScattering.hh"
//for root debug file
/*
  #include "TFile.h"
  #include "TGraph.h"
  #include "TMath.h"
  #include "TObject.h"
  #include "TH1.h"
  #include "Randomize.hh"
*/
#include  "CLHEP/Random/RandGeneral.h"
#include <iostream>
using namespace std;
using namespace CLHEP;
namespace km3net
{

  /**
   * Auxiliary method to describe light scattering in water (Henyey-Greenstein)
   *
   * \param  g      angular dependence parameter
   * \param  x      cosine scattering angle
   * \return        probability
   */

  static inline double henyey_greenstein(const double& g,
                                         const double& x)

  {
    static const double a0 = (1.0 - g*g) / (4*pi);

    const double y = 1.0 + g*g - 2.0*g*x;

    return a0 / (y*sqrt(y));
  }


  /**
   * Auxiliary method to describe light scattering in water (Rayleigh or more precise Einstein-Smoluchowski)
   *
   * \param  a      angular dependence parameter
   * \param  x      cosine scattering angle
   * \return        probability
   */
  static inline double rayleigh(const double& a,
                                const double& x)

  {
    static const double a0 = 1.0 / (1.0 + a/3.0) / (4*pi);

    return a0 * (1.0 + a*x*x);
  }


  /**
   * Auxiliary method to describe light scattering in water (Rayleigh)
   *
   * \param  x      cosine scattering angle
   * \return        probability
   */
  static inline double rayleigh(const double& x)

  {
    return rayleigh(0.853, x);
  }


  /**
   * Model specific function to describe light scattering in water (p00075)
   *
   * \param  x      cosine scattering angle
   * \return        probability
   */
  static inline double p00075(const double& x)
  {
    static const double g = 0.924;
    static const double f = 0.17;

    return f * rayleigh(x)  +  (1.0 - f) * henyey_greenstein(g,x);
  }
  /////////////////////////
  // Class Implementation
  /////////////////////////

  //////////////
  // Operators
  //////////////

  // KM3ParametrisedScattering::operator=(const KM3ParametrisedScattering &right)
  // {
  // }

  /////////////////
  // Constructors
  /////////////////



  KM3ParametrisedScattering::KM3ParametrisedScattering(const G4String& processName, G4ProcessType type)
  : G4VDiscreteProcess(processName, type)
  {
    SetProcessSubType(fOpRayleigh);

    thePhysicsTable = 0;

    DefaultWater = false;

    if (verboseLevel>0) {
      G4cout << GetProcessName() << " is created " << G4endl;
    }

    SetSmithBacker ();

    BuildThePhysicsTable();
  }

  // KM3ParametrisedScattering::KM3ParametrisedScattering(const KM3ParametrisedScattering &right)
  // {
  // }

  ////////////////
  // Destructors
  ////////////////

  KM3ParametrisedScattering::~KM3ParametrisedScattering()
  {
    if (thePhysicsTable!= 0) {
      thePhysicsTable->clearAndDestroy();
      delete thePhysicsTable;
      delete[] fprobList;
    }
  }

  ////////////
  // Methods
  ////////////

  // PostStepDoIt
  // -------------
  //

  void
  KM3ParametrisedScattering::SetSmithBacker (bool B)
  {
    if (!fprobList)
      {
        fprobList = new double[fnBins];
      }
    if (B)
      {
        for (unsigned ii = 0; ii < fnBins; ii++)
          {
            double x = pi*ii/(fnBins-1);
            fprobList[ii]=p00075(cos(x))*sin(x);
            //cout << "x = " << x << " prob " << probList[ii] << endl;
          }
      }
    else
      {
        for (unsigned ii = 0; ii < fnBins; ii++)  {
          double x = pi*ii/(fnBins-1);
          fprobList[ii]=rayleigh(cos(x))*sin(x);
        }
      }
  }

  G4VParticleChange*
  KM3ParametrisedScattering::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
  {
    aParticleChange.Initialize(aTrack);

    const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();

    if (verboseLevel>0) {
      G4cout << "Scattering Photon!" << G4endl;
      G4cout << "Old Momentum Direction: "
             << aParticle->GetMomentumDirection() << G4endl;
      G4cout << "Old Polarization: "
             << aParticle->GetPolarization() << G4endl;
    }
    const G4Material* aMaterial = aTrack.GetMaterial();
    if (aMaterial->GetName() == "AntaresWater" || "NEMOWater") {
      G4double phidiff = 2*pi*G4UniformRand();

      CLHEP::RandGeneral GenDist(fprobList,fnBins);

      double thetadiff = pi*GenDist.shoot();  // shoots values using the engine
      // in the static generator. shoot()
      // provides the same functionality
      // of fire() in this case.
      //double thetadiff = pi;
      G4ThreeVector OldMomentumDirection = aParticle->GetMomentumDirection();
      G4ThreeVector NewMomentumDirection = aParticle->GetMomentumDirection();
      //idea - first to generate vector with |theta_new-theta_old|=thetadiff
      //after rotate a vector with (theta_new,phi) around old vector
      double thetanew = OldMomentumDirection.theta()-thetadiff;
      if (thetanew<0 || thetanew>=pi) thetanew = OldMomentumDirection.theta()+thetadiff;
      NewMomentumDirection.setTheta(thetanew); //change theta on
      NewMomentumDirection.rotate(phidiff,OldMomentumDirection);
      aParticleChange.ProposeMomentumDirection(NewMomentumDirection);
    }
    else if (aMaterial->GetName() == "Water") {   //internal gean4 Rayleght
      // find polar angle w.r.t. old polarization vector

      G4double rand = G4UniformRand();

      G4double CosTheta = std::pow(rand, 1./3.);
      G4double SinTheta = std::sqrt(1.-CosTheta*CosTheta);

      if(G4UniformRand() < 0.5)CosTheta = -CosTheta;

      // find azimuthal angle w.r.t old polarization vector

      rand = G4UniformRand();

      G4double Phi = twopi*rand;
      G4double SinPhi = std::sin(Phi);
      G4double CosPhi = std::cos(Phi);

      G4double unit_x = SinTheta * CosPhi;
      G4double unit_y = SinTheta * SinPhi;
      G4double unit_z = CosTheta;

      G4ThreeVector NewPolarization (unit_x,unit_y,unit_z);

      // Rotate new polarization direction into global reference system

      G4ThreeVector OldPolarization = aParticle->GetPolarization();
      OldPolarization = OldPolarization.unit();

      NewPolarization.rotateUz(OldPolarization);
      NewPolarization = NewPolarization.unit();

      // -- new momentum direction is normal to the new
      // polarization vector and in the same plane as the
      // old and new polarization vectors --

      G4ThreeVector NewMomentumDirection =
        OldPolarization - NewPolarization * CosTheta;

      if(G4UniformRand() < 0.5)NewMomentumDirection = -NewMomentumDirection;
      NewMomentumDirection = NewMomentumDirection.unit();

      aParticleChange.ProposePolarization(NewPolarization);

      aParticleChange.ProposeMomentumDirection(NewMomentumDirection);

      if (verboseLevel>0) {
        G4cout << "New Polarization: "
               << NewPolarization << G4endl;
        G4cout << "Polarization Change: "
               << *(aParticleChange.GetPolarization()) << G4endl;
        G4cout << "New Momentum Direction: "
               << NewMomentumDirection << G4endl;
        G4cout << "Momentum Change: "
               << *(aParticleChange.GetMomentumDirection()) << G4endl;
      }
    }

    /*
      TFile *f1 = new TFile("debug.root", "UPDATE");
      TH1I *histo_theta = (TH1I *)f1->Get("theta_distr");
      if (histo_theta==0) {
      histo_theta = new TH1I("theta_distr","theta_distr",3600,-180.,180.);
      }
      TH1I *histo_step = (TH1I *)f1->Get("step_distr");
      if (histo_step==0) {
      histo_step = new TH1I("step_distr","step_distr",10000,0.,1000.);
      }
      double theta = (aParticle->GetMomentumDirection().angle(*(aParticleChange.GetMomentumDirection())))/pi*180.;
      histo_theta->Fill(theta);
      histo_theta->Write(histo_theta->GetName(),TObject::kOverwrite);
      histo_step->Fill(aStep.GetStepLength()/m);
      histo_step->Write(histo_step->GetName(),TObject::kOverwrite);
      f1->Close();
    */


    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }

  // BuildThePhysicsTable for the Rayleigh Scattering process
  // --------------------------------------------------------
  //
  void KM3ParametrisedScattering::BuildThePhysicsTable()
  {
    //      Builds a table of scattering lengths for each material

    if (thePhysicsTable) return;

    const G4MaterialTable* theMaterialTable=
      G4Material::GetMaterialTable();
    G4int numOfMaterials = G4Material::GetNumberOfMaterials();

    // create a new physics table

    thePhysicsTable = new G4PhysicsTable(numOfMaterials);

    // loop for materials

    for (G4int i=0 ; i < numOfMaterials; i++)
      {
        G4PhysicsOrderedFreeVector* ScatteringLengths =
          new G4PhysicsOrderedFreeVector();

        G4MaterialPropertiesTable *aMaterialPropertiesTable =
          (*theMaterialTable)[i]->GetMaterialPropertiesTable();

        if(aMaterialPropertiesTable){

          G4MaterialPropertyVector* AttenuationLengthVector =
            aMaterialPropertiesTable->GetProperty("RAYLEIGH");

          if(!AttenuationLengthVector){

            if ((*theMaterialTable)[i]->GetName() == "Water")
              {
                // Call utility routine to Generate
                // Rayleigh Scattering Lengths

                DefaultWater = true;

                ScatteringLengths =
                  RayleighAttenuationLengthGenerator(aMaterialPropertiesTable);
              }
          }// else ScatteringLengths = aMaterialPropertiesTable->GetProperty("RAYLEIGH"); //here!!!
        }

        thePhysicsTable->insertAt(i,ScatteringLengths);
      }
  }

  // GetMeanFreePath()
  // -----------------
  //
  G4double KM3ParametrisedScattering::GetMeanFreePath(const G4Track& aTrack,
                                                      G4double ,
                                                      G4ForceCondition* )
  {
    const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
    const G4Material* aMaterial = aTrack.GetMaterial();

    G4double thePhotonEnergy = aParticle->GetTotalEnergy();

    G4double AttenuationLength = DBL_MAX;

    if (aMaterial->GetName() == "Water" && DefaultWater){

      G4bool isOutRange;

      AttenuationLength =
        (*thePhysicsTable)(aMaterial->GetIndex())->
        GetValue(thePhotonEnergy, isOutRange);
    }
    else {

      G4MaterialPropertiesTable* aMaterialPropertyTable =
        aMaterial->GetMaterialPropertiesTable();

      if(aMaterialPropertyTable){
        G4MaterialPropertyVector* AttenuationLengthVector =
          aMaterialPropertyTable->GetProperty("RAYLEIGH");
        if(AttenuationLengthVector){
          AttenuationLength = AttenuationLengthVector ->
            Value(thePhotonEnergy);
        }
        else{
          //               G4cout << "No Rayleigh scattering length specified" << G4endl;
        }
      }
      else{
        //             G4cout << "No Rayleigh scattering length specified" << G4endl;
      }
    }

    return AttenuationLength;
  }

  // RayleighAttenuationLengthGenerator()
  // ------------------------------------
  // Private method to compute Rayleigh Scattering Lengths (for water)
  //
  G4PhysicsOrderedFreeVector*
  KM3ParametrisedScattering::RayleighAttenuationLengthGenerator(G4MaterialPropertiesTable *aMPT)
  {
    // Physical Constants

    // isothermal compressibility of water
    G4double betat = 7.658e-23*m3/MeV;

    // K Boltzman
    G4double kboltz = 8.61739e-11*MeV/kelvin;

    // Temperature of water is 10 degrees celsius
    // conversion to kelvin:
    // TCelsius = TKelvin - 273.15 => 273.15 + 10 = 283.15
    G4double temp = 283.15*kelvin;

    // Retrieve vectors for refraction index
    // and photon energy from the material properties table

    G4MaterialPropertyVector* Rindex = aMPT->GetProperty("RINDEX");

    G4double refsq;
    G4double e;
    G4double xlambda;
    G4double c1, c2, c3, c4;
    G4double Dist;
    G4double refraction_index;

    G4PhysicsOrderedFreeVector *RayleighScatteringLengths =
      new G4PhysicsOrderedFreeVector();
    /*
      TFile *f1 = new TFile("debug.root", "UPDATE");
      TGraph *gr1 = new TGraph();
      gr1->SetName("scat_length");
    */


    if (Rindex ) {

      for (size_t indexit=0 ; indexit < Rindex->GetVectorLength () ; indexit++) {

        e = Rindex->Energy(indexit);

        refraction_index = (*Rindex)[indexit];
        refsq = refraction_index*refraction_index;
        xlambda = h_Planck*c_light/e;

        if (verboseLevel>0) {
          G4cout << Rindex->Energy(indexit) << " MeV\t";
          G4cout << xlambda << " mm\t";
        }

        c1 = 1 / (6.0 * pi);
        c2 = std::pow((2.0 * pi / xlambda), 4);
        c3 = std::pow( ( (refsq - 1.0) * (refsq + 2.0) / 3.0 ), 2);
        c4 = betat * temp * kboltz;

        Dist = 1.0 / (c1*c2*c3*c4);

        if (verboseLevel>0) {
          G4cout << Dist << " mm" << G4endl;
        }
        RayleighScatteringLengths->
          InsertValues(Rindex->Energy(indexit), Dist);
        /*
          int npoint = gr1->GetN();
          gr1->SetPoint(npoint,197.3*MeV*fermi*2.*pi/(Rindex->GetPhotonEnergy())/um, Dist/m);
        */
      }

    }
    /*
      gr1->Write();
      f1->Close();
    */

    return RayleighScatteringLengths;
  }
}
