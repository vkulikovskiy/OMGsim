#ifndef KM3PETZOLD_h
#define KM3PETZOLD_h 1


#include "globals.hh"
#include "templates.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleMomentum.hh"
#include "G4Step.hh"
#include "G4VDiscreteProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4OpticalPhoton.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsOrderedFreeVector.hh"
#include <vector>

namespace km3net
{

  class KM3PetzoldScattering : public G4VDiscreteProcess
  {
    friend class KM3MaterialMessenger;

  public:

    KM3PetzoldScattering (const G4String& processName = "KM3PetzoldScattering",
                          G4ProcessType type = fOptical);

    ~KM3PetzoldScattering ();

  public:
    //process enum
    enum PROCESS
      {
        NOTHING=0,
        KOPELEVITCH,
        EINSTEINSMOLUCHOWSKI
      };

    ////////////
    // Methods
    ////////////

    G4double GetMeanFreePath(const G4Track& atrack,
                             G4double ,
                             G4ForceCondition* );

    G4bool IsApplicable(const G4ParticleDefinition& aParticleType);
    // Returns true -> 'is applicable' only for an optical photon.

    void BuildPhysicsTable(const G4ParticleDefinition& aParticleType);
    // Build table at a right time

    G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                    const G4Step&  aStep);
    // This is the method implementing Rayleigh scattering.

    G4PhysicsTable* GetPhysicsTable() const;
    // Returns the address of the physics table.

    PROCESS GetLastProcess () const {return fLastProcess;}

    void SetESFactor (double fact) {fESfactor=fact;}
    void SetKopelFactor (double fact) {fKopelFactor=fact;}

    double GetESFactor () {return fESfactor;}
    double GetKopelFactor () {return fKopelFactor;}
  private:

    void BuildThePhysicsTable();

    /////////////////////
    // Helper Functions
    /////////////////////

    G4PhysicsOrderedFreeVector* PetzoldAttenuationLengthGenerator();

    ///////////////////////
    // Class Data Members
    ///////////////////////

  private:
    //internal methods developped by V. Kulikovskyi
    void InitializePetzold();
    double GetPetzoldAngle(double r);
    double GetESAngle(double r);
    double GetESVolScat(double lambda);
    double GetKopelevichVolScat(double lambda);
    double GetTotalLength(double lambda);
    double GetAngle(double lambda);
    double GetEta(double *x, double *par);



  protected:

    G4PhysicsTable* thePhysicsTable;
    //  A Physics Table can be either a cross-sections table or
    //  an energy table (or can be used for other specific
    //  purposes).

  private:
    double fESfactor = 1;
    double fKopelFactor = 0.9;


    PROCESS fLastProcess=NOTHING;

    bool fpetzoldinitialised = false;
    G4bool DefaultWater;
    std::vector <double> part_rad;
    std::vector <double> part_beta;
    std::vector <double>  part_acu;



  };


}
#endif
