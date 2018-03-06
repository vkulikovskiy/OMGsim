#include "KM3PetzoldScattering.hh"
#include "G4ios.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4OpProcessSubType.hh"

#include <fstream>


using namespace std;

namespace km3net
{


  KM3PetzoldScattering::KM3PetzoldScattering(const G4String& processName, G4ProcessType type) :
  G4VDiscreteProcess(processName, type)

  {

    SetProcessSubType(fOpRayleigh);

    thePhysicsTable = NULL;

    DefaultWater = false;

    if (verboseLevel>0) {
      G4cout << GetProcessName() << " is created " << G4endl;
    }

  }

  KM3PetzoldScattering::~KM3PetzoldScattering()
  {
    if (thePhysicsTable!= NULL) {
      thePhysicsTable->clearAndDestroy();
      delete thePhysicsTable;
    }
  }

  /*const int ndat_beta = 55;

   */

  void KM3PetzoldScattering::InitializePetzold (){
    if (!fpetzoldinitialised) {
      double tmprad,tmpbeta;

      ifstream ifs((string(PROJECT_SOURCE_DIR)+"/common/data/physics/petzold_phase.data").data ());
      if (!ifs.good ())
        {
          cerr << "The petzold_phase.data couldn't be open" << endl;
          exit (1);
        }
      while (ifs >> tmprad >> tmpbeta)
        {
          part_rad.push_back (tmprad);
          part_beta.push_back (tmpbeta);
        }

      part_rad[0] = part_rad[0]/180.*CLHEP::pi;
      part_beta[0] *= 2*CLHEP::pi*sin(part_rad[0]);
      part_acu.push_back((2.*CLHEP::pi/0.654)*part_beta[0]*part_rad[0]*part_rad[0]);


      for (unsigned ii = 1; ii<part_beta.size (); ii++) {
        part_rad[ii] *= CLHEP::pi/180.;
        part_beta[ii] *= 2*CLHEP::pi*sin(part_rad[ii]);
        part_acu.push_back(part_acu[ii-1]+(part_rad[ii]-part_rad[ii-1])*(part_beta[ii]+part_beta[ii-1])/2.);
      }
      for (unsigned ii = 0; ii<part_acu.size (); ii++) {
        part_acu[ii] /= part_acu[part_acu.size () -1];
      }
      fpetzoldinitialised = true;
    }
  }


  //for randomly generated r from [0,1] get cos(angle) according to tabulated Petzold phase function
  double KM3PetzoldScattering::GetPetzoldAngle(double r) {
    InitializePetzold();

    double angle;
    unsigned k = 0;
    while (k<part_acu.size () && r > part_acu[k]) k++;
    if (k == 0)
      angle = r * part_rad[0] / part_acu[0]; //needs improvement --> linear for now but better shape may be used
    else
      angle = part_rad[k-1] + (r-part_acu[k-1])*(part_rad[k] - part_rad[k-1])/(part_acu[k]-part_acu[k-1]);  //--> the same. it's linear in cumulative pdf, but it is better to have it linear in pdf
    return angle;
  }


  //for randomly generated r from [0,1] get cos(angle) according to Einstein-Smoluchowski scattering
  double KM3PetzoldScattering::GetESAngle(double r) {
    double b = 0.835;
    double p = 1./b;
    double q = (b+3.)*(r-0.5)/b;
    double d = q*q +p*p*p;
    double u1 = -q+sqrt(d);
    double u = pow(fabs(u1),1./3.);
    if (u1 < 0.) u = -u;
    double v1 = -q-sqrt(d);
    double v = pow(fabs(v1),1./3.);
    if (v1 < 0.) v = -v;

    double cosp = u+v;
    if (cosp > 1) cosp = 1;
    if (cosp < -1) cosp = -1;
    return acos(cosp);
  }

  double KM3PetzoldScattering::GetESVolScat(double lambda){
    return fESfactor*(1.7e-3*pow(550*nm/lambda,4.32)); //ES volume scattering function (m^-1)
  }

  double KM3PetzoldScattering::GetKopelevichVolScat(double lambda) {
    return fKopelFactor*(7.5e-3*(1.34*pow(550*nm/lambda,1.7)+0.312*pow(550*nm/lambda,0.3))); //KopelevichES scattering length (m^-1)
  }

  G4PhysicsOrderedFreeVector* KM3PetzoldScattering::PetzoldAttenuationLengthGenerator ()
  {
    G4PhysicsOrderedFreeVector* to_return = new G4PhysicsOrderedFreeVector;
    for (double e=1.*eV; e<8*eV; e+=0.025*eV)
      {
        to_return->InsertValues(e, GetTotalLength (CLHEP::c_light*CLHEP::h_Planck/e/1000*m)*m );
      }


    return to_return;
  }


// GetMeanFreePath()
// -----------------
//
  G4double KM3PetzoldScattering::GetMeanFreePath(const G4Track& aTrack,
                                         G4double ,
                                         G4ForceCondition* )
  {

    const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
    const G4Material* aMaterial = aTrack.GetMaterial();



    G4double thePhotonEnergy = aParticle->GetTotalEnergy();

    G4double AttenuationLength = DBL_MAX;

    //if (DefaultWater){
    if ((*thePhysicsTable)(aMaterial->GetIndex()))
      {

      G4bool isOutRange;

      AttenuationLength =
        (*thePhysicsTable)
	(aMaterial->GetIndex())
	->
        GetValue(thePhotonEnergy, isOutRange);

      }
    else {

      G4MaterialPropertiesTable* aMaterialPropertyTable =
        aMaterial->GetMaterialPropertiesTable();

      if(aMaterialPropertyTable){
        G4MaterialPropertyVector* AttenuationLengthVector =
          aMaterialPropertyTable->GetProperty("KM3PETZOLDSCATTERING");

        if(AttenuationLengthVector){
          AttenuationLength = AttenuationLengthVector ->
            Value(thePhotonEnergy);

        }
        else{
          ;
        }
      }
      else{
        ;
      }
    }


    return AttenuationLength;
  }


  double KM3PetzoldScattering::GetTotalLength(double lambda) { // in meters!
    //cout << fKopelFactor << ' ' <<  fESFactor << endl;

    return 1./(GetESVolScat(lambda)+GetKopelevichVolScat(lambda));
  }

  double KM3PetzoldScattering::GetAngle(double lambda) { // total scattering in meters, two random numbers from [0,1] are needed
    double r1=G4UniformRand(), r2=G4UniformRand();
    double pprob = GetKopelevichVolScat(lambda);
    double esprob = GetESVolScat(lambda);
    double eta = esprob/(pprob+esprob);
    double angle;

    if (r1 > eta)
      {
        angle = GetPetzoldAngle(r2);
        fLastProcess = KOPELEVITCH;
      }
    else
      {
        angle = GetESAngle(r2);
        fLastProcess = EINSTEINSMOLUCHOWSKI;
      }
    return angle;
  }

  //TEST: Drawable eta

  double KM3PetzoldScattering::GetEta(double *x, double *par){
    par=par;
    double lambda = x[0];
    double pprob = GetKopelevichVolScat(lambda);
    double esprob = GetESVolScat(lambda);
    double eta = esprob/(pprob+esprob);
    return eta;
  }

  void KM3PetzoldScattering::BuildThePhysicsTable()
  {
    //      Builds a table of scattering lengths for each material

    if (thePhysicsTable) return;

    const G4MaterialTable* theMaterialTable=
      G4Material::GetMaterialTable();
    G4int numOfMaterials = G4Material::GetNumberOfMaterials();

    // create a new physics table

    thePhysicsTable = new G4PhysicsTable(numOfMaterials);

    // loop for materials

    for (int i=0 ; i < numOfMaterials; i++)
      {
        G4PhysicsOrderedFreeVector* ScatteringLengths = NULL;

        G4MaterialPropertiesTable *aMaterialPropertiesTable =
          (*theMaterialTable)[i]->GetMaterialPropertiesTable();

        if(aMaterialPropertiesTable){

          G4MaterialPropertyVector* AttenuationLengthVector =
            aMaterialPropertiesTable->GetProperty("KM3PetzoldScattering");

          if (AttenuationLengthVector && AttenuationLengthVector->GetVectorLength () < 2)
            {
              DefaultWater = true;

              ScatteringLengths =
                PetzoldAttenuationLengthGenerator();

            }
        }

        thePhysicsTable->insertAt(i,ScatteringLengths);
      }
  }

  void KM3PetzoldScattering::BuildPhysicsTable(const G4ParticleDefinition&)
  {
    if (!thePhysicsTable) BuildThePhysicsTable();
  }


  G4VParticleChange* KM3PetzoldScattering::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
  {

    aParticleChange.Initialize(aTrack);

    const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();

    if (verboseLevel>0) {
      G4cout << "Predifined scattering Photon!" << G4endl;
      G4cout << "Old Momentum Direction: "
             << aParticle->GetMomentumDirection() << G4endl;
      G4cout << "Old Polarization: "
             << aParticle->GetPolarization() << G4endl;
    }
    double energy = aParticle->GetTotalEnergy ();

    double ScatAngle  = GetAngle (CLHEP::c_light*CLHEP::h_Planck/energy/1000*m);

    G4double PhiAngle = twopi*G4UniformRand();

    G4ThreeVector OldMomentumDirection, NewMomentumDirection;
    G4ThreeVector OldPolarization, NewPolarization;

    OldMomentumDirection = aParticle->GetMomentumDirection();
    OldMomentumDirection = OldMomentumDirection.unit();

    NewMomentumDirection.set(sin(ScatAngle)*cos(PhiAngle),
                             sin(ScatAngle)*sin(PhiAngle),
                             cos(ScatAngle));

    NewMomentumDirection.rotateUz(OldMomentumDirection);
    NewMomentumDirection = NewMomentumDirection.unit();

    OldPolarization = aParticle->GetPolarization();
    G4double constant = -1./NewMomentumDirection.dot(OldPolarization);

    NewPolarization = NewMomentumDirection + constant*OldPolarization;
    NewPolarization = NewPolarization.unit();

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

    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);


  }



  inline
  G4bool KM3PetzoldScattering::IsApplicable(const G4ParticleDefinition& aParticleType)
  {
    return ( &aParticleType == G4OpticalPhoton::OpticalPhoton() );
  }

  inline
  G4PhysicsTable* KM3PetzoldScattering::GetPhysicsTable() const
  {
    return thePhysicsTable;
  }

}
