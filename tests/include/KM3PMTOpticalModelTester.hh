#ifndef  __KM3PMTOpticalModelTester_hh__
#define  __KM3PMTOpticalModelTester_hh__

#include <KM3PMTOpticalModel.hh>
#include <KM3Material.hh>

#include <string>
#include <map>

#include <G4ThreeVector.hh>

class TObject;

using namespace std;
namespace km3net
{
  class KM3PMTOpticalModelTester
  {
  public:
    KM3PMTOpticalModelTester ();
    ~KM3PMTOpticalModelTester () {};

    bool LoadMaterial (std::string PhMaterial, std::string InMaterial, std::string OutMaterial);
    bool TestElectronProduction (G4ThreeVector &dir, G4ThreeVector &pol, G4ThreeVector &norm, G4ThreeVector &insideDirection);
    float CalcElectronProductionProb(G4ThreeVector &dir, G4ThreeVector &pol, G4ThreeVector &norm, G4ThreeVector &insideDirection);
    void SaveTreesAndHistograms (std::string FileName);
    float TestCalculateNormalCoefficients (G4ThreeVector &pos, float _wavelength);
    float TestCalculateCoefficients (G4ThreeVector &pos, G4ThreeVector &dir, G4ThreeVector &pol, G4ThreeVector &norm, float wavelength);
    float GetQuantumEfficiency (G4ThreeVector &pos, G4ThreeVector &dir, G4ThreeVector &pol, G4ThreeVector &norm, G4ThreeVector &insideDirection, float wavelength);

    void QuantumEfficiencyMeasurement ();
		float QECheck (double wavelength, double theta = 0);
    void PhotoCathodeEfficiencyReader ();

    double TestElectronEscapingProbability (const G4ThreeVector &dir, const G4ThreeVector &pol, const G4ThreeVector &norm) const;
    bool LastElectronIsEscaped ();
    float GetEfficiency (float wavelength);
    float GetGlassAbsorption (float wavelength, float thickness);
    G4MaterialPropertyVector GetKIndex () {return *(fOpticalModel._kindex_photocathode);}

    void SetKIndex (const G4MaterialPropertyVector& vector) {
      fOpticalModel._kindex_photocathode = new G4MaterialPropertyVector (vector);
    }

    void SetCathodeThickness (const G4MaterialPropertyVector& vector) {
      fOpticalModel._thickness_photocathode = new G4MaterialPropertyVector (vector);
    }

    void SetElectronProductionLevel (int level) {fOpticalModel._electronProductionLevel = level;}
    void SetElectronAbsorptionLength (float length) {fOpticalModel.fElectronAttenuationLength=length;}


    TObject* GetTObject (string name) { return fObjectCol[name];}
    void PutTObject (string name, TObject* object) { fObjectCol[name] = object; return;}
    bool InTObject (string name) {if (fObjectCol.find(name) != fObjectCol.end()) return true; return false;}
    void GetTotalProbabilities (double& Abs, double& Trans, double& Refl) {Abs=fOpticalModel.fA_t; Trans=fOpticalModel.fT_t; Refl=fOpticalModel.fR_t;}
    float fGlassThickness = 1*mm;
    int fElPrecision = 10000;

  private:
    G4MaterialPropertyVector *fGlassAbsLength = NULL;

    KM3PMTOpticalModel fOpticalModel;
    map<string,TObject*> fObjectCol;
  };
}

#endif

