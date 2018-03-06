#ifndef __KM3OPTICALMODULEMGR_H_
#define __KM3OPTICALMODULEMGR_H_


#include <G4SystemOfUnits.hh>
#include "G4Types.hh"

#include <map>
#include <utility>
#include <string>


#include <KM3ParametrizedVolume.hh>

class G4LogicalVolume;
class G4VisAttributes;
class G4VPhysicalVolume;
class G4Step;
class G4MaterialPropertiesTable;


#include <G4MaterialPropertyVector.hh>

#include <KM3PhotoTubeMgr.hh>

namespace km3net
{
  class KM3OpticalModuleMgr: public KM3ParametrizedVolume
  {
  public:
    static const std::string fAbsorberName;
    static const std::string fBentoName;
    static const std::string fPMTContainerName;

  private:
    struct OpticalModuleParameters: public KM3ParametrizedVolume::VolumeParameters
    {
      const KM3PhotoTubeMgr::PhotoTubeParameters *PMTprt = NULL;
      std::string PhotoTubeName;
      //Choice of the type, dom or om
      std::string DetectionUnitType="OM";
      //OM input parameters (ANTARES values)
      double SphereRmax       = 17./2 * inch;
      double GlassThickness   = 1.5 * cm;
      double AbsorberOvertake = 0;
      double GelOvertake      = 14.5* mm; //with more, goes to bigger angles
      double DistSpherePMT     = 0.5 * cm;
      double DistAbsorberGlass =0.25*mm;
      std::string AbsorberSurface = "AbsorberPVC";
      std::string BentoSurface = "";
      //DOM specific values
      bool   AreThereReflectors = true;
      double ReflectorAngle     = 30*deg;
      double ReflectorRadiusMin = 40.25*mm;
      double ReflectorDistGlass = 1*mm;
      double ReflectorThickness = 12*mm;
      double MushroomAngularSize=0.7; //the mushroom will always be at -180 deg
      std::string MushroomMaterial="RoughAluminium";
      std::vector<std::pair<float,float>> PMTsAngles;

      //to be calculated
      double GelAngleOvertake;
      double AbsorberOvertakeAngle;
      double SphereRmin;

      //The materials to get, "classic" OM with big PM by default
      std::string GlassMaterial = "AntaresGlass";
      std::string GelMaterial   = "WackerSilGel612_A100B67";
      std::string AbsMaterial   = "AbsorberPVC";


      void CalculateParameters ()
      {
        SphereRmin       = SphereRmax - GlassThickness;;
        GelAngleOvertake = GelOvertake / PMTprt->RadiusSphere;
        AbsorberOvertakeAngle = AbsorberOvertake / SphereRmin;
      }

    };

  private:
    class KM3OpticalModule: public KM3ParametrizedVolume::KM3Volume
    {
    public:
      KM3OpticalModule (OpticalModuleParameters* PTP);
      ~KM3OpticalModule ();

      G4LogicalVolume* GetDOM ();
      G4LogicalVolume* GetOM ();

      G4LogicalVolume* Construct (const VolumeParameters* = NULL);

    };


  public:
    static KM3OpticalModuleMgr* GetIt ();

    void ReadConfiguration (std::string);


  private:
    static KM3ParametrizedVolume* fInstance;
    KM3Volume* ConstructTheVolume (VolumeParameters*);
    KM3OpticalModuleMgr ();
    ~KM3OpticalModuleMgr ();

  };
}
#endif
