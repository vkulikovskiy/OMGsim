#ifndef __KM3PHOTOTUBEMGR_H_
#define __KM3PHOTOTUBEMGR_H_

#include <G4SystemOfUnits.hh>
#include "G4Types.hh"

#include <map>
#include <string>

#include <KM3ParametrizedVolume.hh>

class G4LogicalVolume;
class G4VisAttributes;
class G4VPhysicalVolume;
class G4Step;
class G4MaterialPropertiesTable;

#include <G4MaterialPropertyVector.hh>

namespace km3net
{
  class KM3PhotoTubeMgr : public KM3ParametrizedVolume
  {
  public:
    static const std::string fPhotoName;
    static const std::string fPhotoRegionName;
    static const std::string fAbsorberName;
    static const std::string fBentoName;
    static const std::string fPMTContainerName;


    struct PhotoTubeParameters: public KM3ParametrizedVolume::VolumeParameters
    {
      //PMT input parameters (ANTARES pmt values)
      double RadiusSphere   = 136.7 * mm;
      double BulbRadiusMin  = 220. / 2 * mm;
      double BulbRadiusMax  = 253. / 2 * mm;
      double BulbHeight     = 80. * mm;
      double ConeHeight     = 40. * mm;
      double TubeRadius     = 84.5 / 2 * mm;
      double GlassThickness = 1 * mm;
      double TotalPMTHeight = 245 * mm;
      double DinodsPos      = -1.5 * cm;
      double DinodsRadius   = 3.25 * cm;
      double DinodsHeight   = 3.5 * cm;

      //to be calculated
      double PhotoAngle;
      double PhotoAngleIn;
      double PhotoHeight;
      double Bulbr;

      double BulbRIn;
      double BulbrIn;

      //The materials to get, "classic" big PM by default
      std::string GlassMaterial        = "PMTGlass";
      std::string PhotoCathodeMaterial = "photocathode10inches";
      std::string InternalReflectorMaterial  = "MirrorAluminium";
      std::string DinodsSurface              = "DinodsSurface"            ;
      std::string CapDinodsSurface           = "CapDinodsSurface"         ;
      std::string GridDinodsSurface          = "GridDinodsSurface"        ;
      std::string PMTSurface                 = ""                         ;
      void CalculateParameters ()
      {
        PhotoAngle     = asin(BulbRadiusMin/(RadiusSphere));
        PhotoAngleIn   = asin((BulbRadiusMin-GlassThickness/sin (PhotoAngle))/(RadiusSphere-GlassThickness));
        PhotoHeight    = RadiusSphere - cos (PhotoAngle) * RadiusSphere;
        Bulbr          = sqrt (BulbHeight*BulbHeight/4/(1-BulbRadiusMin*BulbRadiusMin/BulbRadiusMax/BulbRadiusMax));

        BulbRIn=BulbRadiusMax-GlassThickness;
        BulbrIn=sqrt (BulbHeight*BulbHeight/4/(1-pow(BulbRadiusMin-GlassThickness/sin(PhotoAngle),2)/BulbRIn/BulbRIn));
      }

    };

  private:
    class KM3PhotoTube: public KM3ParametrizedVolume::KM3Volume
    {
    public:
      KM3PhotoTube (PhotoTubeParameters* PTP);
      ~KM3PhotoTube ();

      G4LogicalVolume* Construct (const VolumeParameters* = NULL);

      unsigned char CheckDetection (G4Step*);
      unsigned char CalculateDetection (const G4Step*);
      bool LoadMaterialProps (G4MaterialPropertiesTable*);

      G4MaterialPropertyVector* _efficiency_photocathode;
      G4MaterialPropertyVector* _angular_efficiency;
      G4MaterialPropertyVector* _angular_efficiency_err;

    };


  public:
    static KM3PhotoTubeMgr* GetIt ();
    void ReadConfiguration (std::string);

    bool          IsInAPhotochatode ();
    unsigned char IsDetected (G4Step* STEP);

  private:
    KM3PhotoTubeMgr ();
    ~KM3PhotoTubeMgr ();

  private:
    static KM3ParametrizedVolume* fInstance;
    KM3Volume* ConstructTheVolume (VolumeParameters*);
    std::string fCurrentPMTName;
  };
}

#endif
