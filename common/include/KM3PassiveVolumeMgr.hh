
#ifndef __KM3PASSIVEVOLUMEMGR_H_
#define __KM3PASSIVEVOLUMEMGR_H_


#include <G4SystemOfUnits.hh>
#include "G4Types.hh"

#include <map>
#include <utility>
#include <string>

#include <KM3ParametrizedVolume.hh>
namespace libconfig
{
  class Config;
  class Setting;
}
class G4LogicalVolume;
class G4VisAttributes;
class G4VPhysicalVolume;
class G4Step;
class G4MaterialPropertiesTable;
class G4VSolid;

#include <G4MaterialPropertyVector.hh>

#include <KM3PhotoTubeMgr.hh>

namespace km3net
{
  /**
   * @brief Manage the automatic volume generator, configured in common/data/KM3Vol*.dat, see VolumeParameters
   *
   *
   *
   *
   */
  class KM3PassiveVolumeMgr: public KM3ParametrizedVolume
  {
  public:
    enum EnumVolumeKind
      {
        Ellipsoid,
        Box,
        Con,
        Tub,
        Sphere,
        Torus,
        Hexagone,
        Union,
        Subtraction,
        Intersection
      };
  private:
    /**
     *
     * @brief  Contains the Volume parameter, managed by KM3PassiveVolumeMgr and configured in common/data/KM3Vol*.dat
     *
     * The parameters for the volume properties (in addition of the ones inherited from KM3ParametrizedVolume) are

     *     std::string RefVolume;     : for the boolean kind, the volume on which the operation will be applied (gives the 0)
     *     std::string OperandVolume; : for the boolean kind
     *     double RotTheta;               : for the boolean kind
     *     double RotPhi;               : for the boolean kind
     *
     *     double X; : for the Box and boolean kinds
     *     double Y; : for the Box and boolean kinds
     *     double Z; : for the Box, Tubes, Cone and boolean kinds
     *
     *     double RadiusMin;         : For the Sphere, Tubes, Hexagone (opt) and cone kinds
     *     double RadiusMax;         : For the Sphere, Torus, Tubes, Ellipsoid and cone kinds
     *     double PhiMin   ;         : For the Sphere, Torus, Tubes and cone kind
     *     double PhiMax   ;         : For the Sphere, Torus, Tubes and cone kind
     *     double ThetaMin ;         : For the Sphere kinds
     *     double ThetaMax ;         : For the Sphere kinds
     *
     *     double SecondRadiusMin; : For the Cons and Torus kind
     *     double SecondRadiusMax; : For the Cons and Ellipsoid kinds
     *
     *     double ThirdRadiusMax;    : For the Ellipsoid kind
     *     double ZCutMin = 0;    : For the Ellipsoid kind
     *     double ZCutMax = 0;    : For the Ellipsoid kind
     *
     */
    struct PassiveVolumeParameters: public KM3ParametrizedVolume::VolumeParameters
    {
      std::string Name;         /**< The name used for invokation (eg for the GPS) */
      std::string VolumeKind;   /**< Chose between Ellipsoid,Box,Con,Tub,Sphere,Torus,Union,Subtraction and Intersection for the volume*/

      //predefined values are optional

      std::string RefVolume;     /**< for the boolean kind, the volume on which the operation will be applied (gives the 0)*/
      std::string OperandVolume; /**< for the boolean kind*/
      double RotTheta;               /**< for the boolean kind*/
      double RotPhi;               /**< for the boolean kind*/

      double X; /**< for the Box and boolean kinds */
      double Y; /**< for the Box and boolean kinds */
      double Z; /**< for the Box, Tubes, Cone, Hexagone and boolean kinds */

      double RadiusMin;         /**< For the Sphere, Tubes and cone kinds */
      double RadiusMax;         /**< For the Sphere, Torus, Hexagone, Tubes, Ellipsoid and cone kinds */
      double PhiMin   ;         /**< For the Sphere, Torus, Tubes and cone kind */
      double PhiMax   ;         /**< For the Sphere, Torus, Tubes, Hexagone (opt) and cone kind */
      double ThetaMin ;         /**< For the Sphere kinds */
      double ThetaMax ;         /**< For the Sphere kinds */

      double SecondRadiusMin; /**< For the Cons and Torus kind */
      double SecondRadiusMax; /**< For the Cons and Ellipsoid kinds */

      double ThirdRadiusMax;    /**< For the Ellipsoid kind */
      double ZCutMin = 0;    /**< For the Ellipsoid kind */
      double ZCutMax = 0;    /**< For the Ellipsoid kind */

      //Material
      std::string Material = "G4_AIR";
      std::string SkinMaterial = "";

      void CalculateParameters ()
      {
      }

    };

  private:
    class KM3PassiveVolume: public KM3ParametrizedVolume::KM3Volume
    {
    public:
      KM3PassiveVolume (PassiveVolumeParameters* PTP);
      ~KM3PassiveVolume ();


      G4LogicalVolume* Construct (const VolumeParameters* = NULL);

    private:
      G4LogicalVolume* ConstructEllipsoidVolume ();
      G4LogicalVolume* ConstructBoxVolume ();
      G4LogicalVolume* ConstructConVolume ();
      G4LogicalVolume* ConstructTubVolume ();
      G4LogicalVolume* ConstructSphereVolume ();
      G4LogicalVolume* ConstructTorusVolume ();
      G4LogicalVolume* ConstructHexagoneVolume ();
      G4LogicalVolume* ConstructBooleanVolume ();

      G4LogicalVolume* MakeItLogical (G4VSolid* Solid);
    };


  public:
    static KM3PassiveVolumeMgr* GetIt ();

    void ReadConfiguration (std::string);
    void ReadCurrentConfiguration (libconfig::Setting* stg);

  private:
    void ReadEllipsoidConfiguration(libconfig::Setting* stg);
    void ReadBoxConfiguration(libconfig::Setting* stg);
    void ReadConConfiguration(libconfig::Setting* stg);
    void ReadTubConfiguration(libconfig::Setting* stg);
    void ReadSphereConfiguration(libconfig::Setting* stg);
    void ReadTorusConfiguration(libconfig::Setting* stg);
    void ReadHexagoneConfiguration(libconfig::Setting* stg);
    void ReadBooleanConfiguration(libconfig::Setting* stg);


  private:
    static std::string fname;
    static KM3ParametrizedVolume* fInstance;
    KM3Volume* ConstructTheVolume (VolumeParameters*);
    KM3PassiveVolumeMgr ();
    ~KM3PassiveVolumeMgr ();
  };
}
#endif
