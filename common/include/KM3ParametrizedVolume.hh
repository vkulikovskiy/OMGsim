#ifndef __KM3PARAMETRIZEDVOLUME_H_
#define __KM3PARAMETRIZEDVOLUME_H_

#include <G4SystemOfUnits.hh>
#include "G4Types.hh"

#include <map>
#include <string>


class G4LogicalVolume;
class G4VisAttributes;
class G4VPhysicalVolume;
class G4Step;
class G4MaterialPropertiesTable;
class TFormula;

#include <G4MaterialPropertyVector.hh>

#include <KM3ParametrizedVolume.hh>
#include <libconfig.h++>
namespace km3net
{
  const float inch = 2.54*cm;

  class KM3ParametrizedVolume
  {
  protected:
    static G4VisAttributes  *fgold;
    static G4VisAttributes  *fsilver;
    static G4VisAttributes  *fbrown;
    static G4VisAttributes  *forange;
    static G4VisAttributes  *fcyan;
    static G4VisAttributes  *ftransparent;

    std::string ConfigFilePrefix = " ";
  protected:

    struct VolumeParameters
    {
      std::string Prefix = "";
      //overlaps checking
      bool CheckOverlaps = false;

      bool IsOk = false;

      virtual void CalculateParameters ()
      {
      }

    };

  protected:
    class KM3Volume
    {
    public:
      KM3Volume (VolumeParameters* ptp = NULL) {fTheParameters=ptp;};
      virtual ~KM3Volume ();

      const VolumeParameters* GetParameters () const
      {
        return fTheParameters;
      }

      virtual G4LogicalVolume* Construct (const VolumeParameters* = NULL) =0;
      G4LogicalVolume* GetLogicalVolume ()
      {
        return fTheLogicalVolume;
      }

      VolumeParameters        *fTheParameters;
      G4LogicalVolume         *fTheLogicalVolume=NULL;
    };


  public:
    virtual void SetParametersListFromMac (std::string Path);
    virtual void SetDataFilesSources (std::string Path);
    virtual void LoadDirectory (std::string);
  protected:
    void ReadAStringParameterFromMac (std::string SPar);
    virtual void ReadConfiguration (std::string) =0;
    bool ReadBaseParameters (libconfig::Config*,  VolumeParameters*);
    bool ReadBaseParameters (libconfig::Setting*,  VolumeParameters*);

  public:
    G4LogicalVolume* GetLogicalVolume (std::string);
    const KM3Volume* GetVolume (std::string);

    long GetIntegerParameter  (const libconfig::Setting& STG, float x=0, float y=0, float z=0);
    double GetFloatingParameter (const libconfig::Setting& STG, float x=0, float y=0, float z=0);
    std::string GetStringParameter   (const libconfig::Setting& STG);
    std::vector<double> GetScalarArrayParameter  (const libconfig::Setting& STG) const;

    long GetIntegerParameter  (std::string NAME);
    double GetFloatingParameter (std::string NAME);
    std::string GetStringParameter   (std::string NAME);
    std::string GetFormulaParameter  (std::string NAME);


  protected:
    KM3ParametrizedVolume ();
    ~KM3ParametrizedVolume ();

    virtual KM3Volume* ConstructTheVolume (VolumeParameters*) =0;

    std::string isAFormula (const libconfig::Setting& STG);
    std::string isAFormula (std::string& STR);
    void FormulaLooping ();
    std::vector<double> ReadThisArray (const libconfig::Setting& STG) const;
    std::string ReadThisFormula (const libconfig::Setting& STG);
    std::string ReadThisFormula (std::string& STR);

  protected:

    std::map <std::string, VolumeParameters*>      fParametersDict;
    std::map <std::string, KM3Volume*>             fVolumeDict;
    std::map <const G4LogicalVolume*, std::string> fLogicalVolDict;


    //private:
    static std::map <std::string, long>         fIntegerParameters;
    static std::map <std::string, double>       fFloatingParameters;
    static std::map <std::string, std::string>  fStringParameters;
    static std::map <std::string, std::string>  fFormulaParameters;

    std::map <std::string, std::vector<double> >         fScalarArrayParameters;
    std::map <std::string, std::vector<std::string> >    fStringArrayParameters;
    /*The parameters' names should start with __*/
    void ReadParameterList (std::string);

    TFormula* fFormula;

  };
}

#endif
