// This file is part of the GenericLAND software library.
// $Id: BidoneInputDataReader.hh,v 1.1.1.1 2004/12/21 22:29:48 Bidonesim Exp $
//
// KM3InputDataReader.hh
// v.0 by Glenn Horton-Smith, Feb 12, 1999
// adapted for namespace support in Geant4.1, 4.2 -- 14-Jul-2000

#ifndef KM3InputDataReader_H
#define KM3InputDataReader_H 1

#include "globals.hh"

#include <string>
#include <map>

#include <G4Material.hh>
#include <G4SystemOfUnits.hh>


class G4Material;

namespace km3net{

  class KM3InputDataReader {
  private:
    enum class prop_option
    {
      none,
        wavelength,
        dy_dwavelength,
        energy,
        eV,
        shifter_density,
        constant
        };

  public:


    class MyTokenizer {
    private:
      std::istream *isptr;
    public:
      MyTokenizer(std::istream &is) { isptr= &is; nval=0.0;}

      enum { TT_EOF=-1, TT_STRING='a', TT_NUMBER='0' };

      int ttype;

      G4double nval;
      G4String sval;

      int nextToken(void);

      void dumpOn(std::ostream &os);
    };

    static int ReadMaterials(std::istream &is);

  private:
    static G4Material* CreateMaterial(std::string materialName,
                                      const std::map <std::string, double>& DensityPerComponents,
                                      double densityGperCM3, G4State state=kStateSolid);
  };
}
#endif // KM3InputDataReader_H
