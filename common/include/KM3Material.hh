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
#ifndef KM3Material_HH
#define KM3Material_HH
////////////////////////////////////////////////////////////////////////////////
#include "G4Material.hh"
#include "G4PhysicalConstants.hh"
#include <vector>
#include <map>

class G4OpticalSurface;

namespace km3net
{

  class KM3MaterialMessenger;

  ////////////////////////////////////////////////////////////////////////////////
  //
  class KM3Material {
    friend class  KM3MaterialMessenger;

  public:
    static KM3Material* GetIt ();

  private:

    KM3Material();
    ~KM3Material();
    static KM3Material *fInstance;

  public:

    void AddMaterial(G4String, G4String, G4double, G4String, G4double tem =
                     STP_Temperature, G4double pres = STP_Pressure);
    void AddMaterial(G4Material* aMaterial) {
      Material.push_back(aMaterial);
    }
    ;

    void RemoveScattering ();
    void SetDataFilesSources (std::string);
    G4Material* GetMaterial(G4int i) {
      return Material[i];
    }
    ;
    const G4Material* GetMaterial(G4int i) const {
      return Material[i];
    }
    ;
    static G4Material* GetMaterial(G4String name);
    G4OpticalSurface* GetOpticalSurface(G4String name);
    G4int GetMaterialIndex(G4String);
    G4int GetNbOfMaterial() {
      return Material.size();
    }
    ;
    void DeleteMaterial(G4int);
    void DeleteMaterial(G4String);

    void ListMaterial();
    void DefineMaterials();

  private:

    std::vector<G4Material*> Material;
    std::vector<G4Element*> Element;
    std::vector<G4Isotope*> Isotope;

    std::map<G4String,G4OpticalSurface*> fOpticalSurfaceDict;

  private:
    static const G4String ELU[110];
    static const G4String ELL[110];
    static const G4String EUU[110];
    static const G4double A[110];


    KM3MaterialMessenger* fMaterialMessenger;

  };
}
////////////////////////////////////////////////////////////////////////////////
#endif
