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
////////////////////////////////////////////////////////////////////////////////
//
#include "KM3Material.hh"
#include "KM3MaterialData.hh"
#include "KM3InputDataReader.hh"
#include <KM3MaterialMessenger.hh>

#include "globals.hh"
#include "G4UnitsTable.hh"
#include "G4NistManager.hh"
#include "G4ios.hh"
#include <G4OpticalSurface.hh>

#include <vector>
#include <iomanip>
#include <iostream>
#include <sys/types.h>
#include <dirent.h>

#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>


using namespace std;
using namespace CLHEP;

////////////////////////////////////////////////////////////////////////////////
//
namespace km3net
{


  void ReadMaterialProperties (string SFILE)
  {
    std::ifstream ifs(SFILE.data ());

    if (ifs.fail()) {
      G4cerr << "Error, material properties file could not be opened.\n";
      if (getenv("GLG4DATA") == NULL)
        G4cerr
          << "GLG4DATA environment variable is not set, so I was looking for data/photocathode.dat from the current directory."
          << G4endl;
      else
        G4cerr << "I was looking for photocathode.dat in the GLG4DATA directory, " << getenv("GLG4DATA") << G4endl;
      G4Exception("Error, material properties file could not be opened.\n",
                  "", FatalException, "");
    }

    // now read materials, keeping error count
    KM3InputDataReader::ReadMaterials(ifs);
    // close file
    ifs.close();
  }

  void ScanDirectory (string DIRNAME)
  {
    struct dirent *de=NULL;
    DIR *d=NULL;

    d=opendir(DIRNAME.data ());
    if(d == NULL)
      {
        cerr << "Couldn't open directory" << DIRNAME << endl;
        return;
      }

    // Loop while not NULL
    while( (de = readdir(d)) )
      {
        string filename=de->d_name;
        if (de->d_type == 8 && filename.substr (0,6) == "KM3Mat" && filename.rfind (".dat") == filename.size () - 4)
          {
            clog << "KM3Material: Loading file " << DIRNAME  << '/' << de->d_name<< endl;
            ReadMaterialProperties (DIRNAME+"/"+filename);
          }
      }
    closedir(d);
  }

  KM3Material* KM3Material::fInstance=NULL;

  KM3Material*  KM3Material::GetIt ()
  {
    if (!fInstance)
      fInstance=new KM3Material ();
    return fInstance;
  }

  G4Material* KM3Material::GetMaterial(G4String name) {
    G4Material* toreturn = G4Material::GetMaterial(name);
    if (!toreturn)
      {
        clog << "Looking in G4NistManager for " << name <<  endl;
        toreturn = G4NistManager::Instance()->FindOrBuildMaterial(name);
      }
    return toreturn;
  }

  G4OpticalSurface* KM3Material::GetOpticalSurface (G4String name)
  {
    if (fOpticalSurfaceDict.count (name))
      return fOpticalSurfaceDict [name];
    cout << "creating the surface" << endl;
    G4OpticalSurface* the_surface =new G4OpticalSurface(name + "_surface");
    G4MaterialPropertiesTable* thePropTable = G4Material::GetMaterial(name)->GetMaterialPropertiesTable();
    //thePropTable->DumpTable ();

    the_surface->SetMaterialPropertiesTable(thePropTable);
    //converting the constparameters to surface properties
    //    if (thePropTable->GetProperty ("REFLECTIVITY") != NULL && thePropTable->GetPropertiesCMap ()->size() )
      {

        if (thePropTable->ConstPropertyExists ("dielectric_metal"))
          the_surface->SetType (dielectric_metal);
        else
          the_surface->SetType(dielectric_dielectric);


        cout << "specifiying the reflect" << endl;
        if (thePropTable->ConstPropertyExists ("polishedfrontpainted")) the_surface->SetFinish (polishedfrontpainted);
        else if (thePropTable->ConstPropertyExists ("polished")) the_surface->SetFinish (polished);
        else if (thePropTable->ConstPropertyExists ("polishedbackpainted")) the_surface->SetFinish (polishedbackpainted);
        else if (thePropTable->ConstPropertyExists ("ground")) the_surface->SetFinish (ground);
        else if (thePropTable->ConstPropertyExists ("groundfrontpainted")) the_surface->SetFinish (groundfrontpainted);
        else if (thePropTable->ConstPropertyExists ("groundbackpainted")) the_surface->SetFinish (groundbackpainted);
        else if (thePropTable->ConstPropertyExists ("polishedlumirrorair")) the_surface->SetFinish (polishedlumirrorair);
        else if (thePropTable->ConstPropertyExists ("polishedlumirrorglue")) the_surface->SetFinish (polishedlumirrorglue);
        else if (thePropTable->ConstPropertyExists ("polishedair")) the_surface->SetFinish (polishedair);
        else if (thePropTable->ConstPropertyExists ("polishedteflonair")) the_surface->SetFinish (polishedteflonair);
        else if (thePropTable->ConstPropertyExists ("polishedtioair")) the_surface->SetFinish (polishedtioair);
        else if (thePropTable->ConstPropertyExists ("polishedtyvekair")) the_surface->SetFinish (polishedtyvekair);
        else if (thePropTable->ConstPropertyExists ("polishedvm2000air")) the_surface->SetFinish (polishedvm2000air);
        else if (thePropTable->ConstPropertyExists ("polishedvm2000glue")) the_surface->SetFinish (polishedvm2000glue);
        else if (thePropTable->ConstPropertyExists ("etchedlumirrorair")) the_surface->SetFinish (etchedlumirrorair);
        else if (thePropTable->ConstPropertyExists ("etchedlumirrorglue")) the_surface->SetFinish (etchedlumirrorglue);
        else if (thePropTable->ConstPropertyExists ("etchedair")) the_surface->SetFinish (etchedair);
        else if (thePropTable->ConstPropertyExists ("etchedteflonair")) the_surface->SetFinish (etchedteflonair);
        else if (thePropTable->ConstPropertyExists ("etchedtioair")) the_surface->SetFinish (etchedtioair);
        else if (thePropTable->ConstPropertyExists ("etchedtyvekair")) the_surface->SetFinish (etchedtyvekair);
        else if (thePropTable->ConstPropertyExists ("etchedvm2000air")) the_surface->SetFinish (etchedvm2000air);
        else if (thePropTable->ConstPropertyExists ("etchedvm2000glue")) the_surface->SetFinish (etchedvm2000glue);
        else if (thePropTable->ConstPropertyExists ("groundlumirrorair")) the_surface->SetFinish (groundlumirrorair);
        else if (thePropTable->ConstPropertyExists ("groundlumirrorglue")) the_surface->SetFinish (groundlumirrorglue);
        else if (thePropTable->ConstPropertyExists ("groundair")) the_surface->SetFinish (groundair);
        else if (thePropTable->ConstPropertyExists ("groundteflonair")) the_surface->SetFinish (groundteflonair);
        else if (thePropTable->ConstPropertyExists ("groundtioair")) the_surface->SetFinish (groundtioair);
        else if (thePropTable->ConstPropertyExists ("groundtyvekair")) the_surface->SetFinish (groundtyvekair);
        else if (thePropTable->ConstPropertyExists ("groundvm2000air")) the_surface->SetFinish (groundvm2000air);
        else if (thePropTable->ConstPropertyExists ("groundvm2000glue")) the_surface->SetFinish (groundvm2000glue);


        if (thePropTable->ConstPropertyExists ("unified")) the_surface->SetModel (unified);
        else if (thePropTable->ConstPropertyExists ("LUT")) the_surface->SetModel (LUT);
        else if (thePropTable->ConstPropertyExists ("dichroic")) the_surface->SetModel (dichroic);
        else if (thePropTable->ConstPropertyExists ("glisur")) the_surface->SetModel (glisur);

        if (thePropTable->ConstPropertyExists ("polish"))
          the_surface->SetPolish (thePropTable->GetConstProperty ("polish"));

        cout << the_surface->GetFinish () << " " << the_surface->GetPolish () <<' ' << polished <<" " << the_surface->GetModel () <<  endl;

        if (thePropTable->ConstPropertyExists ("sigmaalpha"))
          the_surface->SetSigmaAlpha (thePropTable->GetConstProperty ("sigmaalpha"));



      }

    fOpticalSurfaceDict [name] = the_surface;
    //the_surface->GetMaterialPropertiesTable () -> DumpTable ();
    return the_surface;
  }

  void KM3Material::SetDataFilesSources (string Path){
    size_t doubledotPos=Path.find (':');
    string subPath=Path;
    while (doubledotPos != string::npos) {
      subPath=Path.substr (0,doubledotPos);
      ScanDirectory (subPath);
      Path.erase (0,doubledotPos+1);
      doubledotPos=Path.find(':');
    }
    if (Path.size ())
      ScanDirectory (subPath);
  }

 KM3Material::KM3Material() {
    Material.clear();
    Element.clear();
    Isotope.clear();

    ScanDirectory ((string(PROJECT_SOURCE_DIR)+"/common/data/"));

    DefineMaterials();

    fMaterialMessenger=new KM3MaterialMessenger (this);

  }
  ////////////////////////////////////////////////////////////////////////////////
  //
  KM3Material::~KM3Material() {
    delete fMaterialMessenger;
  }

  void KM3Material::RemoveScattering ()
  {
  }

  ////////////////////////////////////////////////////////////////////////////////
  //
  void KM3Material::AddMaterial(G4String name, G4String formula, G4double density,
                                G4String state, G4double tem, G4double pres) {
    G4int isotope, Z;
    size_t i;
    for (i = 0; i < Material.size(); i++) {
      if (Material[i]->GetName() == name) {
        G4cerr << " AddMaterial : material " << name << " already exists."
               << G4endl;
        G4cerr <<"--> Command rejected." <<G4endl;
        return;
      }
    }

    char *tokenPtr1 = NULL;
    char *sname = NULL;
    G4String s, s1("0123456789");
    G4String element, isotopename;
    G4int ncomponents, natoms;
    G4double fatoms = 0.;
    size_t ls, id = 0, ll, lr;
    ncomponents = 0;

    sname = new char[strlen(formula) + 1];
    strcpy(sname, formula);
    tokenPtr1 = strtok(sname, "-");

    while (tokenPtr1 != NULL) {
      ncomponents++;
      tokenPtr1 = strtok(NULL, "-");
    }
    delete[] sname;

    G4Material* aMaterial = 0;
    //  G4cout << name <<" "<< formula << " " << density/(g/cm3) << " " << tem <<" " <<pres << G4endl;

    if (state == "") {
      aMaterial = new G4Material(name, density, ncomponents);
    } else if (state == "solid" && tem > 0.) {
      aMaterial = new G4Material(name, density, ncomponents, kStateSolid,
                                 tem);
    } else if (state == "gas" && pres > 0.) {
      aMaterial = new G4Material(name, density, ncomponents, kStateGas, tem,
                                 pres);
    }
    if (aMaterial == 0) {
      G4cerr << " AddMaterial : Name " << name << "." << G4endl;
      G4cerr <<"--> Command failed." <<G4endl;
      return;
    }

    sname = new char[strlen(formula) + 1];
    strcpy(sname, formula);
    tokenPtr1 = strtok(sname, "-");

    while (tokenPtr1 != NULL) {
      isotope = 0;
      //      G4cout << tokenPtr1 << G4endl;
      s = G4String(tokenPtr1);
      ls = s.length();
      ll = s.find("(");
      lr = s.find(")");
      if (ll == lr) {
        id = s.find_first_of(s1);
        element = s.substr(0, id);

        if (element.length() == 1)
          element.insert(0, " ");
        for (i = 0; i < 110; i++) {
          if (element == ELU[i])
            break;
        }
        if (i == 110) {
          for (i = 0; i < 110; i++) {
            if (element == ELL[i])
              break;
          }
          if (i == 110) {
            for (i = 0; i < 110; i++) {
              if (element == EUU[i])
                break;
            }
          }
        }

        if (i == 110) {
          G4cerr << "AddMaterial : Invalid element in material formula."
                 << element << G4endl;
          G4cerr <<"--> Command rejected." <<G4endl;
          //        delete aMaterial;
          //	Material[NbMat] = NULL;
          return;
        }

        Z = i + 1;
        element = ELU[i];
        if (id == std::string::npos) {
          natoms = 1;
        } else {
          natoms = atoi((s.substr(id, ls - id)).c_str());
        }
        if (natoms < 1)
          fatoms = atof((s.substr(id, ls - id)).c_str());
        //	G4cout << "   Elements = " << element << G4endl;
        //G4cout << "   Nb of atoms = " << natoms << G4endl;
      } else {
        element = s.substr(0, ll);
        isotope = atoi((s.substr(ll + 1, lr - ll)).c_str());
        if (element.length() == 1)
          element.insert(0, " ");
        for (i = 0; i < 110; i++) {
          if (element == ELU[i])
            break;
        }
        if (i == 110) {
          for (i = 0; i < 110; i++) {
            if (element == ELL[i])
              break;
          }
          if (i == 110) {
            for (i = 0; i < 110; i++) {
              if (element == EUU[i])
                break;
            }
          }
        }
        if (i == 110) {
          G4cerr << "AddMaterial : Invalid element in material formula."
                 << element << G4endl;
          G4cerr <<"--> Command rejected." <<G4endl;
          //        delete aMaterial;
          //	Material[NbMat] = NULL;
          return;
        }

        Z = i + 1;
        element = ELU[i];
        isotopename = element + s.substr(ll + 1, lr - ll - 1);
        if (lr == std::string::npos) {
          natoms = 1;
        } else {
          natoms = atoi((s.substr(lr + 1, ls - lr)).c_str());
        }
        if (natoms < 1)
          fatoms = atof((s.substr(id, ls - id)).c_str());
        if (fatoms == 0.)
          natoms = 1;
        //
        //	G4cout << "   Elements = " << element << G4endl;
        //   G4cout << "   Isotope Nb = " << isotope << G4endl;
        //	G4cout << "   Nb of atoms = " << natoms << G4endl;
      }
      if (isotope != 0) {
        if (G4Isotope::GetIsotope(isotopename) == NULL) {
          //        G4Isotope* aIsotope = new G4Isotope(isotopename, Z, isotope, A[Z-1]*g/mole);
          G4Isotope* aIsotope = new G4Isotope(isotopename, Z, isotope,
                                              isotope * g / mole);
          G4Element* aElement = new G4Element(isotopename, element, 1);
          aElement->AddIsotope(aIsotope, 100. * perCent);
          Isotope.push_back(aIsotope);
          if (natoms > 0) {
            aMaterial->AddElement(aElement, natoms);
          } else {
            aMaterial->AddElement(aElement, fatoms);
          }
          Element.push_back(aElement);
        } else {
          if (natoms > 0) {
            aMaterial->AddElement(G4Element::GetElement(isotopename),
                                  natoms);
          } else {
            aMaterial->AddElement(G4Element::GetElement(isotopename),
                                  fatoms);
          }
        }
      } else {
        if (G4Element::GetElement(element) == NULL) {
          G4Element* aElement = new G4Element(element, element, Z,
                                              A[Z - 1] * g / mole);
          if (natoms > 0) {
            aMaterial->AddElement(aElement, natoms);
          } else {
            aMaterial->AddElement(aElement, fatoms);
          }
          Element.push_back(aElement);
        } else {
          if (natoms > 0) {
            aMaterial->AddElement(G4Element::GetElement(element),
                                  natoms);
          } else {
            aMaterial->AddElement(G4Element::GetElement(element),
                                  fatoms);
          }
        }
      }
      tokenPtr1 = strtok(NULL, "-");
      //      s.empty();
      //element.erase();
      //
    }

    delete[] sname;
    Material.push_back(aMaterial);
    G4cout << " Material:" << name << " with formula: " << formula << " added! "
           << G4endl;
    G4cout << "     Nb of Material = " << Material.size() << G4endl;
    G4cout << "     Nb of Isotope =  " << Isotope.size() << G4endl;
    G4cout << "     Nb of Element =  " << Element.size() << G4endl;
  }
  ////////////////////////////////////////////////////////////////////////////////
  //
  void KM3Material::DeleteMaterial(G4int j) {
    size_t i(j - 1);
    if (i > Material.size()) {
      G4cerr << "DeleteMaterial : Invalid material index " << j << "."
             << G4endl;
      G4cerr <<"--> Command rejected." <<G4endl;
    } else {
      G4cerr <<"It seems there is no mechanism in G4 for deleting a material yet!"
             <<G4endl;
      G4cerr <<"--> Command rejected." <<G4endl;
    }
  }
  ////////////////////////////////////////////////////////////////////////////////
  //
  void KM3Material::DeleteMaterial(G4String) {
    G4cerr
      << "It seems there is no mechanism in G4 for deleting a material yet!"
      << G4endl;
    G4cerr <<"--> Command rejected." <<G4endl;
  }
  ////////////////////////////////////////////////////////////////////////////////
  //
  G4int KM3Material::GetMaterialIndex(G4String name) {
    size_t i;
    for (i = 0; i < Material.size(); i++) {
      if (Material[i]->GetName() == name)
        break;
    }
    G4int k = G4int(i);
    if (i == Material.size())
      k = -1;
    return k;
  }
  ////////////////////////////////////////////////////////////////////////////////
  //
  void KM3Material::ListMaterial() {
    G4cout << " There are" << std::setw(3) << Material.size()
           << " materials defined." << G4endl;
    for (size_t i = 0; i< Material.size(); i++)
      G4cout <<"     Material Index " <<std::setw(3) <<i+1 <<" "
             <<std::setw(14) <<Material[i]->GetName()
             <<"  density: " <<std::setw(6) <<std::setprecision(3)
             <<G4BestUnit(Material[i]->GetDensity(),"Volumic Mass") <<G4endl;
    G4cout <<G4endl;

  }
  ////////////////////////////////////////////////////////////////////////////////

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void KM3Material::DefineMaterials() {

    //Nothing to do, use the KM3Mat configuration files please!
  }
}
