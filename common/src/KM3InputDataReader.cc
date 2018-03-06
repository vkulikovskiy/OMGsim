// This file is part of the GenericLAND software library.
// $Id: KM3InputDataReader.cc,v 1.1.1.1 2004/12/21 22:29:48 KM3sim Exp $
//
// KM3InputDataReader.cc
// v.0 by Glenn Horton-Smith, Feb 12, 1999

#include "KM3InputDataReader.hh"
#include "G4MaterialPropertiesTable.hh"
#include <KM3Material.hh>
#include <G4NistManager.hh>
#include <ctype.h>

#include <map>
#include <string>
#include <iostream>
#include <limits>
using namespace std;

using namespace CLHEP;

namespace km3net
{

  void ToLower (string& str)
  {
    for (auto& c : str)
      {
        c = tolower (c);
      }
  }

  G4Material* KM3InputDataReader::CreateMaterial (string materialName, const map <string, double>& DensityPerComponents, double densityGperCM3, G4State state)
  {
    //for the case of unique material with default density
    G4NistManager* man = G4NistManager::Instance();
    if (densityGperCM3 <= 0)
      {
        if (DensityPerComponents.size () > 1)
          {
            clog << "KM3InputDataReader::CreateMaterial warning, more than 1 component and non default density!"
                 << materialName << endl;
          }
        G4Material* to_return =man->FindOrBuildMaterial(DensityPerComponents.begin ()->first);
        if (to_return == 0)
          {
            cerr << "error! " << DensityPerComponents.begin ()->first << " does not exist in the database!" << endl;
            exit (1);
          }
        return to_return;
      }

    //first the normalization
    double normfactor=0;
    for (const auto& namedensPair : DensityPerComponents)
      {
        normfactor += namedensPair.second;
      }
    //then the material
    G4Material* CreatedMaterial = new G4Material(materialName,
                                                 densityGperCM3 * g / cm3,
                                                 DensityPerComponents.size (),
                                                 state);
    //finishing with the composing
    for (auto& namedensPair : DensityPerComponents)
      {
        CreatedMaterial->AddMaterial (man->FindOrBuildMaterial(namedensPair.first),
                                      namedensPair.second/normfactor);
      }

    return CreatedMaterial;
  }

  int KM3InputDataReader::
  ReadMaterials(std::istream &is)
  {
    static const char funcname[]= "KM3InputDataReader::ReadMaterials";

    G4Material* currentMaterial = NULL;
    G4MaterialPropertiesTable * currentMPT = NULL;
    G4MaterialPropertyVector *currentPV = NULL;

    MyTokenizer t(is);
    int   wavelength_opt= 0;
    int   errorCount=0;
    float shifter_density = 0;
    bool  constproperty=false;

    string propertyname;

    while (t.nextToken() != MyTokenizer::TT_EOF ) {

      // expect either a pair of numbers or a keyword
      if (t.ttype == MyTokenizer::TT_STRING) {
        if (t.sval == "CREATE")
          {
            string materialName;
            if (t.nextToken () == MyTokenizer::TT_STRING)
              {
                materialName=t.sval.data ();
              }
            else
              {
                G4cerr << funcname << "error interpreting a material file. An name should be given to the created material."<< endl;
                exit (1);
              }
            map <string, double> DensityPerComponents;
            double densityGperCM3 = -1;
            G4State state=kStateSolid;
            t.nextToken ();
            while (t.ttype != MyTokenizer::TT_STRING || t.sval != "CREATE")
              {
                if (t.ttype == MyTokenizer::TT_STRING && t.sval == "COMPONENT")
                  {
                    if (t.nextToken () == MyTokenizer::TT_STRING)
                      DensityPerComponents [t.sval] = 1;
                  }
                else if (t.ttype == MyTokenizer::TT_STRING && t.sval == "COMPONENTS")
                  {
                    while (t.nextToken () != MyTokenizer::TT_STRING || t.sval != "COMPONENTS")
                      {
                        if (t.ttype == MyTokenizer::TT_STRING)
                          {
                            string compName=t.sval;
                            if (t.nextToken () == MyTokenizer::TT_NUMBER)
                              {
                                DensityPerComponents [compName]= t.nval;
                              }
                          }
                      }
                  }
                else if (t.ttype == MyTokenizer::TT_STRING && t.sval == "DENSITY")
                  {
                    if (t.nextToken () == MyTokenizer::TT_NUMBER)
                      {
                        densityGperCM3 = t.nval;
                      }
                  }
                else if (t.ttype == MyTokenizer::TT_STRING && t.sval == "STATE")
                  {
                    if (t.nextToken () == MyTokenizer::TT_STRING && t.sval == "gas")
                      {
                        state=kStateGas;
                      }
                    else if (t.ttype == MyTokenizer::TT_STRING && t.sval == "liquid")
                      {
                        state=kStateLiquid;
                      }
                  }
                t.nextToken ();
              }
            currentMaterial = CreateMaterial (materialName, DensityPerComponents, densityGperCM3, state);
            currentMPT= new G4MaterialPropertiesTable();
            currentMaterial->SetMaterialPropertiesTable(currentMPT);
          }
        else if (t.sval == "MATERIAL") {
          if (t.nextToken() == MyTokenizer::TT_STRING) {
            currentMaterial= KM3Material::GetMaterial(t.sval);

            currentPV= NULL;
            wavelength_opt= 0;
            if (currentMaterial == NULL) {
              currentMPT= NULL;
              errorCount++;    // error message issued in GetMaterial
              G4cerr << " No material " << t.sval << endl;
              return errorCount;
            } else {
              currentMPT= currentMaterial->GetMaterialPropertiesTable();
              if (currentMPT == NULL) {
                currentMPT= new G4MaterialPropertiesTable();
                currentMaterial->SetMaterialPropertiesTable(currentMPT);
              }
            }
          }
          else {
            G4cerr << funcname << " expected string after MATERIAL\n";
            errorCount++;
          }
        }
        else if (t.sval == "PROPERTY") {
          wavelength_opt= 0;
          errorCount=0;
          shifter_density = 0;
          constproperty=false;

          if (t.nextToken() == MyTokenizer::TT_STRING) {
            currentPV= NULL;
            propertyname=t.sval;
            if (currentMPT != NULL && !constproperty) {
              currentPV =
                currentMPT->GetProperty((char *)(const char *)(t.sval));
              if (currentPV == NULL) {
                currentPV= new G4MaterialPropertyVector();
                currentMPT->AddProperty((char *)(const char *)(t.sval),
                                          currentPV);
              }
            }
          }
          else {
            G4cerr << funcname << " expected string after PROPERTY\n";
            errorCount++;
          }
        }
        else if (t.sval == "OPTION") {
          if (t.nextToken() == MyTokenizer::TT_STRING) {
            if (t.sval == "wavelength")
              {
                wavelength_opt= 1;
              }
            else if (t.sval == "dy_dwavelength")
              wavelength_opt= 2;
            else if (t.sval == "energy")
              wavelength_opt= 0;
            else if (t.sval == "eV")
              wavelength_opt= 3;
            else if (t.sval == "shifter_density")
              {
                t.nextToken ();
                if (t.ttype != MyTokenizer::TT_NUMBER)
                  {
                    G4cerr << funcname << " expected a number in mol/m3 after shifter_density " << G4endl;
                    errorCount++;
                  }
                else
                  shifter_density = t.nval;
              }
            else if (t.sval == "constant")
              constproperty=true;
            else {
              G4cerr << funcname << " unknown option " << t.sval << G4endl;
              errorCount++;
            }
          }
          else {
            G4cerr << funcname << " expected string after OPTION\n";
            errorCount++;
          }
        }
        else {
          G4cerr << funcname << " unknown keyword " << t.sval << G4endl;
          errorCount++;
        }
      }
      else if (t.ttype == MyTokenizer::TT_NUMBER) {
        double E_value= t.nval;
        if (constproperty)
          {
            currentMPT->RemoveProperty (propertyname.data());
            currentMPT->AddConstProperty(propertyname.data(), E_value);
          }
        else if (t.nextToken() == MyTokenizer::TT_NUMBER) {
          double p_value= t.nval;
          if (currentMPT != NULL && currentPV != NULL) {
            if (wavelength_opt && wavelength_opt < 3) {
              if (E_value != 0.0) {
                double lam= E_value;
                E_value= 2*pi*hbarc/(lam*nanometer);
                if (wavelength_opt == 2)
                  p_value *= lam / E_value;
              }
              else {
                G4cerr << funcname << " zero wavelength!\n";
                errorCount++;
              }
            } else if (wavelength_opt == 3){
              E_value*=eV;
            }
            if (shifter_density != 0)
              {
                if (p_value > 0)
                  {
                    p_value=1./shifter_density/p_value*cm;
                  }
                else
                  p_value=numeric_limits<double>::max ();
              }
            currentPV->InsertValues(E_value, p_value);
          }
          else {
            G4cerr << funcname << currentMaterial<< " got number pair, but have no pointer to ";
            if (currentMPT == NULL) G4cerr << "MaterialPropertyTable ";
            if (currentPV == NULL) G4cerr << "MaterialPropertyVector ";
            G4cerr << G4endl;
            errorCount++;
          }
        }
        else {
          G4cerr << funcname
                 << " expected second number, but tokenizer state is ";
          t.dumpOn(G4cerr);
          G4cerr << G4endl;
          errorCount++;
        }
      }
      else {
        G4cerr << funcname
               << " expected a number or a string, but tokenizer state is ";
        t.dumpOn(G4cerr);
        G4cerr << G4endl;
        errorCount++;
      }
    }

    return errorCount;
  }


  void KM3InputDataReader::MyTokenizer::
  dumpOn(std::ostream &os)
  {
    os << "KM3InputDataReader::MyTokenizer[ttype="
       << ttype
       << ",nval="
       << nval
       << ",sval="
       << sval
       << "] ";
  }


  int
  KM3InputDataReader::MyTokenizer::
  nextToken(void)
  {
    int i=0;
    G4bool negateFlag=false;
    do {
      i= isptr->get();
      if (i=='+')
        i= isptr->get();
      if (i=='-') {
        i= isptr->get();
        negateFlag= !negateFlag;
      }
      if (i=='#') {  // comment to end of line
        do {
          i= isptr->get();
        } while (i!=EOF && i!='\n');
      }
      if (i== EOF)
        return (ttype=TT_EOF);
    } while (isspace(i));

    if (isdigit(i) || i=='.') {
      nval= 0.0;
      isptr->putback(i);
      (*isptr) >> nval;
      if (negateFlag)
        nval= -nval;
      return (ttype=TT_NUMBER);
    }
    else if (negateFlag) {
      isptr->putback(i);
      return (ttype='-');
    }
    else if (isalpha(i) || i=='_') {
      isptr->putback(i);
      (*isptr) >> sval;
      return (ttype=TT_STRING);
    }
    else if (i=='"') {
      sval="";
      while (true) {
        i= isptr->get();
        while (i=='\\') {
          i= isptr->get();
          sval.append((char)i);
          i= isptr->get();
        }
        if (i==EOF || i=='"')
          break;
        sval.append((char)i);
      }
      return (ttype=TT_STRING);
    }
    else {
      return (ttype=i);
    }
  }
}
