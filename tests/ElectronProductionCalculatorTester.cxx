#include <TH1F.h>

#include <KM3PMTOpticalModelTester.hh>

#include <G4ThreeVector.hh>
#include <G4DynamicParticle.hh>
#include <G4SauterGavrilaAngularDistribution.hh>
#include <G4OpticalPhoton.hh>

#include <iostream>

using namespace std;
using namespace km3net;

TH1F* h=new TH1F("h", "h", 90,0,90*deg);

int main  (int argc, char** argv)
{
  string sOuputTree;
  float x=0,z=-1;
  int level = 1;
  for (int argIt = 1; argIt < argc; argIt++)
    {
      string sarg = argv[argIt];
      if (sOuputTree == "")
        sOuputTree = sarg;

      else if (sarg == "-angle")
        {
          sarg = argv[++argIt];
          float angle = stof (sarg)*deg;
          x= sin (angle);
          z= -cos (angle);
        }
      else if (sarg == "-el")
        {
          sarg = argv[++argIt];
          level=stoi(sarg);
        }
    }

  if (!sOuputTree.size ())
    {
      cerr << "needs a file.root as argument for the output" << endl;
      return 1;
    }

  KM3PMTOpticalModelTester theModelTester;
  cout << "Create model tester" <<endl;
  theModelTester.LoadMaterial ("photocathode10inches", "PMTGlass", "Vacuum");

  //VLA
  G4ThreeVector pos(0,0,-1);
  G4ThreeVector pol(0,0,-1);

  theModelTester.TestCalculateNormalCoefficients (pos, 440*nm);
  theModelTester.SetElectronProductionLevel (level);
  if (level < 2)
    theModelTester.SetElectronAbsorptionLength (22*nm);


  for (int angle=0;angle<90;angle++)
    {
      x= sin (angle*deg);
      z= -cos (angle*deg);
      G4ThreeVector dir(0,x,z);
      //G4ThreeVector pos(0,0,-1);
      G4ThreeVector norm(0,0,1);
      G4ThreeVector insidedir = -norm;


      float escaped=0;
      float nb_elect=50000;
      for (int i=0; i<nb_elect;i++)
        {
          theModelTester.TestCalculateCoefficients (pos, dir, dir, norm, 440*nm);
          if (theModelTester.TestElectronProduction (dir, pol, norm, insidedir))
            escaped++;
        }

      cout << angle+1 << " " <<escaped/nb_elect << endl;
      h->SetBinContent (angle+1,escaped/nb_elect);
    }

  theModelTester.SaveTreesAndHistograms ("/scratch/hugon/out.root");
  h->SaveAs ("/scratch/hugon/bah.root");

  return 0;

  /*/this part was to check the angular emission of the electron. Left for history ;p
  //G4DynamicParticle DynPart(G4OpticalPhoton::Definition (), dir, 2*eV);
  G4SauterGavrilaAngularDistribution SGAnglular;

  for (int i=0; i<10000; i++)
  {
  G4ThreeVector PEemit = SGAnglular.SampleDirection (&DynPart);
  cout << PEemit.x () << " " << PEemit.y () << " " << PEemit.z () << " " << endl;
  h->Fill (acos (PEemit.z ()));
  }



  return 0;*/
}
