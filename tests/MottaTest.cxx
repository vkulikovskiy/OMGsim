#include <TH1F.h>
#include <TFile.h>
#include <TRandom.h>


#include <KM3PMTOpticalModelTester.hh>

#include <G4ThreeVector.hh>
#include <G4DynamicParticle.hh>
#include <G4SauterGavrilaAngularDistribution.hh>
#include <G4OpticalPhoton.hh>

#include <iostream>

using namespace std;
using namespace km3net;

int nbins=900;
TH1F* hElectEscape=new TH1F("hElectEscape", "Electron escaping probability;deg", nbins,0,90);
TH1F* hPhotonAbs=new TH1F("hPhotonAbs", "Electron absorption probability;deg", nbins,0,90);
TH1F* hPhotonRefl=new TH1F("hPhotonRefl", "Electron absorption probability;deg", nbins,0,90);
TH1F* hPhotonTrans=new TH1F("hPhotonTrans", "Electron absorption probability;deg", nbins,0,90);

TH1F* hPhotonAbsP=new TH1F("hPhotonAbsP", "Electron absorption probability;deg", nbins,0,90);
TH1F* hPhotonReflP=new TH1F("hPhotonReflP", "Electron absorption probability;deg", nbins,0,90);
TH1F* hPhotonTransP=new TH1F("hPhotonTransP", "Electron absorption probability;deg", nbins,0,90);

TH1F* hPhotonAbsS=new TH1F("hPhotonAbsS", "Electron absorption probability;deg", nbins,0,90);
TH1F* hPhotonReflS=new TH1F("hPhotonReflS", "Electron absorption probability;deg", nbins,0,90);
TH1F* hPhotonTransS=new TH1F("hPhotonTransS", "Electron absorption probability;deg", nbins,0,90);

TH1F* hPhotonAbsM=new TH1F("hPhotonAbsM", "Electron absorption probability;deg", nbins,0,90);
TH1F* hPhotonReflM=new TH1F("hPhotonReflM", "Electron absorption probability;deg", nbins,0,90);
TH1F* hPhotonTransM=new TH1F("hPhotonTransM", "Electron absorption probability;deg", nbins,0,90);

TH1F* hPhotonAbsR=new TH1F("hPhotonAbsR", "Electron absorption probability;deg", nbins,0,90);
TH1F* hPhotonReflR=new TH1F("hPhotonReflR", "Electron absorption probability;deg", nbins,0,90);
TH1F* hPhotonTransR=new TH1F("hPhotonTransR", "Electron absorption probability;deg", nbins,0,90);




TRandom TR;

void RandomPolarization (const G4ThreeVector &dir, G4ThreeVector &pol)
{
  double theta = dir.getTheta();
  double phi = dir.getPhi();
  //pol.set(0,1,0);
  //pol.set(1,0,0);
  double poltheta = G4UniformRand()*twopi;
  pol.set(cos(poltheta),sin(poltheta),0);
  pol.rotateY(theta);
  pol.rotateZ(phi);
}

void PPolarization (const G4ThreeVector &dir, G4ThreeVector &pol)
{
  pol.set (0, 1, 0);
  pol=dir.cross (pol);
  pol = pol.unit();
}
void SPolarization (const G4ThreeVector &dir, G4ThreeVector &pol)
{
  pol.set (1, 0, 0);
  pol=dir.cross (pol);
  pol = pol.unit();
}

void MixedPolarization (const G4ThreeVector &dir, G4ThreeVector &pol)
{
  G4ThreeVector pols, polp;
  SPolarization(dir,pols);
  PPolarization(dir,polp);
  pol = (pols+polp);
  pol = pol.unit();

}

int main  (int argc, char** argv)
{
  string sOutputTree="MottaTest.root";
  float x=0,z=-1;
  int level = 1;
  bool MottaCalc = true;
  for (int argIt = 1; argIt < argc; argIt++)
    {
      string sarg = argv[argIt];
      if (sarg == "-f") {
        sarg = argv[++argIt];
        sOutputTree = sarg;
      }
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
      else if (sarg == "-Motta")
        {
          sarg = argv[++argIt];
          cout << "sarg " << sarg << endl;
          istringstream(sarg) >> std::boolalpha >> MottaCalc;
        }

    }

  if (!sOutputTree.size ())
    {
      cerr << "needs a file.root as argument for the output" << endl;
      return 1;
    }

  if (MottaCalc) cout << "Motta calculation with the scintillator (for reflection no first reflection on scintillator)." << endl;
  else cout << "Calculations for A,R,T for the glass+photocathode+vacuum layer only." << endl;

  KM3PMTOpticalModelTester theModelTester;
  cout << "Create model tester" <<endl;
  theModelTester.LoadMaterial ("PhotoMotta", "PMTGlassMotta", "Vacuum");

  //VLA
  G4ThreeVector pos(0,0,-1);
  G4ThreeVector pol(0,0,1);
  theModelTester.TestCalculateNormalCoefficients (pos, 442*nm);
  theModelTester.SetElectronProductionLevel (level);

  theModelTester.SetElectronAbsorptionLength (220000*nm);


  double angle;
  bool nofirstref = true;
  double F;
  for (double angle1=hPhotonAbsS->GetBinCenter(0);angle1<hPhotonAbsS->GetXaxis()->GetXmax();angle1+=hPhotonAbsS->GetBinWidth(0))
    {
      angle = sin(angle1*deg)*1.50/1.51;
      if (angle < 1) {
        angle = asin(angle)/deg;
        x= sin (angle*deg);
        z= -cos (angle*deg);
        G4ThreeVector dir(x,0,z);
        SPolarization (dir, pol);
        G4ThreeVector norm(0,0,1);
        G4ThreeVector insideDirection (0,0,-1);

        float escaped=0;
        float nb_elect=50000;
        theModelTester.TestCalculateCoefficients (pos, dir, pol, norm, 442*nm);
        double Abs,Trans,Refl;
        double AbsT,TransT,ReflT;
        theModelTester.GetTotalProbabilities (Abs, Trans, Refl);


        if (MottaCalc) {
          F = pow(sin(angle1*deg-angle*deg)/sin(angle1*deg+angle*deg),2);
          Trans= Trans*(1-F)/(1-F*Refl);
          Refl = F+Refl*(1.-F)*(1.-F)/(1-F*Refl);
          Abs = 1-Trans-Refl;
          Refl -= F;
        }


        hPhotonAbsS->Fill (angle1,Abs);
        hPhotonTransS->Fill (angle1,Trans);
        hPhotonReflS->Fill (angle1,Refl);


        AbsT=Abs; TransT=Trans; ReflT=Refl;
        PPolarization (dir, pol);
        theModelTester.TestCalculateCoefficients (pos, dir, pol, norm, 442*nm);
        theModelTester.GetTotalProbabilities (Abs, Trans, Refl);

        if (MottaCalc) {
          F = pow(tan(angle1*deg-angle*deg)/tan(angle1*deg+angle*deg),2);
          Trans= Trans*(1-F)/(1-F*Refl);
          Refl = F+Refl*(1.-F)*(1.-F)/(1-F*Refl);
          Abs = 1-Trans-Refl;
          Refl -= F;
        }

        hPhotonAbsP->Fill (angle1,Abs);
        hPhotonTransP->Fill (angle1,Trans);
        hPhotonReflP->Fill (angle1,Refl);

        AbsT+=Abs; TransT+=Trans; ReflT+=Refl;
        hPhotonAbs->Fill (angle1,AbsT/2);
        hPhotonTrans->Fill (angle1,TransT/2);
        hPhotonRefl->Fill (angle1,ReflT/2);


        MixedPolarization (dir, pol);
        theModelTester.TestCalculateCoefficients (pos, dir, pol, norm, 442*nm);
        theModelTester.GetTotalProbabilities (Abs, Trans, Refl);
        hPhotonAbsM->Fill (angle1,Abs);
        hPhotonTransM->Fill (angle1,Trans);
        hPhotonReflM->Fill (angle1,Refl);

        double steps = 100;
        for (int ii = 0; ii < steps; ii++) {
          RandomPolarization(dir,pol);
          theModelTester.TestCalculateCoefficients (pos, dir, pol, norm, 442*nm);
          theModelTester.GetTotalProbabilities (Abs, Trans, Refl);
          hPhotonAbsR->Fill(angle1,Abs/steps);
          hPhotonTransR->Fill(angle1,Trans/steps);
          hPhotonReflR->Fill(angle1,Refl/steps);
        }

        MixedPolarization (dir, pol);
        float Pescaped = theModelTester.TestElectronEscapingProbability (dir, pol,norm);
        //cout << angle+1 << " " << Pescaped << endl;
        hElectEscape->Fill (angle1,Pescaped);
      }
    }

  //theModelTester.SaveTreesAndHistograms ("/scratch/hugon/out.root");
  TFile tf(sOutputTree.data (),"recreate");
  tf.Add (hElectEscape);
  tf.Add (hPhotonAbs);  tf.Add (hPhotonTrans);  tf.Add (hPhotonRefl);
  tf.Add (hPhotonAbsS);  tf.Add (hPhotonTransS);  tf.Add (hPhotonReflS);
  tf.Add (hPhotonAbsP);  tf.Add (hPhotonTransP);  tf.Add (hPhotonReflP);
  tf.Add (hPhotonAbsM);  tf.Add (hPhotonTransM);  tf.Add (hPhotonReflM);
  tf.Add (hPhotonAbsR);  tf.Add (hPhotonTransR);  tf.Add (hPhotonReflR);
  tf.Write ();
  tf.Close ();

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
