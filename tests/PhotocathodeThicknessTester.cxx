#include <TH1I.h>
#include "TGraph.h"
#include "TFile.h"

#include <KM3PMTOpticalModelTester.hh>

#include <TerminalHistogram.hh>
#include <iostream>



//to set the random engine
#include "Randomize.hh"
#include <sys/time.h>
#include <unistd.h>
#include <random>

using namespace std;
using namespace km3net;

KM3PMTOpticalModelTester theModelTester;

bool PercentTo (float percent, TH1* histo, float value, int binmin, int binmax)
{
  /*
    float total = 0;
    for (int bin = binmin; bin<binmax;bin++)
    total+=abs(histo->GetBinContent (bin)-value);

    total/=(binmax-binmin);
    cout << "atotal " << total << endl;
    return (total < percent);
  */


  float worse = 0;
  float diff;
  for (int bin = binmin; bin<binmax;bin++) {
    diff = abs(histo->GetBinContent (bin)-value);
    if (worse>diff) worse= diff;
  }

  cout << "worse " << worse << endl;
  return (worse < percent);
}

char CompareQuantumEfficiency (float wavelength)
{
  G4ThreeVector pos(0,0,1);
  G4ThreeVector norm(0,0,1);
  G4ThreeVector insideDirection (0,0,-1);
  G4ThreeVector dir(0,0,-1);
  G4ThreeVector pol(0,0,1);

  float compareit (theModelTester.GetEfficiency(wavelength) -
                   theModelTester.GetQuantumEfficiency (pos,dir,pol, norm, insideDirection, wavelength));

  cout << "compareit " << compareit << endl;

  if (compareit + 0.01 < 0)
    return -1;
  if (compareit - 0.01 > 0)
    return 1;
  return 0;

}

int main (int argc, char** argv)
{

	float minwavelength = 400;
	float maxwavelength = 800;
	float mincathodethick = 1;
  float maxcathodethick = 1000;
	float theta = 0;

  int level = 5;
  float precision = 0.005;
  for (int argIt = 1; argIt < argc; argIt++)
    {
      string sarg = argv[argIt];

      if (sarg == "-w")
        {
          sarg = argv[++argIt];
          minwavelength  = stof (sarg);
					maxwavelength  = minwavelength;
        }
			if (sarg == "-theta")
			{
				sarg = argv[++argIt];
				theta  = stof (sarg);
			}
      else if (sarg == "-level")
        {
          sarg = argv[++argIt];
          level=stoi(sarg);
        }
      else if (sarg == "-cathode")
        {
          sarg = argv[++argIt];
          mincathodethick=stoi(sarg);
					maxcathodethick = mincathodethick;
        }
      else if (sarg == "-precision")
        {
          sarg = argv[++argIt];
          precision = stof (sarg);
        }
    }


  // automatic (time-based) random seeds and filenames for each run
  CLHEP::HepRandom::setTheEngine(new CLHEP::DualRand);
  timeval tim;
  gettimeofday(&tim, NULL);
  long ttime = tim.tv_sec xor tim.tv_usec;
  G4int pid = getpid();
  long Seed = ttime xor (pid << 8);
  long seeds[2];
  seeds[0] =  Seed;
  seeds[1] =  (long)(Seed*G4UniformRand());
  CLHEP::HepRandom::getTheEngine()->setSeeds(seeds,-1);
  CLHEP::HepRandom::showEngineStatus();

  bool calibrateQE = false;
  float glassThickness = 1*mm;

  theModelTester.LoadMaterial ("photocathode3inches", "PMTGlass_R12199_02", "Vacuum");
  theModelTester.fGlassThickness =glassThickness;
  theModelTester.SetElectronProductionLevel(level);
  theModelTester.SetElectronAbsorptionLength (1000*nm);

  if (precision != 0) theModelTester.fElPrecision = 25./(precision*precision);  //this is the number of simulated e-, precision is aprrox 5~sigma and sigma is sqrt(N_e)
  else theModelTester.fElPrecision = 0;

  TerminalHistogram *measuredQE =0;TerminalHistogram *simulatedQE =0; TerminalHistogram *ratioQE = 0;


  double optcathodethick = 0;
  double optqe = 0;
  double qe;

  for (float wavelength = minwavelength; wavelength <= maxwavelength; wavelength+=50) {
    optcathodethick = 0;
    optqe = 0;
    for (double cathodethick = mincathodethick; cathodethick <= maxcathodethick; cathodethick+=1) {
      G4MaterialPropertyVector newThicknessVector;
      newThicknessVector.InsertValues (0, cathodethick*1e-6);
      newThicknessVector.InsertValues (188, cathodethick*1e-6);
      theModelTester.SetCathodeThickness(newThicknessVector);
      qe = theModelTester.QECheck(wavelength,theta);
      if (qe < 0) {
        //cout << "WARNING: " << wavelength << " " << cathodethick << " " << qe << endl;
      }
      if (qe > optqe) {
        optqe = qe;
        optcathodethick = cathodethick;
      }
      //cout << wavelength << " " << cathodethick << " " << qe << endl;
    }
    cout << wavelength << " " << optcathodethick << " " << optqe << endl;
  }
	return 0;

}
