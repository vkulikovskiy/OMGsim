#include "TGraph.h"
#include "TFile.h"
#include <TH1I.h>

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
  G4ThreeVector insideDirection(0,0,-1);
  G4ThreeVector dir(0,0,-1);
  G4ThreeVector pol(0,0,1);
	MixedPolarization(dir,pol);

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

	CLHEP::HepRandom::setTheEngine(new CLHEP::DualRand);
	//CLHEP::HepRandom::restoreEngineStatus();
	//G4int ttime = time(0);
	//G4cout << "ttime=" << ttime << G4endl;
	timeval tim;
	gettimeofday(&tim, NULL);
	long ttime = tim.tv_sec xor tim.tv_usec;
	G4cout << "ttime=" << ttime << G4endl;

	G4int pid = getpid();
	G4cout << "pid=" << pid << G4endl;

	long Seed = ttime xor (pid << 8);
	// automatic (time-based) random seeds and filenames for each run
	G4cout << "******************" << G4endl;
	G4cout << "*** AUTOSEED ON ***" << G4endl;
	G4cout << "*******************" << G4endl;
	long seeds[2];
	seeds[0] =  Seed;
	seeds[1] =  (long)(Seed*G4UniformRand());
	G4cout << "seed1: " << seeds[0] << "; seed2: " << seeds[1] << G4endl;

	CLHEP::HepRandom::getTheEngine()->setSeeds(seeds,-1);
	CLHEP::HepRandom::showEngineStatus();


  bool calibrateQE = false;
  //float glassThickness = 1*mm;  //for big 
	float glassThickness = 2*mm;
	float precision = 0.05;

  for (int argIt = 1; argIt < argc; argIt++)
    {
      string sarg = argv[argIt];
      if (sarg == "-cal")
        calibrateQE = true;
      if (sarg == "-PhThick")
        {
          sarg = argv[++argIt];
          glassThickness = stof (sarg)*mm;
        }
      if (sarg == "-precision")
        {
          sarg = argv[++argIt];
          precision = stof (sarg);
        }
    }

  cout << "Create model tester" <<endl;
  //theModelTester.LoadMaterial ("photocathode3inches", "PMTGlass", "Vacuum");
	theModelTester.LoadMaterial ("photocathode3inches", "PMTGlass_R12199_02", "Vacuum");
  theModelTester.fGlassThickness =glassThickness;
	theModelTester.SetElectronProductionLevel(3);
	theModelTester.SetElectronAbsorptionLength (1000*nm);

  if (precision != 0) theModelTester.fElPrecision = 25./(precision*precision);  //this is the number of simulated e-, precision is approx 5~sigma and sigma is sqrt(N_e)
  else {
		cout << "setting theModelTester.fElPrecision to 0" << endl;
		theModelTester.fElPrecision = 0;
		precision = 0.01;   //calibration precision arbitrary set to 1% 
	}


  TerminalHistogram *measuredQE =0;TerminalHistogram *simulatedQE =0; TerminalHistogram *ratioQE = 0;

	TGraph *gratioQE;


/*
//BEGIN VLA
  double _wavelength = 300;
	TGraph *g_QE = new TGraph();
	ostringstream oss;
	oss << "QE" << _wavelength;
	g_QE->SetName(oss.str().c_str());

	G4MaterialPropertyVector oldKindexVector = theModelTester.GetKIndex ();
	G4ThreeVector pos(0,0,1);
  G4ThreeVector norm(0,0,-1);
  G4ThreeVector dir(0,0,-1);
  G4ThreeVector pol(0,0,1);
	for (double kindex = 0.1; kindex < 10; kindex+=0.01) {
		G4MaterialPropertyVector newKindexVector;
		for (double wavelength = 250; wavelength <= 700; wavelength+=10) {
			float energy = twopi*hbarc/wavelength/nm;
			float newKindex = 0;
			if (wavelength != _wavelength) newKindex = oldKindexVector.Value (energy);
			else newKindex =  kindex;
			//cout << wavelength << " " << newKindex << endl;
			newKindexVector.InsertValues (energy, newKindex);
		}
		theModelTester.SetKIndex (newKindexVector);
		g_QE->SetPoint(g_QE->GetN(),kindex, theModelTester.GetQuantumEfficiency (pos,dir,pol, norm, _wavelength));
	}
	TFile *fout = new TFile("debug.root","RECREATE");
	g_QE->Write();
	fout->Close();
//END VLA
*/


	TerminalHistogram* simulatedQE_prev = NULL;

	cout << "starting calc " << endl;

	bool continueopt;

  do {
		cout << "DEBUG: inside the loop. " << endl;
		continueopt = false;
		/*
    if (measuredQE != 0)
      {
        measuredQE->Reset ();
        simulatedQE->Reset ();
      }
			*/
    theModelTester.QuantumEfficiencyMeasurement ();
    theModelTester.PhotoCathodeEfficiencyReader ();

    measuredQE  = (TerminalHistogram*)(theModelTester.GetTObject ("MaterialEfficiencyHisto"));
    simulatedQE = (TerminalHistogram*)(theModelTester.GetTObject ("QuantumEfficiencyHisto"));

		cout << "DEBUG: Materials are read." << endl;


		if (!theModelTester.InTObject("RatioQE")) {
			ratioQE = new TerminalHistogram ("RatioQE","QE sim / QE meas;nm",measuredQE->GetNbinsX(),measuredQE->GetXaxis()->GetXmin(),measuredQE->GetXaxis()->GetXmax());
			//ratioQE = (TerminalHistogram*)measuredQE->Clone("RatioQE");
			gratioQE = new TGraph();
			gratioQE->SetName("gRatioQE");
			theModelTester.PutTObject("gRatioQE",gratioQE);
			theModelTester.PutTObject("RatioQE",ratioQE);
    }
    else
      ratioQE = (TerminalHistogram *)(theModelTester.GetTObject("RatioQE"));
		ratioQE->Reset();

    G4MaterialPropertyVector oldKindexVector = theModelTester.GetKIndex ();
    G4MaterialPropertyVector newKindexVector;
		double step = 1.1;
    cout << "wavelength kindex " << endl;
    for (int bin = 1; bin <= measuredQE->GetNbinsX (); bin++)
      {
        float newbin=0;
        if (measuredQE->GetBinContent (bin) != 0)
          newbin=simulatedQE->GetBinContent (bin)/measuredQE->GetBinContent (bin);
        //cout << newbin << endl;
        ratioQE->SetBinContent (bin, newbin);
        float wavelength = ratioQE->GetBinLowEdge (bin);
        float energy = twopi*hbarc/wavelength/nm;
        float newKindex = 0;

				if (newbin != 0) gratioQE->SetPoint(gratioQE->GetN(),energy/eV,1./newbin);

        newKindex = oldKindexVector.Value (energy);

        if (wavelength > 250 && newbin != 0)
          {
            if (newbin < 1-precision){
							if (simulatedQE_prev!=0) {
								if (simulatedQE->GetBinContent (bin) / simulatedQE_prev->GetBinContent (bin) <= 1. + precision/5.) {    //stop to increase (stays stable at the level of 1 sigma fluctuation)
									cout << "*";  //means that optimisation stopped since QE started to decrease
								} else {
									newKindex*=step;
									cout << "!"; //means that optimisation is going on
									continueopt = true;
								}
							} else {
								newKindex*=step;
								cout << "!"; //means that optimisation is going on
								continueopt = true;
							}
						}
						if (newbin > 1+precision) {
							//cout << "check: " << newKindex << " " << newbin << " " << newKindex*(1/newbin) << endl;
							cout << "!"; //means that optimisation is going on
              newKindex = newKindex*(1/newbin);
							continueopt = true;
						}
          }


        newKindexVector.InsertValues (energy, newKindex);
				cout << wavelength << " " << newKindex << endl;

        //cout <<wavelength << " old " << oldKindexVector.Value (energy) << " " <<  newKindex << endl;
      }
    theModelTester.SetKIndex (newKindexVector);
		cout << "Ratio of calculated QE to QE in the data" << endl;
    ratioQE->Draw ();

		if (simulatedQE_prev != NULL) delete simulatedQE_prev;
		simulatedQE_prev=(TerminalHistogram*)simulatedQE->Clone();

	}while (continueopt && calibrateQE);
	if (calibrateQE) cout << "precision " << precision  << " reached" << endl;
	theModelTester.SaveTreesAndHistograms ("testerout.root");
  // ( !PercentTo (0.02,measuredQE, 1, 11, 42) && calibrateQE);
  //}while ( !PercentTo (0.02,measuredQE, 1.1, 11, 42) && calibrateQE);



  return 0;
}
