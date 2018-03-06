#define NoiseCoincidences_cxx
#define RunInfo_cxx

#include "NoiseCoincidences.h"
#include "RunInfo.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <bitset>
#include <iostream>
#include <vector>
#include <TVector3.h>

#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TString.h>

using namespace std;

TH2F* hpairlist;
TH1F *hcoincidences,*hcoincidencesMin,*hcoincidencesMax;
TH1F* hPhotonsInPM;
TGraph* gAngularCoincidences=0;

//new version
float PMTPositions[31][2] = {
  {0.,0.}
  ,{0.56,0.524}, {0.56,1.571},{0.56,2.62},{0.56,3.665},{0.56,4.712},{0.56,5.76}
  ,{0.98,0.},{0.98,1.047},{0.98,2.094},{0.98,3.142}, {0.98,4.1889},{0.98,5.236}
  ,{1.27,0.524}, {1.27,1.571},{1.27,2.62},{1.27,3.665},{1.27,4.712},{1.27,5.76}
  ,{1.872,5.236},{1.872,0.},{1.872,1.047},{1.872,2.094},{1.872,3.142},{1.872,4.1889}
  ,{2.162,5.76},{2.162,0.524},{2.162,1.571},{2.162,2.62},{2.162,3.665},{2.162,4.712}
};
                            //*/
  /*  //old version
  {{0.,0.}
   ,{0.56,0.}, {0.56,1.047},{0.56,2.094},{0.56,3.142},{0.56,4.1889},{0.56,5.236}
   ,{0.98,0.524},{0.98,1.571},{0.98,2.62},{0.98,3.665}, {0.98,4.712},{0.98,5.76}
   ,{1.27,0.}, {1.27,1.047},{1.27,2.094},{1.27,3.142},{1.27,4.1889},{1.27,5.236}
   ,{1.872,0.524},{1.872,1.571},{1.872,2.62},{1.872,3.665},{1.872,4.712},{1.872,5.76}
   ,{2.162,0.},{2.162,1.047},{2.162,2.094},{2.162,3.142},{2.162,4.1889},{2.162,5.236}
  };*/



unsigned myPow(unsigned x, unsigned p) {
  if (p == 0) return 1;
  if (p == 1) return x;
  return x * myPow(x, p-1);
}

unsigned ConvertPMMask (unsigned PM)
{
  if (PM == 0)
    return 0;
  return myPow (2,PM-1);
}

unsigned GetPairId (unsigned PM1, unsigned PM2)
{
  if (PM1 > PM2)
    swap (PM1,PM2);

  // if (PM2 == 21)
  //   swap (PM1,PM2);


  //       unsigned j=PM2-PM1;

  //       unsigned tmp=0;
  //       for (int it=0;it<j-1;it++)
  //         tmp+=it;

  //       return (j-1)*31-tmp+PM1;


  //return (j-1)*31-tmp+PM1+1;
  unsigned PMR=(31-PM2);
  return 31*PMR-(PMR+1)*PMR/2+(PM2-PM1);

}

#include <random>
float RandomFactor=1;
std::default_random_engine generator;
std::normal_distribution<double> distribution(0,1.8);
std::uniform_real_distribution<double> Udistribution(0.0,1.0);

void FillPairHistogramm (const float* Times)
{
  for (int first=1; first<31;first++)
    for (int second=first+1; second<32;second++)
      {
        if (Times[first] < 0. || Times[second] < 0)
          continue;
        double dTTS = distribution(generator)-distribution(generator);

        hpairlist->Fill(dTTS+Times[first]-Times[second], GetPairId (first, second));
        //hpairlist->Fill(Times[second]-Times[first], GetPairId (second, first));
      }


}

TH1D* hpairlistProj;

void GetAngularCoincidences ()
{
  gAngularCoincidences=new TGraph;
  gAngularCoincidences->SetTitle("Angular coincidences;deg;Hz");
  gAngularCoincidences->SetName("AngularCoincidences");
  int pointcounter=0;

  hpairlistProj=hpairlist->ProjectionY("hpairlistProj");

  for (int PM1=1; PM1<31;PM1++)
    for (int PM2=PM1+1; PM2<32;PM2++)
      {
        TVector3 tvPM1(0,0,1),tvPM2(0,0,1);
        tvPM1.SetMagThetaPhi(1,PMTPositions[PM1-1][0],PMTPositions[PM1-1][1]);
        tvPM2.SetMagThetaPhi(1,PMTPositions[PM2-1][0],PMTPositions[PM2-1][1]);
        float angle=tvPM2.Angle(tvPM1);
        int bin=hpairlistProj->FindBin (GetPairId(PM1,PM2));
        gAngularCoincidences->SetPoint(pointcounter++, angle*180/M_PI, hpairlistProj->GetBinContent (bin));
        //cout << PM1 << ' ' << PM2<< " " << hpairlistProj->GetBinContent(bin)<<endl;
      }
}

void NoiseCoincidences::InitVar ()
{
  ftphotons=0;
  fcoincidences[0] = 0; fcoincidences[1] = 0; fcoincidences[2] = 0;
  fPMMask[0]=0;  fPMMask[1]=0;  fPMMask[2]=0;
  memset (fHitTimePerPM,0xaa,sizeof(fHitTimePerPM)); //put everything to less than 0
}

void NoiseCoincidences::FillHistograms ()
{
  hPhotonsInPM->Fill(ftphotons);
  FillPairHistogramm (fHitTimePerPM);
  hcoincidencesMin->Fill (fcoincidences[0]);
  hcoincidences->Fill (fcoincidences[1]);
  hcoincidencesMax->Fill (fcoincidences[2]);

	//Vla
	//if (fcoincidences[1] > 3) 
	//cout << "coinc: " << fcoincidences[1] << endl;

}

void NoiseCoincidences::AddEvent ()
{
  ftphotons++;
  unsigned CurrentPMMask=ConvertPMMask(Vertex_PMID);
	//Vla - ANTARES
	//unsigned CurrentPMMask=ConvertPMMask(Vertex_OMID);

  
  if (Hit_Limit == 1 || Hit_Limit == 2 || Hit_Limit == 3) {
		if (!(fPMMask[2] & CurrentPMMask)) {
      fcoincidences[2]++;
      fPMMask[2]|=CurrentPMMask;
		}
		if (Hit_Limit == 2 || Hit_Limit == 3) {
			if (!(fPMMask[1] & CurrentPMMask)) {
				 fcoincidences[1]++;
				fPMMask[1]|=CurrentPMMask;
				fHitTimePerPM[Vertex_PMID]=Hit_TotalTime;
			}
			if (Hit_Limit == 3 && !(fPMMask[0] & CurrentPMMask)) {
				fcoincidences[0]++;
				fPMMask[0]|=CurrentPMMask;
			}
		}
	}

	
/*
  if (Hit_Limit == 3)
    {
			if (!(fPMMask[0] & CurrentPMMask)) {
				fcoincidences[0]++;  
				fPMMask[0]|=CurrentPMMask;
			}
      if (!(fPMMask[1] & CurrentPMMask)) {
        fcoincidences[1]++;  
        fPMMask[1]|=CurrentPMMask;
				fHitTimePerPM[Vertex_PMID]=Hit_TotalTime;
      }
      if (!(fPMMask[2] & CurrentPMMask)) {
        fcoincidences[2]++;  
        fPMMask[2]|=CurrentPMMask;
      }
    }else if (Hit_Limit == 2)
    {
      if (!(fPMMask[1] & CurrentPMMask)) {
        fcoincidences[1]++;  
        fPMMask[1]|=CurrentPMMask;
        fHitTimePerPM[Vertex_PMID]=Hit_TotalTime;
      }
      if (!(fPMMask[2] & CurrentPMMask)) {
        fcoincidences[2]++;  
        fPMMask[2]|=CurrentPMMask;
      }
    } else if (Hit_Limit == 1 && !(fPMMask[2] & CurrentPMMask))
    {
      fcoincidences[2]++;
      fPMMask[2]|=CurrentPMMask;
    }
*/

		//cout << "limit, omid " << (int)Hit_Limit << " " << Vertex_OMID << endl;

}
unsigned NoiseCoincidences::Loop()
{
	cout << "RandomFactor " << RandomFactor << endl;
   if (fChain == 0) return 0;

   Long64_t nentries = fChain->GetEntriesFast();

   float previousTime;
   unsigned previousPMID;

   unsigned PMMask=0;
   unsigned CurrentPMMask=0;

   unsigned tphotons=0;
   unsigned HitCounter=0;

   float HitTimePerPM[32];
   memset (HitTimePerPM,0xaa,sizeof(HitTimePerPM)); //put everything to less than 0

   ULong64_t CurrentHitNumber=0xffffffff;
   unsigned short coincidences[3]={1,1,1};

   hcoincidences=new TH1F("hcoincidences","Folds",30,0,30);
	 hcoincidences->Sumw2();
   hcoincidencesMax=new TH1F("hcoincidencesMax","Folds",30,0,30);
   hcoincidencesMin=new TH1F("hcoincidencesMin","Folds",30,0,30);
   hPhotonsInPM=new TH1F("hPhotonsInPM","Folds",2000,0,2000);
   hpairlist=new TH2F("hpairlist","coincidences time per pair id",140,-35,35,500,0,500);

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries
          ;jentry++) {
     if (jentry%10000 == 0)
       {
         cout << "\r" <<(float)jentry/nentries;
         cout.flush ();
       }
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
			//cout << HitNumber << " OM " << Vertex_OMID << " PMT " << Vertex_PMID << " limit " << (int)Hit_Limit  << endl;
      if (Hit_Limit == 4 || (RandomFactor != 1 && Udistribution(generator) > RandomFactor)) continue;



      if (HitNumber != CurrentHitNumber)
        {

          FillHistograms ();
          InitVar ();
          AddEvent ();

          HitCounter++;
          CurrentHitNumber=HitNumber;
        }
      else
        {
          AddEvent ();
        }


      /*      if (HitNumber != CurrentHitNumber)
        {
          hPhotonsInPM->Fill(tphotons);
          FillPairHistogramm (HitTimePerPM);

          tphotons=0;
          HitCounter++;
          //cout << coincidences << endl;
          hcoincidencesMin->Fill (coincidences[0]);
          hcoincidences->Fill (coincidences[1]);
          hcoincidencesMax->Fill (coincidences[2]);

          coincidences[0] = 1; coincidences[1] = 1; coincidences[2] = 1;
          CurrentHitNumber=HitNumber;
          PMMask=ConvertPMMask(Vertex_PMID);

          memset (HitTimePerPM,0xaa,sizeof(HitTimePerPM)); //put everything to less than 0
          if (Hit_Limit > 1)
            HitTimePerPM[Vertex_PMID]=Hit_TotalTime;


          //previousTime=0;
          //previousPMID=Vertex_PMID;
        }
      else
        {
          tphotons++;
          CurrentPMMask=ConvertPMMask(Vertex_PMID);
          //cout << "CurrentPMMask "<< hex << CurrentPMMask << " " <<dec<< Vertex_PMID << endl;
          if (! (PMMask & CurrentPMMask))
            {

              if (Hit_Limit > 1)
                HitTimePerPM[Vertex_PMID]=Hit_TotalTime;
              //previousTime=Hit_TotalTime;
              //previousPMID=Vertex_PMID;

              if (Hit_Limit == 3)
                {
                  coincidences[0]++;
                  coincidences[1]++;
                  coincidences[2]++;
                }
              else if (Hit_Limit == 2)
                {
                  coincidences[1]++;
                  coincidences[2]++;
                }
              else if (Hit_Limit == 1)
                {
                  coincidences[2]++;
                }
              PMMask = CurrentPMMask | PMMask;
            }
        }*/

   }
   return HitCounter;

}

double RunInfo::Loop()
{
  if (fChain == 0) return 0;

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  double NbTriggTotal=0;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    NbTriggTotal+=double(NBTrigg)/TargetVolume; //number of generated events (decays of K40 for example) divided on generated volume (in m^3)
    //NbTriggTotal+=double(NBTrigg);            //number of generated events for radioactivity in the glass study
    //cout <<"norm ici " << NbTriggTotal << endl;
  }
  return NbTriggTotal;

}

int main (int argc, char** argv)
{
  string sinputfile;
  string soutputfile;
  for (int i=1;i<argc;i++)
    {
      string sarg=argv[i];
      if (sarg == "-rf" || sarg == "--random-factor")
        {
          RandomFactor = atof (argv[++i]);
          continue;
        }
      if (! sinputfile.size ())
        sinputfile = argv[i];
      else if (! soutputfile.size ())
        soutputfile = argv[i];
      else
        cout << argv[i] << " ignored" << endl;
    }

  NoiseCoincidences nc(sinputfile.data ());
  unsigned NbHit=nc.Loop ();
  RunInfo ri (sinputfile.data ());
  double NBTrigg=ri.Loop ();

  cout << "norm " <<NBTrigg << endl;

  //unsigned NbNoCoin=NBTrigg-NbHit;
  //hcoincidences->Fill (0., NbNoCoin);

  //hcoincidences->Scale (13750./NBTrigg);
  //hcoincidencesMin->Scale (13750./NBTrigg);
  //hcoincidencesMax->Scale (13750./NBTrigg);
	
	//Vla: let's don't scale since for SN it is different
	hcoincidences->Scale (1/NBTrigg);   //so it is now in number of coincidences per decay * m^3 (generated volume)
	hcoincidencesMin->Scale (1/NBTrigg);
	hcoincidencesMax->Scale (1/NBTrigg);
	
	cout << "DEBUG: " << hcoincidences->GetBinContent(1) << " " << hcoincidences->GetBinContent(2) << " " << hcoincidences->GetBinContent(3) << endl;

  hpairlist->Scale (13750./NBTrigg);

  GetAngularCoincidences ();

  TFile output(soutputfile.data (),"recreate");
  output.Add(hcoincidences);
  output.Add(hcoincidencesMin);
  output.Add(hcoincidencesMax);
  output.Add(hpairlist);
  output.Add(hPhotonsInPM);
  output.Add(gAngularCoincidences);
  output.Add(hpairlistProj);

  output.Write ();
  output.Close ();

}
