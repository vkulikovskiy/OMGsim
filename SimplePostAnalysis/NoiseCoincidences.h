//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Sep 28 11:13:28 2015 by ROOT version 5.34/19
// from TTree Hit/Hits in a TTree
// found on file: Water.root
//////////////////////////////////////////////////////////

#ifndef NoiseCoincidences_h
#define NoiseCoincidences_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <string>
// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class NoiseCoincidences {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Declaration of leaf types
  Float_t         Vertex_VertexX;
  Float_t         Vertex_VertexY;
  Float_t         Vertex_VertexZ;
  Float_t         Vertex_VertexT;
  Float_t         Vertex_VertexP;
  Float_t         Vertex_VertexDirX;
  Float_t         Vertex_VertexDirY;
  Float_t         Vertex_VertexDirZ;
  Float_t         Vertex_VertexEnergy;
  Float_t         Vertex_OMHitX;
  Float_t         Vertex_OMHitY;
  Float_t         Vertex_OMHitZ;
  Float_t         Vertex_OMHitT;
  Float_t         Vertex_OMHitP;
  Float_t         Vertex_PhotoX;
  Float_t         Vertex_PhotoY;
  Float_t         Vertex_PhotoZ;
  Float_t         Vertex_PhotoT;
  Float_t         Vertex_PhotoP;
  Short_t         Vertex_OMID;
  Short_t         Vertex_PMID;
  Float_t         Hit_TotalTime;
  Float_t         Hit_Energy;
	//this one appeared in the latest KM3Sim versions:
	Bool_t          Polarity;
  UChar_t         Hit_Limit;
  ULong64_t       HitNumber;

  // List of branches
  TBranch        *b_Vertex;   //!
  TBranch        *b_Hit;   //!
  TBranch        *b_EventNb;   //!

  std::string    fFileName;

  NoiseCoincidences(std::string);
  virtual ~NoiseCoincidences();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual unsigned Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);

 private:
  void FillHistograms ();
  void InitVar ();
  void AddCoincidence (unsigned POS, unsigned CurrentPMMask);
  void AddEvent ();

  unsigned short fcoincidences[3];
  unsigned       ftphotons;
  float          fHitTimePerPM[32];
  unsigned       fPMMask[3];
};

#endif

#ifdef NoiseCoincidences_cxx
NoiseCoincidences::NoiseCoincidences(std::string filename) : fChain(0), fFileName(filename)
{
  TTree* tree=0;
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename.data ());
  if (!f || !f->IsOpen()) {
    f = new TFile(filename.data ());
  }
  f->GetObject("Hit",tree);
  Init(tree);
}

NoiseCoincidences::~NoiseCoincidences()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t NoiseCoincidences::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t NoiseCoincidences::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void NoiseCoincidences::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("Vertex", &Vertex_VertexX, &b_Vertex);
  fChain->SetBranchAddress("Hit", &Hit_TotalTime, &b_Hit);
  fChain->SetBranchAddress("HitNumber", &HitNumber, &b_EventNb);
  Notify();
}

Bool_t NoiseCoincidences::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void NoiseCoincidences::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}
Int_t NoiseCoincidences::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef NoiseCoincidences_cxx
