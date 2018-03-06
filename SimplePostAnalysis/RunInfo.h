//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Sep 30 14:45:45 2015 by ROOT version 5.34/19
// from TTree RunInfo/Infos about the simulation.
// found on file: Air.root
//////////////////////////////////////////////////////////

#ifndef RunInfo_h
#define RunInfo_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <string>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class RunInfo {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   ULong64_t       NBTrigg; //total number of K40 decay
   Float_t         TargetVolume; //volume in which the K40 are generated (NB the world is 10% bigger)

   // List of branches
   TBranch        *b_NbTrigg;   //!
   TBranch        *b_TargetVolume;   //!

   RunInfo(std::string);
   virtual ~RunInfo();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual double   Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef RunInfo_cxx
RunInfo::RunInfo(std::string filename) : fChain(0)
{
  TTree* tree=0;
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename.data ());
  if (!f || !f->IsOpen()) {
    f = new TFile(filename.data ());
  }
  f->GetObject("RunInfo",tree);
  Init(tree);
}

RunInfo::~RunInfo()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t RunInfo::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t RunInfo::LoadTree(Long64_t entry)
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

void RunInfo::Init(TTree *tree)
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

   fChain->SetBranchAddress("NBTrigg", &NBTrigg, &b_NbTrigg);
   fChain->SetBranchAddress("TargetVolume", &TargetVolume, &b_TargetVolume);
   Notify();
}

Bool_t RunInfo::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void RunInfo::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t RunInfo::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef RunInfo_cxx
