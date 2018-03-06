#ifndef __TerminalHistogram_hh__
#define __TerminalHistogram_hh__

#include <TH1D.h>

class TerminalHistogram : public TH1D
{
public:
  TerminalHistogram (const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup);
  virtual ~TerminalHistogram () {};
  void Draw (Option_t* o=0);

};





#endif
