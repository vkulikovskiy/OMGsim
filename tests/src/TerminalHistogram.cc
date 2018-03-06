#include <TerminalHistogram.hh>

#include <iostream>
#include <iomanip>
#include <string>

#include <sys/ioctl.h>
#include <unistd.h>

using namespace std;

TerminalHistogram::TerminalHistogram (const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup) : TH1D (name, title, nbinsx, xlow, xup){}

void TerminalHistogram::Draw (Option_t* o)
{
  if (o)
    {
      //TODO manage some possible options, but remove the warning for now
      o=o;
    }
  struct winsize w;
  ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);

  int left_reserve = (to_string (TH1D::GetXaxis()->GetXmax ()).size () + 1);
  int tmp = (to_string (TH1D::GetXaxis()->GetXmin ()).size () + 1);
  if (tmp >left_reserve) left_reserve = tmp;

  string sreserve(left_reserve,' ');

	int prec = 5;
  float max_horizontal=(w.ws_col-left_reserve*1.2-prec-2)/TH1D::GetMaximum ();  //5 digits for number print - Vla

  cout << "\r" << setw(w.ws_col) << right << TH1D::GetMaximum () << left << sreserve << "|" << TH1D::GetMinimum () << endl;
  //cout << sreserve << TH1D::GetMinimum () << flush;m


  for (int i=1; i<=TH1D::GetNbinsX (); i++)
    {
      cout << sreserve << "|";
      for (int nbdashes=0; nbdashes<max_horizontal*TH1D::GetBinContent (i); nbdashes++)
        cout << "-";
      //cout << "\r" << flush << setw(left_reserve-1) << right <<TH1D::GetBinCenter (i);
			cout << setprecision(prec) << TH1D::GetBinContent (i) << "\r" << flush << setw(left_reserve-1) << right << TH1D::GetBinLowEdge (i);
      cout << endl;

    }

}
