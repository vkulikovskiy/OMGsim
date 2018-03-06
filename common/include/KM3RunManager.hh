#ifndef KM3RUNMANAGER_HH
#define KM3RUNMANAGER_HH

#include <G4RunManager.hh>
#include <KM3RunMessenger.hh>

namespace km3net
{
  class KM3RunManager : public G4RunManager
  {
  public:
    KM3RunManager ();
    virtual ~KM3RunManager ();

  private:
    KM3RunMessenger* fMessenger;

  };


}

#endif
