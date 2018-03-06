#ifndef WPropSteppingAction_h
#define WPropSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include <KM3SteppingAction.hh>

namespace km3net
{
  class WPropSteppingAction : public KM3SteppingAction
  {
  public:
    WPropSteppingAction();
    virtual ~WPropSteppingAction();

    virtual void UserSteppingAction(const G4Step*);

  public:
    double fMinZ=100;

  };
}

#endif
