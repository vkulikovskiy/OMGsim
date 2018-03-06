#ifndef KM3TrackingAction_h
#define KM3TrackingAction_h

#include <G4UserTrackingAction.hh>


namespace km3net
{
  class KM3TrackingAction : public G4UserTrackingAction
  {
  public:
    KM3TrackingAction ();
    ~KM3TrackingAction ();

    virtual void PreUserTrackingAction(const G4Track* aTrack);
    virtual void PostUserTrackingAction(const G4Track* aTrack);
  };


}
#endif
