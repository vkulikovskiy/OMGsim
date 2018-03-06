#include "KM3TrackInformation.hh"
#include "G4ios.hh"
#include <G4VPhysicalVolume.hh>
#include <G4NavigationHistory.hh>
#include <iostream>
using namespace std;

namespace km3net {
  G4Allocator<KM3TrackInformation> aTrackInformationAllocator;

  KM3TrackInformation::KM3TrackInformation()
  {
    originalTrackID = 0;
    originalParticleName = "";
    originalParticleType = "";
    originalPosition = G4ThreeVector(0.,0.,0.);
    originalMomentum = G4ThreeVector(0.,0.,0.);
    originalEnergy = 0.;
    originalTime = 0.;
  }

  KM3TrackInformation::KM3TrackInformation(const KM3TrackInformation* INFO){
    originalTrackID = INFO->GetOriginalTrackID();
    originalParticleName = INFO->GetOriginalParticleName();
    originalParticleType = INFO->GetOriginalParticleType();
    originalPosition = INFO->GetOriginalPosition();
    originalLocalPosition = INFO->GetOriginalPosition();
    originalPreVolume = INFO->GetOriginalPrePV();
    originalPostVolume = INFO->GetOriginalPostPV();
    originalMomentum = INFO->GetOriginalMomentum();
    originalEnergy = INFO->GetOriginalEnergy();
    originalTime = INFO->GetOriginalTime();
    originalPolarization = INFO->GetOriginalPolarization ();
  }

  KM3TrackInformation::KM3TrackInformation(const G4Track* aTrack)
  {

    originalTrackID = aTrack->GetTrackID();
    originalParticleName = aTrack->GetDefinition()->GetParticleName();
    originalParticleType = aTrack->GetDefinition()->GetParticleType();
    originalPosition = aTrack->GetPosition();


    if (aTrack->GetStep () != NULL) {
      originalLocalPosition =
        aTrack->GetStep ()->GetPostStepPoint ()->
        GetTouchable()->GetHistory ()->GetTopTransform().TransformPoint(originalPosition);
      originalPreVolume = aTrack->GetStep ()->GetPreStepPoint()->GetPhysicalVolume();
      originalPostVolume = aTrack->GetStep ()->GetPostStepPoint()->GetPhysicalVolume();
    }

    originalMomentum = aTrack->GetMomentum();
    originalEnergy = aTrack->GetTotalEnergy();
    originalTime = aTrack->GetGlobalTime();

    originalPolarization = aTrack->GetPolarization ();

  }

  KM3TrackInformation::~KM3TrackInformation(){;}

  void KM3TrackInformation::Print() const
  {
    G4cout
      << "name " << originalParticleName << " Original track ID " << originalTrackID
      << " at " << originalPosition << G4endl;
  }
}
