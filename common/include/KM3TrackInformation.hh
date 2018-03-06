#ifndef KM3TrackInformation_h
#define KM3TrackInformation_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Allocator.hh"
#include "G4VUserTrackInformation.hh"
#include <string>


class G4VPhysicalVolume;

namespace km3net{
  class KM3TrackInformation : public G4VUserTrackInformation
  {
  public:
    KM3TrackInformation();
    KM3TrackInformation(const G4Track* aTrack);
    KM3TrackInformation(const KM3TrackInformation* INFO);
    virtual ~KM3TrackInformation();

    inline void *operator new(size_t);
    inline void operator delete(void *aTrackInfo);
    inline int operator ==(const KM3TrackInformation& right) const
    {return (this==&right);}

    bool HasAMother () {return originalTrackID != 0;}

    void Print() const;

  private:
    G4int                 originalTrackID;
    std::string           originalParticleName;
    std::string           originalParticleType;
    G4Step                fLastStep;
    G4ThreeVector         originalPosition;
    G4ThreeVector         originalLocalPosition;
    G4ThreeVector         originalMomentum;
    G4ThreeVector         originalPolarization;
    G4double              originalEnergy;
    G4double              originalTime;

    G4VPhysicalVolume*    originalPreVolume = NULL;
    G4VPhysicalVolume*    originalPostVolume = NULL;

    const G4Track*              fTrack=0;
    const G4Step*               fStep=0;

  public:
    void SetTrack (const G4Track* track) {fTrack=track;}
    void SetStep  (const G4Step* step) {fStep=step;}

    inline G4int GetOriginalTrackID() const {return originalTrackID;}
    inline std::string GetOriginalParticleName() const {return originalParticleName;}
    inline std::string GetOriginalParticleType() const {return originalParticleType;}
    inline G4ThreeVector GetOriginalPosition() const {return originalPosition;}
    inline G4ThreeVector GetOriginalLocalPosition() const {return originalLocalPosition;}
    inline G4ThreeVector GetOriginalMomentum() const {return originalMomentum;}
    inline G4ThreeVector GetOriginalPolarization() const {return originalPolarization;}
    inline G4double GetOriginalEnergy() const {return originalEnergy;}
    inline G4double GetOriginalTime() const {return originalTime;}
    inline G4VPhysicalVolume* GetOriginalPrePV() const {return originalPreVolume;}
    inline G4VPhysicalVolume* GetOriginalPostPV() const {return originalPostVolume;}

    inline const G4Track* GetTrack () const {return fTrack;};
    inline const G4Step* GetStep () const {return fStep;};

  };

  extern G4Allocator<KM3TrackInformation> aTrackInformationAllocator;

  inline void* KM3TrackInformation::operator new(size_t)
  { void* aTrackInfo;
    aTrackInfo = (void*)aTrackInformationAllocator.MallocSingle();
    return aTrackInfo;
  }

  inline void KM3TrackInformation::operator delete(void *aTrackInfo)
  { aTrackInformationAllocator.FreeSingle((KM3TrackInformation*)aTrackInfo);}
}
#endif
