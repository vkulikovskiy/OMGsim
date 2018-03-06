#include <WPropSteppingAction.hh>
#include <WPropAnalysisManager.hh>
#include <G4Step.hh>
#include <G4ThreeVector.hh>

#include <cstdlib>

#include <iostream>

using namespace std;

namespace km3net
{

  WPropSteppingAction::WPropSteppingAction()
  {
  }

  WPropSteppingAction::~WPropSteppingAction()
  {
  }


  void WPropSteppingAction::UserSteppingAction(const G4Step* Step)
  {

    if (!Step->GetPostStepPoint()->GetPhysicalVolume())
      return;

    const G4ThreeVector& pos=Step->GetPostStepPoint()->GetPosition ();

    if (pos.z () > 0.9*fMinZ)
      {
        if (pos.z () > fMinZ)
          Step->GetTrack()->SetTrackStatus(fStopAndKill);
        return;
      }

    G4String prestepName=Step->GetPreStepPoint()->GetPhysicalVolume()->GetName();
    G4String poststepName=Step->GetPostStepPoint()->GetPhysicalVolume()->GetName();



    const G4VProcess* theProc = Step->GetPreStepPoint()->GetProcessDefinedStep ();

    if (theProc && theProc->GetProcessName () == "KM3PetzoldScattering")
      {
        const KM3PetzoldScattering* thepet= (const KM3PetzoldScattering*)(theProc);
        WPropAnalysisManager::getIt ()->IncrementProcessWeight (thepet->GetLastProcess ());
      }

    if (prestepName != poststepName && prestepName.substr (1) == "Floor" && poststepName.substr (1) == "Floor")
      {
        int upfloor=Step->GetPostStepPoint()->GetPhysicalVolume()->GetCopyNo ();//atoi (poststepName.substr (0,1).data ());
        int downfloor=Step->GetPreStepPoint()->GetPhysicalVolume()->GetCopyNo ();//atoi (prestepName.substr (0,1).data ());
        WPropAnalysisManager::getIt ()->AddDetectedParticule ((KM3TrackInformation*)Step->GetTrack()->GetUserInformation(), (upfloor > downfloor) ? upfloor : downfloor);
      }

    return;



  }
}
