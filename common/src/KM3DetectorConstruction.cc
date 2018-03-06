//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
#include "KM3DetectorConstruction.hh"
#include "KM3DetectorMessenger.hh"
#include "G4UImanager.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"

#include "KM3Material.hh"
#include "KM3Material.hh"


#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UserLimits.hh"

#include "G4ios.hh"
#include <sstream>

#include <KM3Material.hh>
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4Hype.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4UserLimits.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4Ellipsoid.hh"
#include "G4NistManager.hh"

#include "G4ProductionCuts.hh"
#include "KM3PMTOpticalModel.hh"
#include "KM3InputDataReader.hh"

#include <sys/types.h>
#include <dirent.h>
#include <libconfig.h++>
#include <KM3ParametrizedVolume.hh>
#include <KM3PhotoTubeMgr.hh>
#include <KM3PassiveVolumeMgr.hh>
#include <KM3OpticalModuleMgr.hh>

using namespace std;
using namespace libconfig;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace km3net
{

  const string KM3DetectorConstruction::fOMBasePrefix("OM",2);
  G4double KM3DetectorConstruction::fTargetLength=500*m;

  KM3DetectorConstruction* KM3DetectorConstruction::fInstance=0;

  KM3DetectorConstruction* KM3DetectorConstruction::GetIt ()
  {
    if (!fInstance)
      fInstance=new KM3DetectorConstruction ();
    return (KM3DetectorConstruction*)fInstance;
  }


  KM3DetectorConstruction::KM3DetectorConstruction() {
    detectorMessenger = new KM3DetectorMessenger(this);

    fMaterialsManager = KM3Material::GetIt ();
    KM3PassiveVolumeMgr::GetIt ();
    KM3PhotoTubeMgr::GetIt ();
    KM3OpticalModuleMgr::GetIt ();

    ConfigFilePrefix="KM3Det";

    LoadDirectory ((string(PROJECT_SOURCE_DIR)+"/common/data/"));

    logicTarget = NULL;
    fTargetLength = 5.0 * m;
    //--------- Sizes of the principal geometrical components (solids)  ---------
    fSphereRmax = 13. / 2 * inch;
    FStoreyVolume = 0;

    fMaterialsManager->ListMaterial();
    TargetMater = fMaterialsManager->GetMaterial("AntaresWater");
  }

  G4VPhysicalVolume* KM3DetectorConstruction::PlaceThem (KM3ParametrizedVolume* theVolToPlace, PositionParameters& PosPar)
  {
    G4LogicalVolume* WhereToPlace=0;

    KM3PassiveVolumeMgr* theVolmgr = 0;
    if (fVolumeList.size ())
      {
        theVolmgr = KM3PassiveVolumeMgr::GetIt ();
      }
    G4RotationMatrix* theRot = new G4RotationMatrix;

		//PMT1 [0.9824,-1.5708]/
		//PMT2[1.2713,-3.1416]
		//PMT13 [1.8703,1.5708]
		double theta = 1.8703;
		double phi = 1.5708;

		G4RotationMatrix* theRotPMT3 = new G4RotationMatrix();
		theRotPMT3->rotateZ(-phi);     //this is done first, counter clockwise, rotates axes too
		theRotPMT3->rotateY(-theta);   //this will be done second
		theRotPMT3->invert();
		cout << "inverted rotation " << theRotPMT3->getPhi() << " " << theRotPMT3->getTheta() << " " << theRotPMT3->getPsi() << endl;

    theRot->set(PosPar.Phi, PosPar.Theta, PosPar.Psi);
    G4ThreeVector thePos(PosPar.X, PosPar.Y, PosPar.Z);

    G4LogicalVolume* logicVolume=theVolToPlace->GetLogicalVolume (PosPar.Name);
    if (!logicVolume)
      {
        cerr << "ERROR:KM3DetectorConstruction didn't place " << PosPar.Name <<
          " The corresponding conf file might not exist." << endl;
        return NULL;
      }
    if (theVolmgr && PosPar.MotherVolume.size ())
      {
        WhereToPlace=theVolmgr->GetLogicalVolume (PosPar.MotherVolume);
      }
    else
      {
        WhereToPlace = logicTarget;
      }
    if (!WhereToPlace)
      {
        cerr << "ERROR:KM3DetectorConstruction, volume "<< PosPar.MotherVolume<< " can't be found"<< endl;
        return NULL;
      }
    cout << "Placing " << logicVolume->GetName () << " " << logicVolume << " in " << WhereToPlace->GetName () << " " << WhereToPlace <<  " at " << thePos.x() << " " << thePos.y() << ' ' << thePos.z () << endl;
    return new G4PVPlacement(theRot, thePos, logicVolume, PosPar.Name.data (),
                      WhereToPlace, true, PosPar.ID);

  }
  G4VPhysicalVolume* KM3DetectorConstruction::Construct ()
  {

    //------------------------------
    // World
    //------------------------------
    solidWorld = new G4Box("world", fTargetLength * 1.5,
                           fTargetLength * 1.5 , fTargetLength * 1.5);
    logicWorld = new G4LogicalVolume(solidWorld,
                                     TargetMater, "World", 0, 0, 0);
    logicWorld->SetVisAttributes(G4VisAttributes::Invisible);
    //  Must place the World Physical volume unrotated at (0,0,0).
    //
    physiWorld = new G4PVPlacement(0,               // no rotation
                                   G4ThreeVector(), // at (0,0,0)
                                   logicWorld,      // its logical volume
                                   "World",         // its name
                                   0,               // its mother  volume
                                   false,           // no boolean operations
                                   0);              // no field specific to volume
    //------------------------------
    // Target
    //------------------------------
    solidTarget =  new G4Sphere("Target", 0, fTargetLength,
                                0, twopi, 0,
                                pi);
    /*new G4Box("Target", fTargetLength / 2., fTargetLength / 2.,
      fTargetLength / 2.);*/
    G4cout << "create logic target" << G4endl;
    logicTarget = new G4LogicalVolume(solidTarget,
                                      TargetMater, "Target", 0, 0, 0);
    G4VisAttributes* TargetVisAtt = new G4VisAttributes(
                                                        G4Colour(1.0, 1.0, 1.0));
    logicTarget->SetVisAttributes(TargetVisAtt);

    G4RotationMatrix *targetRot=new G4RotationMatrix(0,0,0);

    targetRot->set (FOrientationPhi, FOrientationTheta, FOrientationPsi);

    targetRegion = new G4Region("Target");


    for (auto& VolPar: fVolumeList [fDetectorTypeName])
      {
        fPlacedVolumeDict[VolPar.Name] = PlaceThem (KM3PassiveVolumeMgr::GetIt (), VolPar);
      }

    for (auto& VolPar: fVolumeList [fDetectorTypeName])
      {
        if (VolPar.BorderSurfaceMaterial.size () == 0 || VolPar.BorderSurfaceWith.size () ==0)
          continue;

        if (G4VPhysicalVolume* VolSecond = fPlacedVolumeDict[VolPar.BorderSurfaceWith])
          if (G4VPhysicalVolume* VolBase = fPlacedVolumeDict [VolPar.Name])
            if (G4OpticalSurface* the_surface = KM3Material::GetIt()->GetOpticalSurface (VolPar.BorderSurfaceMaterial))
              {
                new G4LogicalBorderSurface(VolPar.Name + VolPar.BorderSurfaceWith +
                                           VolPar.BorderSurfaceMaterial + "_border", VolBase, VolSecond,the_surface);
              }
      }
    for (auto& PMPar: fPMList [fDetectorTypeName])
      {
        PlaceThem ( KM3PhotoTubeMgr::GetIt (), PMPar);
      }

    for (auto& OMPar: fOMList [fDetectorTypeName])
      {
        PlaceThem ( KM3OpticalModuleMgr::GetIt (), OMPar);
      }



    physiTarget = new G4PVPlacement(targetRot, G4ThreeVector(), logicTarget, "Target",
                                    logicWorld, false, 0);

    return physiWorld;

  }

  void KM3DetectorConstruction::SetDataFilesSources (string Path) {
    KM3Material::GetIt () ->SetDataFilesSources (Path);
    KM3PassiveVolumeMgr::GetIt () ->SetDataFilesSources (Path);
    KM3PhotoTubeMgr::GetIt () ->SetDataFilesSources (Path);
    KM3OpticalModuleMgr::GetIt () ->SetDataFilesSources (Path);
    KM3ParametrizedVolume::SetDataFilesSources (Path);
  }

  void KM3DetectorConstruction::ReadThisConfiguration (Setting* cfg,
                                                       vector<PositionParameters>&  ConfList)
  {
    PositionParameters thePar;
    Setting
      &sX     = cfg->lookup("X"),
      &sY     = cfg->lookup("Y"),
      &sZ     = cfg->lookup("Z"),
      &sTheta = cfg->lookup("Theta"),
      &sPhi   = cfg->lookup("Phi"),
      &sName  = cfg->lookup("Name"),
      &sID    = cfg->lookup("ID");

    Setting *sPsi=0;
    if (cfg->exists("Psi"))
      sPsi = &(cfg->lookup("Psi"));


    vector <double>
      vX = GetScalarArrayParameter (sX),
      vY = GetScalarArrayParameter (sY),
      vZ = GetScalarArrayParameter (sZ),
      vTheta = GetScalarArrayParameter (sTheta),
      vPhi  = GetScalarArrayParameter (sPhi),
      vPsi;
    if (sPsi != 0)
        vPsi  = GetScalarArrayParameter (*sPsi);

    if (!vX.size ())
      thePar.X     =GetFloatingParameter (sX);
    if (!vY.size ())
      thePar.Y     =GetFloatingParameter (sY);
    if (!vZ.size ())
      thePar.Z     =GetFloatingParameter (sZ);
    if (!vTheta.size ())
      thePar.Theta =GetFloatingParameter (sTheta);
    if (!vPhi.size ())
      thePar.Phi   =GetFloatingParameter (sPhi);
    if (!vPsi.size () && sPsi != 0)
      thePar.Psi   =GetFloatingParameter (*sPsi);

    bool Multiplicity = vX.size() | vY.size() | vZ.size() | vTheta.size() | vPhi.size () | vPsi.size ();

    thePar.Name  =GetStringParameter   (sName);
    thePar.ID    =GetIntegerParameter  (sID);

    if (cfg->exists("MotherVolume"))
      thePar.MotherVolume=GetStringParameter (cfg->lookup("MotherVolume"));
    if (cfg->exists("BorderSurfaceWith"))
      thePar.BorderSurfaceWith=GetStringParameter (cfg->lookup("BorderSurfaceWith"));
    if (cfg->exists("BorderSurfaceMaterial"))
      thePar.BorderSurfaceMaterial=GetStringParameter (cfg->lookup("BorderSurfaceMaterial"));

    if (Multiplicity)
      {
        unsigned counter = 0;
        bool LetsDoIt=true;
        while (LetsDoIt)
          {
            LetsDoIt=false;
            if (counter < vX.size ())
              {
                thePar.X     = vX[counter];
                LetsDoIt=true;
              }

            if (counter < vY.size ())
              {
                thePar.Y     = vY[counter];
                LetsDoIt=true;
              }

            if (counter < vZ.size ())
              {
                thePar.Z     = vZ[counter];
                LetsDoIt=true;
              }

            if (counter < vTheta.size ())
              {
                thePar.Theta = vTheta[counter];
                LetsDoIt=true;
              }

            if (counter < vPhi.size ())
              {
                thePar.Phi   = vPhi[counter];
                LetsDoIt=true;
              }
            if (counter < vPsi.size ())
              {
                thePar.Psi   = vPsi[counter];
                LetsDoIt=true;
              }

            if (LetsDoIt)
              ConfList.push_back(thePar);

            thePar.ID++;

            counter++;
          }

        return;
      }

    ConfList.push_back(thePar);
  }

  void KM3DetectorConstruction::ReadConfiguration (string ConfigFile)
  {
    Config cfg;
    cfg.readFile (ConfigFile.data ());
    cfg.setAutoConvert (true);

    Setting& rootStg=cfg.getRoot();

    string Name;

    if (!rootStg.lookupValue ("Name", Name))
      return;

    int configFileSize=rootStg.getLength ();
    for (int cfIt=0; cfIt < configFileSize; cfIt++)
      {
        if (!rootStg[cfIt].isGroup ())
          continue;

        string currentName=rootStg[cfIt].getName ();
        if (currentName.find("Volume") != string::npos)
          ReadThisConfiguration (&rootStg[cfIt], fVolumeList[Name]);
        else if (currentName.find("OM") != string::npos)
          ReadThisConfiguration (&rootStg[cfIt], fOMList[Name]);
        else if (currentName.find("PM") != string::npos)
          ReadThisConfiguration (&rootStg[cfIt], fPMList[Name]);
      }
    return;
    /*keep the old method for now




    int i=1;
    string RootName="Volume1";
    while (cfg.exists (RootName.data ()))
      {
        ReadThisConfiguration (&cfg, RootName, fVolumeList[Name]);
        ostringstream oss;
        oss << ++i;
        RootName = string("Volume") + oss.str ();
      }

    i=1;
    RootName="PM1";
    while (cfg.exists (RootName.data ()))
      {
        ReadThisConfiguration (&cfg, RootName, fPMList[Name]);
        ostringstream oss;
        oss << ++i;
        RootName = string("PM") + oss.str ();
      }

    i=1;
    RootName="OM1";
    while (cfg.exists (RootName.data ()))
      {

        ReadThisConfiguration (&cfg, RootName, fOMList[Name]);
        ostringstream oss;
        oss << ++i;
        RootName = string("OM") + oss.str ();
      }
    */
  }



  void KM3DetectorConstruction::SetTargetLength(G4double value)
  {
    G4cout << "Target length has been set to " << value/m << " m" << G4endl;
    fTargetLength = value;
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  KM3DetectorConstruction::~KM3DetectorConstruction() {
    G4cout
      << "KM3DetectorConstruction: the Optical Module orientation angle was "
      << FOrientationTheta / deg << " deg" << G4endl;
    delete detectorMessenger;
    cout << "detectorMessenger deleted correctly" << endl;
  }

  void KM3DetectorConstruction::DontDraw (G4String DD)
  {
    fDD=fDD+DD;
  }

  void KM3DetectorConstruction::SetMeanCapRadius (G4double R)
  {
    fMeanCapRadius=R;
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


  void KM3DetectorConstruction::SetTargetMaterial(G4String materialName) {
    // search the material by its name
    cout << "---> " << materialName << endl;
    G4Material* pttoMaterial = KM3Material::GetMaterial(materialName);
    if (pttoMaterial) {
      TargetMater = pttoMaterial;
      if (logicTarget)
        logicTarget->SetMaterial(pttoMaterial);
      cout << "\n----> The target as been changed to "
             << fTargetLength / cm << " cm of " << materialName << G4endl;
    }
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void KM3DetectorConstruction::ProcessLogicalVolume(G4LogicalVolume* volume) {
    G4String vname = volume->GetName();
    G4String vmname = volume->GetMaterial()->GetName();
    //G4cout << "Volume: " << vname << "\t Material: " << vmname << G4endl;
    G4cout << "<Volume name=\"" << vname << "\" material=\"" << vmname << "\">"
           << G4endl;
    //G4cout << "Childs: " << volume->GetNoDaughters() << G4endl;
    for (int i = 0; i < volume->GetNoDaughters(); i++) {
      //G4cout << "Child #" << i << G4endl;
      ProcessLogicalVolume(volume->GetDaughter(i)->GetLogicalVolume());
    }
    G4cout << "</Volume>" << G4endl;
  }
}
