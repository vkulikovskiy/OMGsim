#include <KM3PhotoTubeMgr.hh>
#include <sys/types.h>
#include <dirent.h>
#include <libconfig.h++>
#include <KM3Material.hh>

#include <G4LogicalBorderSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4OpticalSurface.hh>

#include <G4VSolid.hh>
#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4Sphere.hh>
#include <G4Cons.hh>
#include <G4UnionSolid.hh>
#include <G4Ellipsoid.hh>
#include <G4SubtractionSolid.hh>
#include <G4IntersectionSolid.hh>

#include <G4LogicalVolume.hh>

#include <G4PVPlacement.hh>

#include <G4ThreeVector.hh>
#include <G4Transform3D.hh>
#include <G4VisAttributes.hh>

#include <G4UserLimits.hh>

#include <G4Step.hh>
#include <G4MaterialPropertiesTable.hh>
#include <G4TransportationManager.hh>
#include <G4RunManager.hh>

#include <KM3PMTOpticalModel.hh>
#include <KM3PhysicsList.hh>
#include <Randomize.hh>

using namespace std;
using namespace libconfig;


namespace km3net
{
  KM3PhotoTubeMgr::KM3ParametrizedVolume* KM3PhotoTubeMgr::fInstance = NULL;

  const string KM3PhotoTubeMgr::fPhotoRegionName  ("NoticeablePhotoCathodeRegion");
  const string KM3PhotoTubeMgr::fPhotoName        ("NoticeablePhotoCathode");
  const string KM3PhotoTubeMgr::fAbsorberName     ("NoticeableOpticalModuleAbsorber");
  const string KM3PhotoTubeMgr::fBentoName        ("NoticeableOpticalBentoGlass");
  const string KM3PhotoTubeMgr::fPMTContainerName ("NoticeablePMTContainer");

  KM3PhotoTubeMgr* KM3PhotoTubeMgr::GetIt ()
  {
    if (!fInstance)
      fInstance=new KM3PhotoTubeMgr ();
    return (KM3PhotoTubeMgr*)fInstance;
  }

  KM3PhotoTubeMgr::KM3PhotoTubeMgr ()
  {
    ConfigFilePrefix="KM3PMT";
    LoadDirectory ((string(PROJECT_SOURCE_DIR)+"/common/data/"));
  }

  KM3ParametrizedVolume::KM3Volume* KM3PhotoTubeMgr::ConstructTheVolume (VolumeParameters* PTP)
  {
    return new KM3PhotoTube ((PhotoTubeParameters*)PTP);
  }

  KM3PhotoTubeMgr::~KM3PhotoTubeMgr ()
  {
  }

  void KM3PhotoTubeMgr::ReadConfiguration (std::string ConfigFile)
  {
    Config cfg;
    cfg.setAutoConvert (true);
    cfg.readFile (ConfigFile.data ());
    string name;
    if (cfg.lookupValue("Name", name))
      {
        fParametersDict[name]=new PhotoTubeParameters;
      }else
      {
        cerr << "no name parameter. Will probably crash."<<endl;
        return;
      }

    PhotoTubeParameters *thePTP=(PhotoTubeParameters*)(fParametersDict[name]);

    if (!ReadBaseParameters (&cfg, thePTP))
      {
        cerr << "ReadBaseParameters: error: misse base parameters for " << name << endl;
        return;
      }

    thePTP->IsOk =
      cfg.exists("RadiusSphere") &&
      cfg.exists("BulbRadiusMin") &&
      cfg.exists("BulbRadiusMax") &&
      cfg.exists("BulbHeight") &&
      cfg.exists("ConeHeight") &&
      cfg.exists("TubeRadius") &&
      cfg.exists("GlassThickness") &&
      cfg.exists("TotalPMTHeight") &&
      cfg.exists("DinodsPos") &&
      cfg.exists("DinodsRadius") &&
      cfg.exists("DinodsHeight");

    thePTP->RadiusSphere = GetFloatingParameter
      (cfg.lookup ("RadiusSphere"));
    thePTP->BulbRadiusMin = GetFloatingParameter
      (cfg.lookup ("BulbRadiusMin"));
    thePTP->BulbRadiusMax = GetFloatingParameter
      (cfg.lookup ("BulbRadiusMax"));
    thePTP->BulbHeight = GetFloatingParameter
      (cfg.lookup ("BulbHeight"));
    thePTP->ConeHeight = GetFloatingParameter
      (cfg.lookup ("ConeHeight"));
    thePTP->TubeRadius = GetFloatingParameter
      (cfg.lookup ("TubeRadius"));
    thePTP->GlassThickness = GetFloatingParameter
      (cfg.lookup ("GlassThickness"));
    thePTP->TotalPMTHeight = GetFloatingParameter
      (cfg.lookup ("TotalPMTHeight"));
    thePTP->DinodsPos = GetFloatingParameter
      (cfg.lookup ("DinodsPos"));
    thePTP->DinodsRadius = GetFloatingParameter
      (cfg.lookup ("DinodsRadius"));
    thePTP->DinodsHeight = GetFloatingParameter
      (cfg.lookup ("DinodsHeight"));


    if (!thePTP->IsOk)
      {
        clog << "KM3PhotoTubeMgr:WARNING in file: "<< ConfigFile <<  ": it does not contain all the parameters. The produced PM volume could be NULL, leading to an eventual segmentation fault." <<endl <<
          "Prefix "        << thePTP->Prefix<< endl <<
          "RadiusSphere "  << thePTP->RadiusSphere<< endl <<
          "BulbRadiusMin " << thePTP->BulbRadiusMin<< endl <<
          "BulbRadiusMax " << thePTP->BulbRadiusMax<< endl <<
          "BulbHeight "    << thePTP->BulbHeight<< endl <<
          "ConeHeight "    << thePTP->ConeHeight<< endl <<
          "TubeRadius "  << thePTP->TubeRadius<< endl <<
          "GlassThickness "<< thePTP->GlassThickness<< endl <<
          "TotalPMTHeight "<< thePTP->TotalPMTHeight<< endl <<
          "DinodsPos "     << thePTP->DinodsPos<< endl <<
          "DinodsRadius "  << thePTP->DinodsRadius<< endl <<
          "DinodsHeight "  << thePTP->DinodsHeight<< endl;
      }

    cfg.lookupValue("PhotoCathodeMaterial", thePTP->PhotoCathodeMaterial);
    cfg.lookupValue("InternalReflectorMaterial", thePTP->InternalReflectorMaterial);
    cfg.lookupValue("DinodsSurface", thePTP->DinodsSurface);
    cfg.lookupValue("CapDinodsSurface", thePTP->CapDinodsSurface);
    cfg.lookupValue("GridDinodsSurface", thePTP->GridDinodsSurface);
    cfg.lookupValue("PMTSurface"       , thePTP->PMTSurface);

		cfg.lookupValue("GlassMaterial", thePTP->GlassMaterial);  //VLA
		cout << "Glass Value: " << thePTP->GlassMaterial << endl;

    thePTP->CalculateParameters ();
  }


  KM3PhotoTubeMgr::KM3PhotoTube::KM3PhotoTube (PhotoTubeParameters* PTP): KM3ParametrizedVolume::KM3Volume(PTP)
  {
    cout << "creating it " << PTP->Prefix << endl;
  }

  KM3PhotoTubeMgr::KM3PhotoTube::~KM3PhotoTube ()
  {
  }

  bool KM3PhotoTubeMgr::IsInAPhotochatode ()
  {
    if (G4TransportationManager::GetTransportationManager()->
        GetNavigatorForTracking()->CreateTouchableHistory ()->GetHistoryDepth () < 1)
      return false;

    const G4LogicalVolume* theLogVol = G4TransportationManager::GetTransportationManager()->
      GetNavigatorForTracking()->CreateTouchableHistory ()->GetVolume (1)->GetLogicalVolume ();

    if (fLogicalVolDict.count (theLogVol))
      {
        fCurrentPMTName = fLogicalVolDict[theLogVol];
        return true;
      }
    return false;

  }

  unsigned char KM3PhotoTubeMgr::IsDetected (G4Step* STEP)
  {

    G4Track* theTrack = STEP->GetTrack();
    if (theTrack->GetDefinition()->GetParticleType() != "opticalphoton" || theTrack->GetTrackStatus() != fStopAndKill)
      return false;

    bool isInPh = IsInAPhotochatode ();
    if (!isInPh)
      return false;

    G4bool isOnBoundary =
      (STEP->GetPostStepPoint()->GetStepStatus() == fGeomBoundary);
    if (!isOnBoundary)
      return false;

    return ((KM3PhotoTube*)(fVolumeDict[fCurrentPMTName]))-> CalculateDetection (STEP);
  }

  unsigned char KM3PhotoTubeMgr::KM3PhotoTube::CheckDetection (G4Step* STEP)
  {
    G4VPhysicalVolume* thePrePV=STEP->GetPreStepPoint()->GetPhysicalVolume();
      G4VPhysicalVolume* thePostPV=STEP->GetPostStepPoint()->GetPhysicalVolume();

      G4LogicalSurface* Surface = G4LogicalBorderSurface::GetSurface(thePrePV, thePostPV);
      if (!Surface)
        Surface = G4LogicalBorderSurface::GetSurface(thePostPV, thePrePV);
      if (!Surface)
        Surface = G4LogicalSkinSurface::GetSurface(thePostPV->GetLogicalVolume ());
    if (!Surface)
      Surface = G4LogicalSkinSurface::GetSurface(thePrePV->GetLogicalVolume ());

    if (!Surface ||
        !LoadMaterialProps (dynamic_cast <G4OpticalSurface*> (Surface->GetSurfaceProperty())
                            ->GetMaterialPropertiesTable()) )
      return 0;

    unsigned char detectedGamma = CalculateDetection (STEP);

    return detectedGamma;
  }


  unsigned char KM3PhotoTubeMgr::KM3PhotoTube::CalculateDetection(const G4Step* STEP)
  {
    KM3PMTOpticalModel* theModel = ((KM3PhysicsList*)(G4RunManager::GetRunManager()->GetUserPhysicsList ()))->GetOpticalModel ();

    G4ThreeVector globalPoint = STEP->GetPreStepPoint()->GetPosition ();
    G4ThreeVector localPoint=STEP->GetPostStepPoint()->GetTouchable()->GetHistory()->GetTopTransform().TransformPoint(globalPoint);
    double Energy = STEP->GetPreStepPoint()->GetKineticEnergy();
    unsigned char detectedGamma = 0;
    double Efficiency=_efficiency_photocathode->Value (Energy);
    double eta_E = Efficiency/theModel->CalculateNormalCoefficients ();//\/OpticalModuleMaker::GetCurrentAbsorption ();

    double teta = localPoint.theta ();//acos(costeta);

    double eta_teta=_angular_efficiency->Value(teta/3*2.1);
    double eta_teta_err=_angular_efficiency_err->Value(teta/3*2.1);


    double rndm = G4UniformRand();
    if(rndm < eta_E)
      {
        rndm=G4UniformRand();
        if(rndm < eta_teta*(1.-eta_teta_err))
          {
            detectedGamma = 3;
          }
        else if(rndm < eta_teta)
          {
            detectedGamma = 2;
          }
        else if (rndm < eta_teta*(1+eta_teta_err))
          {
            detectedGamma = 1;
          }
        else
          {
            detectedGamma = 4;
          }

      }

    return detectedGamma;
  }

  bool KM3PhotoTubeMgr::KM3PhotoTube::LoadMaterialProps (G4MaterialPropertiesTable* aMaterialPropertiesTable)
  {
    if (!aMaterialPropertiesTable) {
      return false;
    }

    _angular_efficiency_err= aMaterialPropertiesTable->GetProperty("ANGULAR_EFFICIENCY_ERROR");
    if (_angular_efficiency == NULL)
      {
        return false;
      }

    _angular_efficiency= aMaterialPropertiesTable->GetProperty("ANGULAR_EFFICIENCY");
    if (_angular_efficiency == NULL)
      {
        return false;
      }

    _efficiency_photocathode= aMaterialPropertiesTable->GetProperty("EFFICIENCY");
    if (_efficiency_photocathode == NULL)
      {
        return false;
      }

    G4MaterialPropertyVector* thickness_photocathode= aMaterialPropertiesTable->GetProperty("THICKNESS");
    if (thickness_photocathode == NULL)
      {
        return false;
      }

    G4MaterialPropertyVector* kindex_photocathode= aMaterialPropertiesTable->GetProperty("KINDEX");
    if (kindex_photocathode == NULL)
      {
        return false;
      }

    G4MaterialPropertyVector* rindex_photocathode= aMaterialPropertiesTable->GetProperty("RINDEX");
    if (rindex_photocathode == NULL)
      {
        return false;
      }

    return true;

  }

  G4LogicalVolume* KM3PhotoTubeMgr::KM3PhotoTube::Construct(const VolumeParameters* PTP)
  {
    if (fTheLogicalVolume)
      return fTheLogicalVolume;

    if (PTP)
      fTheParameters = new PhotoTubeParameters (*((PhotoTubeParameters*)PTP));

    const PhotoTubeParameters &thePTP=*((PhotoTubeParameters*)fTheParameters);

    KM3Material* materialsManager = KM3Material::GetIt ();

    //create the fullfilled volumes and unit it to create the container


    //the photocathode
    G4Sphere* solidPMTPhotocathodeExt =
      new G4Sphere(thePTP.Prefix + "PMTGlassExt", 0.,
                   thePTP.RadiusSphere, 0, 2 * pi, 0, thePTP.PhotoAngle);
    G4Ellipsoid* solidPMTEllipsPhotoExt =
      new G4Ellipsoid(thePTP.Prefix + "PhotoEllipseExt", thePTP.BulbRadiusMax, thePTP.BulbRadiusMax,
                      thePTP.Bulbr, 0, thePTP.BulbHeight/2); //hack to draw in geant4.6.3. The G4Union avoid overlaps.
    G4BooleanSolid* solidFullPhotocathode =
      new G4UnionSolid(thePTP.Prefix + fPhotoRegionName,
                       solidPMTPhotocathodeExt, solidPMTEllipsPhotoExt,
                       NULL, G4ThreeVector (0,0,thePTP.RadiusSphere-thePTP.PhotoHeight - thePTP.BulbHeight/2));

    //the Glass
    G4Ellipsoid* solidPMTEllipsGlassExt =
      new G4Ellipsoid(thePTP.Prefix + "GlassEllipseExt", thePTP.BulbRadiusMax, thePTP.BulbRadiusMax,
                      thePTP.Bulbr, -thePTP.BulbHeight/2, 0.);
    G4Cons* solidConeExt = new G4Cons(thePTP.Prefix + "ConeExt",
                                      0., thePTP.TubeRadius,
                                      0., thePTP.BulbRadiusMin,
                                      thePTP.ConeHeight/2,
                                      0., 2 * pi);
    G4Tubs* solidTubeExt = new G4Tubs(thePTP.Prefix + "TubeExt",
                                      0., thePTP.TubeRadius,
                                      (thePTP.TotalPMTHeight - thePTP.PhotoHeight - thePTP.BulbHeight - thePTP.ConeHeight)/2,
                                      0, 2 * pi);

		cout << "Tube height is: " << (thePTP.TotalPMTHeight - thePTP.PhotoHeight - thePTP.BulbHeight - thePTP.ConeHeight) << endl;

    G4BooleanSolid* solidFullGlass = new G4UnionSolid (thePTP.Prefix + "FullGlassStep1",
                                                       solidPMTEllipsGlassExt, solidConeExt,
                                                       NULL, G4ThreeVector (0,0,-thePTP.BulbHeight/2-thePTP.ConeHeight/2));
    solidFullGlass = new G4UnionSolid (thePTP.Prefix + "FullGlass",
                                       solidFullGlass, solidTubeExt,
                                       NULL, G4ThreeVector (0,0, - thePTP.TotalPMTHeight/2 + ( thePTP.PhotoHeight - thePTP.ConeHeight)/2));


    //the container
    G4BooleanSolid* solidPMTContainer = new G4UnionSolid ( thePTP.Prefix + "PMTContainer",
                                                           solidFullPhotocathode,solidFullGlass,
                                                           NULL, G4ThreeVector (0,0,thePTP.RadiusSphere-thePTP.PhotoHeight - thePTP.BulbHeight/2));


    //then create the logical volume to set everything in
    G4LogicalVolume* logicPMTContainer = new G4LogicalVolume(solidPMTContainer,
                                                             materialsManager->GetMaterial("Air"), thePTP.Prefix + fPMTContainerName);

    //create the internal volumes and add it to substract to the ext volumes
    G4Sphere* solidPMTGlassInt =
      new G4Sphere(thePTP.Prefix + "PMTGlassInt", 0.,
                   thePTP.RadiusSphere - thePTP.GlassThickness, 0, 2 * pi, 0, thePTP.PhotoAngleIn);
    G4Ellipsoid* solidPMTEllipsInt =
      new G4Ellipsoid(thePTP.Prefix + "EllipseInt", thePTP.BulbRIn, thePTP.BulbRIn,
                      thePTP.BulbrIn, -thePTP.BulbHeight/2, thePTP.BulbHeight/2);
    G4Cons* solidConeInt =
      new G4Cons(thePTP.Prefix + "ConeiNT",
                 0., thePTP.TubeRadius - thePTP.GlassThickness,
                 0., thePTP.BulbRadiusMin - (thePTP.GlassThickness/sin (thePTP.PhotoAngle)),
                 thePTP.ConeHeight/2,
                 0., 2 * pi);
    G4Tubs* solidTubeInt =
      new G4Tubs(thePTP.Prefix + "TubeiNT",
                 0., thePTP.TubeRadius - thePTP.GlassThickness,
                 (thePTP.TotalPMTHeight - thePTP.PhotoHeight - thePTP.BulbHeight - thePTP.ConeHeight)/2 - thePTP.GlassThickness,
                 0, 2 * pi);


    G4BooleanSolid* SolidInternalPMT =
      new G4UnionSolid (thePTP.Prefix + "InternalPMTStep1",
                        solidPMTGlassInt, solidPMTEllipsInt,
                        NULL, G4ThreeVector (0,0,thePTP.RadiusSphere-thePTP.PhotoHeight - thePTP.BulbHeight/2));
    SolidInternalPMT =
      new G4UnionSolid (thePTP.Prefix + "InternalPMTStep2",
                        SolidInternalPMT, solidConeInt,
                        NULL, G4ThreeVector (0,0,thePTP.RadiusSphere-thePTP.PhotoHeight-thePTP.BulbHeight-thePTP.ConeHeight/2));
    SolidInternalPMT =
      new G4UnionSolid (thePTP.Prefix + "InternalPMT",
                        SolidInternalPMT, solidTubeInt,
                        NULL, G4ThreeVector (0,0,thePTP.RadiusSphere - thePTP.TotalPMTHeight/2 + ( - thePTP.PhotoHeight - thePTP.BulbHeight - thePTP.ConeHeight)/2 + thePTP.GlassThickness));

    //now the substractions

    G4BooleanSolid* solidPhotocathode =
      new G4SubtractionSolid(thePTP.Prefix + fPhotoName,
                             solidFullPhotocathode,SolidInternalPMT,
                             NULL,
                             G4ThreeVector(0,0,0));

    G4BooleanSolid* solidGlass =
      new G4SubtractionSolid(thePTP.Prefix + "Glass",
                             solidFullGlass,SolidInternalPMT,
                             NULL,
                             G4ThreeVector(0,0,-thePTP.RadiusSphere+thePTP.PhotoHeight + thePTP.BulbHeight/2));

    //create dinods solids and grid
    G4Tubs* solidDinods = new G4Tubs(thePTP.Prefix + "Dinods", thePTP.DinodsRadius - 1*mm, thePTP.DinodsRadius,
                                     thePTP.DinodsHeight, 0, 2 * pi);
    G4Tubs* solidCapDinods = new G4Tubs(thePTP.Prefix + "CapDinods", 0, thePTP.DinodsRadius - 1*mm,
                                        0.5*mm, 0, 2 * pi);
    G4Tubs* solidGridDinods = new G4Tubs(thePTP.Prefix + "GridDinods", // /2.6 for km3net...giusto
                                         0., (thePTP.DinodsRadius - 1*mm)/2. * 1.1222, //1.1222 = mean radius of a square, usable for the scans
                                         0.1*mm,
                                         0, 2 * pi);

    G4Tubs* solidBackCapDinods = new G4Tubs(thePTP.Prefix + "BackCapDinods", 0, thePTP.DinodsRadius,
                                     0.5*mm, 0, 2 * pi);

    //create the three logic volumes and place them

    G4LogicalVolume* logicPhotocathode =
      new G4LogicalVolume(solidPhotocathode,
                          materialsManager->GetMaterial(thePTP.GlassMaterial), solidPhotocathode->GetName ());
    G4LogicalVolume* logicGlass =
      new G4LogicalVolume(solidGlass,
                          materialsManager->GetMaterial(thePTP.GlassMaterial), thePTP.Prefix + thePTP.GlassMaterial);
    G4LogicalVolume* logicVaccum =
      new G4LogicalVolume(SolidInternalPMT,
                          materialsManager->GetMaterial("PMTVacuum"), thePTP.Prefix + "VaccumOfPmt");
													//materialsManager->GetMaterial("PMTVacuum"), thePTP.Prefix + "NOREFLEXVaccumInPmt");
    G4LogicalVolume* logicDinods =
      new G4LogicalVolume(solidDinods,
                          materialsManager->GetMaterial("AbsorberPVC"), thePTP.Prefix + "Dinods");
    G4LogicalVolume* logicCapDinods =
      new G4LogicalVolume(solidCapDinods,
                          materialsManager->GetMaterial("AbsorberPVC"), thePTP.Prefix + "CapDinods");
    G4LogicalVolume* logicGridDinods =
      new G4LogicalVolume(solidGridDinods,
                          materialsManager->GetMaterial("AbsorberPVC"), thePTP.Prefix + "GridDinods");

    G4LogicalVolume* logicBackCapDinods =
      new G4LogicalVolume(solidBackCapDinods,
                          materialsManager->GetMaterial("AbsorberPVC"), thePTP.Prefix + "BackCapDinods");


    //Set some cool colors
    logicVaccum->SetVisAttributes(ftransparent);
    logicPMTContainer->SetVisAttributes(ftransparent);
    logicPhotocathode->SetVisAttributes(fgold);
    logicGlass->SetVisAttributes(fsilver);
    logicDinods->SetVisAttributes(fsilver);
    logicCapDinods->SetVisAttributes(fgold);

    //Then create all of the surfaces, reflectives and photocathode

    G4OpticalSurface* our_Mirror_opsurf = KM3Material::GetIt ()->GetOpticalSurface (thePTP.InternalReflectorMaterial);

    G4OpticalSurface* our_Dinods_opsurf = KM3Material::GetIt ()->GetOpticalSurface (thePTP.DinodsSurface);

    G4OpticalSurface* our_CapDinods_opsurf = KM3Material::GetIt ()->GetOpticalSurface (thePTP.CapDinodsSurface);

    G4OpticalSurface* our_GridDinods_opsurf = KM3Material::GetIt ()->GetOpticalSurface (thePTP.GridDinodsSurface);

		G4OpticalSurface* our_BackCapDinods_opsurf = KM3Material::GetIt ()->GetOpticalSurface (thePTP.GridDinodsSurface);  //absorbing back cap

    G4OpticalSurface* Photocathode_opsurf = new G4OpticalSurface(thePTP.Prefix + "Photocathode_opsurf");
    Photocathode_opsurf->SetType(dielectric_dielectric);
    Photocathode_opsurf->SetMaterialPropertiesTable(G4Material::GetMaterial(thePTP.PhotoCathodeMaterial)->GetMaterialPropertiesTable());


    //Apply to the reflective, photcathode and dinods elements

    /*G4LogicalSkinSurface* MirrorGlass =
      new G4LogicalSkinSurface(thePTP.Prefix + "MirrorGlass",
      logicGlass,
      our_Mirror_opsurf);//*/
    new G4LogicalSkinSurface(thePTP.Prefix + "MirrorDinods",
                             logicDinods,
                             our_Dinods_opsurf);
    new G4LogicalSkinSurface(thePTP.Prefix + "MirrorCapDinods",
                             logicCapDinods,
                             our_CapDinods_opsurf);
    new G4LogicalSkinSurface(thePTP.Prefix + "MirrorGridDinods",
                             logicGridDinods,
                             our_GridDinods_opsurf);

    new G4LogicalSkinSurface(thePTP.Prefix + "MirrorBackCapDinods",
                             logicBackCapDinods,
                             our_BackCapDinods_opsurf);
    if (thePTP.PMTSurface.size ())
      {
        new G4LogicalSkinSurface("PMTSurface", logicPhotocathode, KM3Material::GetIt()->GetOpticalSurface (thePTP.PMTSurface));
      }


    G4VPhysicalVolume* physiPhotocathode =
      new G4PVPlacement(0,
                        G4ThreeVector(0., 0., 0.),
                        logicPhotocathode, logicPhotocathode->GetName (),
                        logicPMTContainer,
                        false, 0, !thePTP.CheckOverlaps);

    //decrepated
    //fPhotoCathodeList.push_back (physiPhotocathode);

    G4VPhysicalVolume* physiGlass =
      new G4PVPlacement(0,
                        G4ThreeVector(0., 0., thePTP.RadiusSphere-thePTP.PhotoHeight - thePTP.BulbHeight/2),
                        logicGlass, thePTP.Prefix + "Glass",
                        logicPMTContainer,
                        false, 0, !thePTP.CheckOverlaps);

    G4VPhysicalVolume* physiVaccum =
      new G4PVPlacement(0,
                        G4ThreeVector(0., 0., 0.),
                        logicVaccum, thePTP.Prefix + "VaccumInPmt",
												//logicVaccum, thePTP.Prefix + "NOREFLEXVaccumInPmt",
                        logicPMTContainer,
                        false, 0, !thePTP.CheckOverlaps);
    new G4PVPlacement(0,
                      G4ThreeVector(0., 0., thePTP.DinodsPos),
                      logicDinods, thePTP.Prefix + "Dinods",
                      logicVaccum,
                      false, 0, !thePTP.CheckOverlaps);
    new G4PVPlacement(0,
                      G4ThreeVector(0., 0., thePTP.DinodsPos+thePTP.DinodsHeight+0.5*mm), //to be on top of the dynode tube (0.5*mm is the cap size)
                      logicCapDinods, thePTP.Prefix + "CapDinods",
                      logicVaccum,
                      false, 0, !thePTP.CheckOverlaps);
    new G4PVPlacement(0,
                      G4ThreeVector(0., 0., thePTP.DinodsPos+thePTP.DinodsHeight+0.5*mm*2 + 0.1*mm), //to be on top of the cap
                      logicGridDinods, thePTP.Prefix + "GridDinods",
                      logicVaccum,
                      false, 0, !thePTP.CheckOverlaps);
    new G4PVPlacement(0,
                      G4ThreeVector(0., 0., thePTP.DinodsPos-thePTP.DinodsHeight-0.5*mm), //to be on top of the cap
                      logicBackCapDinods, thePTP.Prefix + "BackCapDinods",
                      logicVaccum,
                      false, 0, !thePTP.CheckOverlaps);

    new G4LogicalBorderSurface(thePTP.Prefix + "photocathodeGlass", physiPhotocathode, physiVaccum, Photocathode_opsurf);
    new G4LogicalBorderSurface(thePTP.Prefix + "photocathodeGlass", physiVaccum, physiPhotocathode, Photocathode_opsurf);//*/
    new G4LogicalBorderSurface(thePTP.Prefix + "MirrorGlass", physiVaccum, physiGlass, our_Mirror_opsurf);//*/
    new G4LogicalBorderSurface(thePTP.Prefix + "MirrorGlass", physiGlass, physiVaccum, our_Mirror_opsurf);//*/

    //create photoregion to apply homemade physics in photocathode

    fTheLogicalVolume = logicPMTContainer;
    return logicPMTContainer;
  }

}
