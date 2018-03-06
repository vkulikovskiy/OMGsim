#include <KM3OpticalModuleMgr.hh>
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
  KM3OpticalModuleMgr::KM3ParametrizedVolume* KM3OpticalModuleMgr::fInstance = NULL;

  const string KM3OpticalModuleMgr::fAbsorberName     ("NoticeableOpticalModuleAbsorber");
  const string KM3OpticalModuleMgr::fBentoName        ("NoticeableOpticalBentoGlass");
  const string KM3OpticalModuleMgr::fPMTContainerName ("NoticeablePMTContainer");

  KM3OpticalModuleMgr* KM3OpticalModuleMgr::GetIt ()
  {
    if (!fInstance)
      fInstance=new KM3OpticalModuleMgr ();
    return (KM3OpticalModuleMgr*)fInstance;
  }

  KM3OpticalModuleMgr::KM3OpticalModuleMgr ()
  {
    ConfigFilePrefix="KM3OM";

    LoadDirectory ((string(PROJECT_SOURCE_DIR)+"/common/data/"));
  }

  KM3OpticalModuleMgr::KM3Volume* KM3OpticalModuleMgr::ConstructTheVolume (VolumeParameters* PTP)
  {
    return new KM3OpticalModule ((OpticalModuleParameters*)PTP);
  }

  void KM3OpticalModuleMgr::ReadConfiguration (string ConfigFile)
  {
    Config cfg;
    cfg.setAutoConvert (true);
    cfg.readFile (ConfigFile.data ());
    string name;

    if (cfg.lookupValue("Name", name))
      {

        fParametersDict[name]=new OpticalModuleParameters;
      }else
      {
        cerr << "KM3OpticalModuleMgr  no Name parameter for " << ConfigFile<< endl;
        return;
      }
    OpticalModuleParameters *thePTP=(OpticalModuleParameters*)(fParametersDict[name]);

    if (!ReadBaseParameters (&cfg, thePTP))
      {
        cerr << "ReadBaseParameters: error: misse base parameters for " << name << endl;
        return;
      }


    cfg.lookupValue("DetectionUnitType", thePTP->DetectionUnitType);
    thePTP->IsOk =
      cfg.exists("SphereRmax") &&
      cfg.exists("GlassThickness") &&
      cfg.exists("DistSpherePMT") &&
      cfg.exists("PhotoTubeName");

    if (!thePTP->IsOk)
      {
        clog << "KM3OpticalModuleMgr:WARNING in file: "<< ConfigFile <<  ": it does not contain all the parameters. The produced PM volume could be NULL, leading to an eventual segmentation fault." <<endl;
      }


    thePTP->SphereRmax = GetFloatingParameter
      (cfg.lookup ("SphereRmax"));
    thePTP->GlassThickness = GetFloatingParameter
      (cfg.lookup ("GlassThickness"));
    thePTP->DistSpherePMT = GetFloatingParameter
      (cfg.lookup ("DistSpherePMT"));
    cfg.lookupValue("PhotoTubeName"    , thePTP->PhotoTubeName    );

    if (thePTP->DetectionUnitType == "OM")
      {
        thePTP->IsOk &=
          cfg.exists("AbsorberOvertake") &&
          cfg.exists("GelOvertake"     );

        thePTP->AbsorberOvertake = GetFloatingParameter
          (cfg.lookup ("AbsorberOvertake"));
        thePTP->GelOvertake = GetFloatingParameter
          (cfg.lookup ("GelOvertake"));
      }
    else if (thePTP->DetectionUnitType == "DOM")
      {
        cfg.lookupValue("AreThereReflectors", thePTP->AreThereReflectors);

	thePTP->IsOk &=
          cfg.exists("AreThereReflectors") &&
          (!thePTP->AreThereReflectors ||
           (cfg.exists("ReflectorAngle") &&
            cfg.exists("ReflectorRadiusMin") &&
            cfg.exists("ReflectorDistGlass") &&
            cfg.exists("ReflectorThickness")));

        if (thePTP->AreThereReflectors)
          {
            thePTP->ReflectorAngle = GetFloatingParameter
              (cfg.lookup ("ReflectorAngle"));
            thePTP->ReflectorRadiusMin = GetFloatingParameter
              (cfg.lookup ("ReflectorRadiusMin"));
            thePTP->ReflectorDistGlass = GetFloatingParameter
              (cfg.lookup ("ReflectorDistGlass"));
            thePTP->ReflectorThickness = GetFloatingParameter
              (cfg.lookup ("ReflectorThickness"));;
          }

        //optional values
        cfg.lookupValue("MushroomMaterial", thePTP->MushroomMaterial);
        if (cfg.exists ("MushroomAngularSize"))
          thePTP->MushroomAngularSize = GetFloatingParameter
            (cfg.lookup ("MushroomAngularSize"));

        thePTP->ReflectorAngle*=deg; //todo change it for radian

        Setting &angles=cfg.lookup("PMTPositions");
        for (int row=0; row < angles.getLength (); row++)
          {
            thePTP->PMTsAngles.push_back( {float(angles[row][0]), float(angles[row][1])});
          }
      }

    cfg.lookupValue("GelMaterial"    , thePTP->GelMaterial);
    cfg.lookupValue("GlassMaterial"    , thePTP->GlassMaterial);
    cfg.lookupValue("AbsMaterial"    , thePTP->AbsMaterial);

    cfg.lookupValue("AbsorberSurface"    , thePTP->AbsorberSurface);
    cfg.lookupValue("BentoSurface"       , thePTP->BentoSurface);



    if (!thePTP->IsOk)
      {
        clog << "KM3OpticalModuleMgr:WARNING in file: "<< ConfigFile <<  ": it does not contain all the parameters. The produced PM volume could be NULL, leading to an eventual segmentation fault." <<endl;
      }

    thePTP->PMTprt = (KM3PhotoTubeMgr::PhotoTubeParameters*)(KM3PhotoTubeMgr::GetIt ()->GetVolume (thePTP->PhotoTubeName)->GetParameters ());
    thePTP->CalculateParameters ();

  }

  KM3OpticalModuleMgr::~KM3OpticalModuleMgr ()
  {
  }

  KM3OpticalModuleMgr::KM3OpticalModule::KM3OpticalModule (OpticalModuleParameters *PTP): KM3ParametrizedVolume::KM3Volume(PTP)
  {
    cout << "creating it " << PTP->Prefix << endl;
  }

  KM3OpticalModuleMgr::KM3OpticalModule::~KM3OpticalModule ()
  {
  }


  G4LogicalVolume* KM3OpticalModuleMgr::KM3OpticalModule::Construct(const VolumeParameters* PTP)
  {
    if (fTheLogicalVolume)
      return fTheLogicalVolume;

    if (PTP)
      fTheParameters = new OpticalModuleParameters (*((OpticalModuleParameters*)PTP));

    const OpticalModuleParameters &ThePTP=*((OpticalModuleParameters*)fTheParameters);

    if (ThePTP.DetectionUnitType == "DOM")
      return fTheLogicalVolume = GetDOM ();
    else if (ThePTP.DetectionUnitType == "OM")
      return fTheLogicalVolume = GetOM ();

    cerr << "KM3OpticalModuleMgr:ERROR, wrong specification for the DetectionUnitType: \"" << ThePTP.DetectionUnitType <<"\" should be DOM or OM" << endl;
    return NULL;
  }

  G4LogicalVolume* KM3OpticalModuleMgr::KM3OpticalModule::GetDOM ()
  {
    KM3Material* materialsManager = KM3Material::GetIt ();

    G4String XPPrefix="Antares";

    const OpticalModuleParameters &ThePTP=*((OpticalModuleParameters*)fTheParameters);

    double PMRPosition=ThePTP.SphereRmin - ThePTP.PMTprt->RadiusSphere - ThePTP.DistSpherePMT;

    G4LogicalVolume *logicPMT =
      KM3PhotoTubeMgr::GetIt ()->GetLogicalVolume (ThePTP.PhotoTubeName);

    G4LogicalVolume* logicPlasticCore;

    //return logicPMT;

    //return logicPMT;
    G4VSolid* solidOM = new G4Sphere(ThePTP.Prefix + fBentoName, 0, ThePTP.SphereRmax, 0, twopi, 0,
                                     pi);
    G4VSolid* solidContainer = new G4UnionSolid (ThePTP.Prefix + "OMContainer",
                                                 solidOM,
                                                 new G4Sphere(ThePTP.Prefix + fBentoName, 0, ThePTP.SphereRmax+1*um, 0, twopi,
                                                              pi/2-0.047, //containing the scotch add it in the conf files
                                                              2*0.047)
                                                 );



    G4VSolid* solidGelCore = new G4Sphere(ThePTP.Prefix + "GelCore1", 0, ThePTP.SphereRmin, 0, twopi, 0,
                                          pi);

    G4VSolid* solidOMAirGap = new G4Tubs (ThePTP.Prefix + "OMAirGap",
                                          ThePTP.SphereRmin, ThePTP.SphereRmax,
                                          0.5*um,
                                          0, 2*pi);

    G4VSolid* solidPlasticCore   = new G4Sphere(ThePTP.Prefix + "PlasticCore1", PMRPosition*1.1, ThePTP.SphereRmin - ThePTP.DistAbsorberGlass,
                                                0, twopi, 0, pi);

    G4VSolid* solidMushroom = new G4Sphere(ThePTP.Prefix + "Mushroom",
                                           PMRPosition*1.1, ThePTP.SphereRmin - ThePTP.DistAbsorberGlass,
                                           0, twopi, 0,ThePTP.MushroomAngularSize);



    G4VSolid* SolidTubeGel = new G4Tubs(ThePTP.Prefix + "TubeGel",
                                        0., ThePTP.ReflectorRadiusMin,         //from article; For the reflector
                                        2*ThePTP.SphereRmin - PMRPosition,
                                        0, 2 * pi);

    //ConRadiusMax --  bigger radius of the reflector ring
    double ConRadiusMax =                      
      ThePTP.ReflectorRadiusMin
      + ThePTP.ReflectorThickness*tan(ThePTP.ReflectorAngle);
    G4VSolid* SolidBigTubeGel =
      new G4Tubs(ThePTP.Prefix + "TubeGel",
                 0., ConRadiusMax,
                 ThePTP.SphereRmin,
                 0, 2 * pi);


    G4Cons* SolidConeSub =
      new G4Cons(ThePTP.Prefix + "reflectorSub",
                 0, ThePTP.ReflectorRadiusMin,
                 0, ConRadiusMax,
                 ThePTP.ReflectorThickness/2,
                 0., 2 * pi);

    //to remove
    G4Box* solidBoahCutter = new G4Box(ThePTP.Prefix + "GelCutter",
                                       1000,
                                       1,
                                       1);

    G4VSolid* solidConRefl =
      new G4Cons(ThePTP.Prefix + "reflector",
                 ThePTP.ReflectorRadiusMin-0.001,ThePTP.ReflectorRadiusMin,
                 ConRadiusMax-0.001, ConRadiusMax,
                 ThePTP.ReflectorThickness/2,
                 0., 2 * pi);

    solidConRefl = new G4SubtractionSolid("reflector",solidConRefl,solidBoahCutter,0,G4ThreeVector());

    G4RotationMatrix PMOrientation (0,0,0);
    G4ThreeVector PMPosition(0,0,PMRPosition);
    G4ThreeVector ConPosition, BigTubePosition;

    //materialsManager->GetMaterial (ThePTP.GlassMaterial)->GetMaterialPropertiesTable ()->DumpTable (); exit (0);

    G4LogicalVolume* logicOMContainer=new G4LogicalVolume (solidContainer, materialsManager->GetMaterial (ThePTP.AbsMaterial), ThePTP.Prefix + "OMContainer");
    G4LogicalVolume* logicOM=new G4LogicalVolume (solidOM, materialsManager->GetMaterial (ThePTP.GlassMaterial), ThePTP.Prefix + fBentoName);
    G4LogicalVolume* logicGelCore=new G4LogicalVolume (solidGelCore, materialsManager->GetMaterial (ThePTP.GelMaterial), ThePTP.Prefix + "Gel");
    G4LogicalVolume* logicConRefl=new G4LogicalVolume (solidConRefl, materialsManager->GetMaterial (ThePTP.AbsMaterial), ThePTP.Prefix + "Reflector");


    int number=0;
    for (auto PMTsAngles : ThePTP.PMTsAngles)
      {
        PMPosition.setRThetaPhi (PMRPosition,
                                 PMTsAngles.first,
                                 PMTsAngles.second);

        //The distance from Reflector Ring to the glass is given as a minimal distance.
        //We need to recalculate it as a horisontal shift.
        double ReflectorDistGlassHor = sqrt(ThePTP.SphereRmin*ThePTP.SphereRmin-ConRadiusMax*ConRadiusMax)
            - sqrt((ThePTP.SphereRmin-ThePTP.ReflectorDistGlass)*(ThePTP.SphereRmin-ThePTP.ReflectorDistGlass)-ConRadiusMax*ConRadiusMax);
        cout << "ReflectorDistGlassHor is " << ReflectorDistGlassHor << endl;
        cout << "DistAbsorberGlass is " << ThePTP.DistAbsorberGlass << endl;
                    

        ConPosition.setRThetaPhi (ThePTP.SphereRmin*cos(asin(ConRadiusMax/ThePTP.SphereRmin))
                                  -ThePTP.ReflectorThickness/2-ReflectorDistGlassHor,
                                  PMTsAngles.first,
                                  PMTsAngles.second);

        BigTubePosition.setRThetaPhi (ThePTP.SphereRmin*cos(asin(ConRadiusMax/ThePTP.SphereRmin))
                                      +ThePTP.SphereRmin-ReflectorDistGlassHor,
                                      PMTsAngles.first,
                                      PMTsAngles.second);


        PMOrientation.set (PMTsAngles.second+M_PI/2,
                           PMTsAngles.first,0);

        SolidTubeGel = new G4Tubs(ThePTP.Prefix + "TubeGel",
                                  0., ThePTP.ReflectorRadiusMin,         //from article; For the reflector
                                  2*ThePTP.SphereRmin - PMRPosition,
                                  0, 2 * pi);
        number++;
        ostringstream oss;
        oss<<number;
        G4String Snumber; Snumber=Snumber+oss.str ();

        solidPlasticCore = new G4SubtractionSolid (ThePTP.Prefix + solidPlasticCore->GetName()+Snumber, solidPlasticCore, SolidTubeGel,
                                                   &PMOrientation, PMPosition);

        solidPlasticCore = new G4SubtractionSolid (ThePTP.Prefix + solidPlasticCore->GetName(), solidPlasticCore, SolidConeSub,
                                                   &PMOrientation, ConPosition);
        solidPlasticCore = new G4SubtractionSolid (ThePTP.Prefix + solidPlasticCore->GetName(), solidPlasticCore, SolidBigTubeGel,
                                                   &PMOrientation, BigTubePosition);

        if (ThePTP.AreThereReflectors)
          {
            new G4PVPlacement(new G4RotationMatrix(PMOrientation), ConPosition, logicConRefl, ThePTP.Prefix + "Reflector" + Snumber, logicGelCore, true, number, ThePTP.CheckOverlaps);
          }

        new G4PVPlacement(new G4RotationMatrix(PMOrientation), PMPosition, logicPMT, ThePTP.Prefix + fPMTContainerName + Snumber, logicGelCore, true, number, ThePTP.CheckOverlaps);
        cout << "pm number " << Snumber << ' ' << number << " at orientation " << PMOrientation.theta ()*180/M_PI
             << ' ' << PMOrientation.phi () * 180/M_PI<< endl
             << "at position " << PMPosition.x () << ';'<<PMPosition.y () << ';' << PMPosition.z () << endl;
      }





    logicPlasticCore=new G4LogicalVolume (solidPlasticCore, materialsManager->GetMaterial (ThePTP.AbsMaterial), ThePTP.Prefix + fAbsorberName);

    G4LogicalVolume* logicMushroom = new G4LogicalVolume(solidMushroom, materialsManager->GetMaterial (ThePTP.MushroomMaterial), ThePTP.Prefix + "Mushroom");
    G4LogicalVolume* logicOMAirGap = new G4LogicalVolume (solidOMAirGap, materialsManager->GetMaterial ("Vacuum"), ThePTP.Prefix + "OMAirGap");

    new G4PVPlacement(0, G4ThreeVector(), logicPlasticCore, ThePTP.Prefix + "PlasticCore", logicGelCore, false, 0, ThePTP.CheckOverlaps);
    new G4PVPlacement(0, G4ThreeVector(), logicOMAirGap, ThePTP.Prefix + "OMAirGap", logicOM, false, 0, ThePTP.CheckOverlaps);
    new G4PVPlacement(0, G4ThreeVector(), logicGelCore, ThePTP.Prefix + "GelCore", logicOM, false, 0, ThePTP.CheckOverlaps);
    new G4PVPlacement(0, G4ThreeVector(), logicOM, ThePTP.Prefix + fBentoName, logicOMContainer, false, 0, ThePTP.CheckOverlaps);
    new G4PVPlacement(0, G4ThreeVector(), logicMushroom, ThePTP.Prefix + "Mushroom", logicPlasticCore, false, 0, ThePTP.CheckOverlaps);

    logicGelCore->SetVisAttributes (fsilver);
    G4VisAttributes *att=new G4VisAttributes(logicPlasticCore->GetVisAttributes ());
    att->SetForceLineSegmentsPerCircle (1);
    logicPlasticCore->SetVisAttributes (att);


    G4OpticalSurface* ConeMirror_opsurf = materialsManager->GetOpticalSurface ("MirrorAluminium");
    G4OpticalSurface* MushroomSurf = materialsManager->GetOpticalSurface ("RoughAluminium");
    G4OpticalSurface* AirGapSurf = materialsManager->GetOpticalSurface ("RoughAirGap");

    new G4LogicalSkinSurface(ThePTP.Prefix + "AirGapSurf"  , logicOMAirGap, AirGapSurf);
    new G4LogicalSkinSurface(ThePTP.Prefix + "MushroomSurf", logicMushroom, MushroomSurf);
    new G4LogicalSkinSurface(ThePTP.Prefix + "ReflectSurf" , logicConRefl, ConeMirror_opsurf);
    if (ThePTP.AbsorberSurface.size ())
      {
        new G4LogicalSkinSurface("AbsorberSurface", logicPlasticCore, KM3Material::GetIt()->GetOpticalSurface (ThePTP.AbsorberSurface));
      }

    if (ThePTP.BentoSurface.size ())
      {
        new G4LogicalSkinSurface("BentoSurface", logicOM, KM3Material::GetIt()->GetOpticalSurface (ThePTP.BentoSurface));
      }

    return logicOMContainer;
  }


  G4LogicalVolume* KM3OpticalModuleMgr::KM3OpticalModule::GetOM () {

    KM3Material* materialsManager = KM3Material::GetIt ();

    const OpticalModuleParameters &ThePTP=*((OpticalModuleParameters*)fTheParameters);
    //PMtube first
    G4LogicalVolume *logicPMT =
      KM3PhotoTubeMgr::GetIt ()->GetLogicalVolume (ThePTP.PhotoTubeName);

    //solids makings for the outer module, absober and optical gel
    G4Sphere* solidOM = new G4Sphere(ThePTP.Prefix + fBentoName, 0, ThePTP.SphereRmax, 0, twopi, 0,
                                     pi);
    G4Sphere* solidOMAir = new G4Sphere(ThePTP.Prefix + "OM", 0, ThePTP.SphereRmin, 0, twopi, 0,
                                        pi);
    /*G4Sphere* solidBento = new G4Sphere(ThePTP.Prefix + fBentoName, fSphereRmin, fSphereRmax, 0,
      twopi, 0, pi);*/
    G4Sphere* solidAbsorber = new G4Sphere(ThePTP.Prefix + fAbsorberName+"full",
                                           ThePTP.SphereRmin - 1*mm, ThePTP.SphereRmin,
                                           0., twopi, pi/2 - ThePTP.AbsorberOvertakeAngle  , pi);
    G4Sphere* solidGel = new G4Sphere(ThePTP.Prefix + "Gelextern",
                                      0., ThePTP.SphereRmin,
                                      0, twopi, 0, pi/2);
    G4Box* solidGelCutter = new G4Box(ThePTP.Prefix + "GelCutter",
                                      ThePTP.SphereRmax,
                                      ThePTP.SphereRmax,
                                      2*ThePTP.SphereRmax-ThePTP.SphereRmin-ThePTP.DistSpherePMT-ThePTP.PMTprt->PhotoHeight-ThePTP.PMTprt->BulbHeight/2-ThePTP.GelOvertake);


    /*//substraction to remove overlaps
      G4SubtractionSolid* Absorbersubtraction =
      new G4SubtractionSolid(ThePTP.Prefix + fAbsorberName,
      solidAbsorber, logicPMT->GetSolid(),
      0, G4ThreeVector(0,0,ThePTP.SphereRmin-ThePTP.PMTprt->RadiusSphere-ThePTP.DistSpherePMT));
    */
    G4SubtractionSolid* Gelsubtraction =
      new G4SubtractionSolid(ThePTP.Prefix + "GelSub1",
                             solidGel, solidGelCutter,
                             0, G4ThreeVector(0,0,0));
    Gelsubtraction =
      new G4SubtractionSolid(ThePTP.Prefix + "GelSub",
                             Gelsubtraction, logicPMT->GetSolid(),
                             0, G4ThreeVector(0,0,ThePTP.SphereRmin-ThePTP.PMTprt->RadiusSphere-ThePTP.DistSpherePMT));



    //logical and placements of everything
    string fXPPrefix = "Antares";
    G4LogicalVolume* logicOM =
      new G4LogicalVolume(solidOM,
                          materialsManager->GetMaterial(ThePTP.GlassMaterial), ThePTP.Prefix + fBentoName);
    G4LogicalVolume* logicOMAir =
      new G4LogicalVolume(solidOMAir,
                          materialsManager->GetMaterial("Air"), ThePTP.Prefix + "OM");

    /*    G4LogicalVolume* logicBento =
          new G4LogicalVolume(solidBento,
          materialsManager->GetMaterial(fXPPrefix + "Glass"), prefix + fBentoName);
    */
    G4LogicalVolume* logicAbsorber =
      new G4LogicalVolume(solidAbsorber,
                          materialsManager->GetMaterial(ThePTP.AbsMaterial), ThePTP.Prefix + fAbsorberName);


    G4LogicalVolume* logicGel =
      new G4LogicalVolume(Gelsubtraction,
                          materialsManager->GetMaterial(ThePTP.GelMaterial), ThePTP.Prefix + "Gel");

    G4VPhysicalVolume *physiOMAir =
      new G4PVPlacement(0, G4ThreeVector(),
                        logicOMAir, ThePTP.Prefix + "OMAir", logicOM, false, 0, ThePTP.CheckOverlaps);

    new G4PVPlacement(0, G4ThreeVector(),
                      logicAbsorber, ThePTP.Prefix + fAbsorberName, logicOMAir, false, 0, ThePTP.CheckOverlaps);


    G4VPhysicalVolume *physiGel=
      new G4PVPlacement(0, G4ThreeVector(0,0,0),
                        logicGel, ThePTP.Prefix + "Gel", logicOMAir, false, 0, ThePTP.CheckOverlaps);




    new G4PVPlacement(0,
                      G4ThreeVector(0., 0., ThePTP.SphereRmin-ThePTP.PMTprt->RadiusSphere-ThePTP.DistSpherePMT),
                      logicPMT,
                      ThePTP.Prefix + fPMTContainerName,
                      logicOMAir,
                      false, 0, ThePTP.CheckOverlaps);

    //create the absorber surface
    if (ThePTP.AbsorberSurface.size ())
      {
        new G4LogicalSkinSurface("AbsorberSurface", logicAbsorber, KM3Material::GetIt()->GetOpticalSurface (ThePTP.AbsorberSurface));
      }

    //create the gel/air interface
    G4OpticalSurface* GelSurf = new G4OpticalSurface(ThePTP.Prefix + "GelSurf");
    G4MaterialPropertiesTable* PropGel = new G4MaterialPropertiesTable();
    GelSurf->SetType(dielectric_dielectric);
    GelSurf->SetModel (glisur);
    GelSurf->SetFinish(ground);	// needed for mirror
    GelSurf->SetPolish(.8);

    GelSurf->SetMaterialPropertiesTable(PropGel);

    new G4LogicalBorderSurface(ThePTP.Prefix + "GelAirSurf", physiOMAir, physiGel, GelSurf);
    new G4LogicalBorderSurface(ThePTP.Prefix + "AirGelSurf", physiGel, physiOMAir, GelSurf);



    return logicOM;

  }


}
