#include <KM3PassiveVolumeMgr.hh>
#include <KM3PhotoTubeMgr.hh>
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
#include <G4Ellipsoid.hh>
#include <G4Torus.hh>
#include <G4Polyhedra.hh>
#include <G4SubtractionSolid.hh>
#include <G4IntersectionSolid.hh>
#include <G4UnionSolid.hh>

#include <G4LogicalVolume.hh>

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
  KM3PassiveVolumeMgr::KM3ParametrizedVolume* KM3PassiveVolumeMgr::fInstance = NULL;

  string KM3PassiveVolumeMgr::fname = "";

  KM3PassiveVolumeMgr* KM3PassiveVolumeMgr::GetIt ()
  {
    if (!fInstance)
      fInstance=new KM3PassiveVolumeMgr ();
    return (KM3PassiveVolumeMgr*)fInstance;
  }

  KM3PassiveVolumeMgr::KM3PassiveVolumeMgr ()
  {
    ConfigFilePrefix="KM3Vol";

    LoadDirectory ((string(PROJECT_SOURCE_DIR)+"/common/data/"));
  }

  KM3PassiveVolumeMgr::KM3Volume* KM3PassiveVolumeMgr::ConstructTheVolume (VolumeParameters* PTP)
  {
    return new KM3PassiveVolume ((PassiveVolumeParameters*)PTP);
  }

  void KM3PassiveVolumeMgr::ReadConfiguration (string ConfigFile)
  {
    Config cfg;
    cfg.setAutoConvert (true);
    cfg.readFile (ConfigFile.data ());

    Setting& rootStg=cfg.getRoot();

    if (rootStg.lookupValue("Name", fname))
      {
        fParametersDict[fname]=new PassiveVolumeParameters;
        ReadCurrentConfiguration (&rootStg);
        if (!fParametersDict[fname]->IsOk)
          fParametersDict.erase(fname);
        return;
      }

    int configFileSize=rootStg.getLength ();
    bool AtLeastOne = false;
    for (int cfIt=0; cfIt < configFileSize; cfIt++)
      {
        if (!rootStg[cfIt].isGroup ())
          continue;

        if (rootStg[cfIt].lookupValue("Name", fname))
          {
            fParametersDict[fname]=new PassiveVolumeParameters;
            ReadCurrentConfiguration (&(rootStg[cfIt]));
            if (!fParametersDict[fname]->IsOk)
              fParametersDict.erase(fname);
            else
              AtLeastOne = true;
          }

      }
    if (AtLeastOne)
      return;

    cerr << "KM3PassiveVolumeMgr  no Name parameter for " << ConfigFile<< endl;
    return;

  }

  void KM3PassiveVolumeMgr::ReadCurrentConfiguration (Setting* SET)
  {
    Setting& stg=*SET;
    PassiveVolumeParameters *thePTP=(PassiveVolumeParameters*)(fParametersDict[fname]);

    if (!ReadBaseParameters (&stg, thePTP))
      {
        cerr << "ReadBaseParameters: error: miss base parameters for " << fname << endl;
        thePTP->IsOk=false;
        return;
      }


    thePTP->Name=fname;
    //Checkout The kind of volume
    thePTP->IsOk =
      stg.lookupValue("VolumeKind"       , thePTP->VolumeKind       );

    if (! thePTP->IsOk)
      {
        cerr << "KM3PassiveVolumeMgr: config volume kind error \""<<thePTP->VolumeKind<<"\", please enter a valid string." << endl;
        return;
      }

    //parse to the volume generation
    if (thePTP->VolumeKind == "Ellipsoid") ReadEllipsoidConfiguration (&stg);
    else if (thePTP->VolumeKind == "Box") ReadBoxConfiguration (&stg);
    else if (thePTP->VolumeKind == "Con") ReadConConfiguration (&stg);
    else if (thePTP->VolumeKind == "Tub") ReadTubConfiguration (&stg);
    else if (thePTP->VolumeKind == "Sphere") ReadSphereConfiguration (&stg);
    else if (thePTP->VolumeKind == "Torus") ReadTorusConfiguration (&stg);
    else if (thePTP->VolumeKind == "Hexagone") ReadHexagoneConfiguration (&stg);
    else if (thePTP->VolumeKind == "Union"       ||
             thePTP->VolumeKind == "Subtraction" ||
             thePTP->VolumeKind == "Intersection"  ) ReadBooleanConfiguration (&stg);

    if (! thePTP->IsOk)
      cerr << "KM3PassiveVolumeMgr: parsing error for the volume " << thePTP->VolumeKind << endl;
    stg.lookupValue("Material"    , thePTP->Material);
    stg.lookupValue("SkinMaterial"    , thePTP->SkinMaterial);

    if (!thePTP->IsOk)
      {
        clog << "KM3PassiveVolumeMgr:WARNING in configuration: "<< stg.getName() <<  ": it does not contain all the parameters. The produced PM volume could be NULL, leading to an eventual segmentation fault." <<endl;
      }

  }

  void KM3PassiveVolumeMgr::ReadEllipsoidConfiguration(Setting* stg)
  {
    PassiveVolumeParameters *thePTP=(PassiveVolumeParameters*)(fParametersDict[fname]);
    thePTP->IsOk =
      stg->exists("RadiusMax") &&
      stg->exists("SecondRadiusMax") &&
      stg->exists("ThirdRadiusMax") &&
      stg->exists("ZCutMin") &&
      stg->exists("ZCutMax") ;

    if (!thePTP->IsOk)
      return;

    thePTP->RadiusMax = GetFloatingParameter
      (stg->lookup ("RadiusMax"));
    thePTP->SecondRadiusMax = GetFloatingParameter
      (stg->lookup ("SecondRadiusMax"));
    thePTP->ThirdRadiusMax = GetFloatingParameter
      (stg->lookup ("ThirdRadiusMax"));
    thePTP->ZCutMin = GetFloatingParameter
      (stg->lookup ("ZCutMin"));
    thePTP->ZCutMax = GetFloatingParameter
      (stg->lookup ("ZCutMax"));
  }
  void KM3PassiveVolumeMgr::ReadBoxConfiguration(Setting* stg)
  {
    PassiveVolumeParameters *thePTP=(PassiveVolumeParameters*)(fParametersDict[fname]);
    thePTP->IsOk =
      stg->exists("X") &&
      stg->exists("Y") &&
      stg->exists("Z") ;

    if (!thePTP->IsOk)
      return;

    thePTP->X = GetFloatingParameter
      (stg->lookup ("X"));
    thePTP->Y = GetFloatingParameter
      (stg->lookup ("Y"));
    thePTP->Z = GetFloatingParameter
      (stg->lookup ("Z"));

  }
  void KM3PassiveVolumeMgr::ReadConConfiguration(Setting* stg)
  {
    PassiveVolumeParameters *thePTP=(PassiveVolumeParameters*)(fParametersDict[fname]);
    thePTP->IsOk =
      stg->exists("RadiusMin") &&
      stg->exists("RadiusMax") &&
      stg->exists("SecondRadiusMin") &&
      stg->exists("SecondRadiusMax") &&
      stg->exists("PhiMin") &&
      stg->exists("PhiMax") &&
      stg->exists("Z") ;

    if (!thePTP->IsOk)
      return;

    thePTP->RadiusMin = GetFloatingParameter
      (stg->lookup ("RadiusMin"));
    thePTP->RadiusMax = GetFloatingParameter
      (stg->lookup ("RadiusMax"));
    thePTP->SecondRadiusMin = GetFloatingParameter
      (stg->lookup ("SecondRadiusMin"));
    thePTP->SecondRadiusMax = GetFloatingParameter
      (stg->lookup ("SecondRadiusMax"));
    thePTP->PhiMin = GetFloatingParameter
      (stg->lookup ("PhiMin"));
    thePTP->PhiMax = GetFloatingParameter
      (stg->lookup ("PhiMax"));
    thePTP->Z = GetFloatingParameter
      (stg->lookup ("Z"));

  }
  void KM3PassiveVolumeMgr::ReadTubConfiguration(Setting* stg)
  {
    PassiveVolumeParameters *thePTP=(PassiveVolumeParameters*)(fParametersDict[fname]);
    thePTP->IsOk =
      stg->exists("RadiusMin") &&
      stg->exists("RadiusMax") &&
      stg->exists("PhiMin") &&
      stg->exists("PhiMax") &&
      stg->exists("Z") ;

    if (!thePTP->IsOk)
      return;

    thePTP->RadiusMin = GetFloatingParameter
      (stg->lookup ("RadiusMin"));
    thePTP->RadiusMax = GetFloatingParameter
      (stg->lookup ("RadiusMax"));
    thePTP->PhiMin = GetFloatingParameter
      (stg->lookup ("PhiMin"));
    thePTP->PhiMax = GetFloatingParameter
      (stg->lookup ("PhiMax"));
    thePTP->Z = GetFloatingParameter
      (stg->lookup ("Z"));

  }
  void KM3PassiveVolumeMgr::ReadSphereConfiguration(Setting* stg)
  {
    PassiveVolumeParameters *thePTP=(PassiveVolumeParameters*)(fParametersDict[fname]);
    thePTP->IsOk =
      stg->exists("RadiusMin") &&
      stg->exists("RadiusMax") &&
      stg->exists("PhiMin") &&
      stg->exists("PhiMax") &&
      stg->exists("ThetaMin") &&
      stg->exists("ThetaMax") ;

    if (!thePTP->IsOk)
      return;

    thePTP->RadiusMin = GetFloatingParameter
      (stg->lookup ("RadiusMin"));
    thePTP->RadiusMax = GetFloatingParameter
      (stg->lookup ("RadiusMax"));
    thePTP->PhiMin = GetFloatingParameter
      (stg->lookup ("PhiMin"));
    thePTP->PhiMax = GetFloatingParameter
      (stg->lookup ("PhiMax"));
    thePTP->ThetaMin = GetFloatingParameter
      (stg->lookup ("ThetaMin"));
    thePTP->ThetaMax = GetFloatingParameter
      (stg->lookup ("ThetaMax"));
  }
  void KM3PassiveVolumeMgr::ReadTorusConfiguration(Setting* stg)
  {
    PassiveVolumeParameters *thePTP=(PassiveVolumeParameters*)(fParametersDict[fname]);
    thePTP->IsOk =
      stg->exists("RadiusMin") &&
      stg->exists("RadiusMax") &&
      stg->exists("SecondRadiusMin") &&
      stg->exists("PhiMin") &&
      stg->exists("PhiMax") ;

    if (!thePTP->IsOk)
      return;

    thePTP->RadiusMin = GetFloatingParameter
      (stg->lookup ("RadiusMin"));
    thePTP->RadiusMax = GetFloatingParameter
      (stg->lookup ("RadiusMax"));
    thePTP->SecondRadiusMin = GetFloatingParameter
      (stg->lookup ("SecondRadiusMin"));
    thePTP->PhiMin = GetFloatingParameter
      (stg->lookup ("PhiMin"));
    thePTP->PhiMax = GetFloatingParameter
      (stg->lookup ("PhiMax"));

  }

  void KM3PassiveVolumeMgr::ReadHexagoneConfiguration(Setting* stg)
  {
    PassiveVolumeParameters *thePTP=(PassiveVolumeParameters*)(fParametersDict[fname]);
    thePTP->IsOk =
      stg->exists("RadiusMax") &&
      stg->exists("Z");

    if (!thePTP->IsOk)
      return;
    thePTP->RadiusMin = 0;
    thePTP->RadiusMax = GetFloatingParameter
      (stg->lookup ("RadiusMax"));
    thePTP->RadiusMin = GetFloatingParameter
      (stg->lookup ("RadiusMin"));

    thePTP->SecondRadiusMin=thePTP->RadiusMin;
    thePTP->SecondRadiusMax=thePTP->RadiusMax;
    if (stg->exists ("SecondRadiusMax"))
      thePTP->SecondRadiusMax = GetFloatingParameter
        (stg->lookup ("SecondRadiusMax"));
    if (stg->exists ("SecondRadiusMin"))
      thePTP->SecondRadiusMin = GetFloatingParameter
        (stg->lookup ("SecondRadiusMin"));


    thePTP->Z = GetFloatingParameter
      (stg->lookup ("Z"));
    if (stg->exists ("PhiMin"))
      thePTP->PhiMin = GetFloatingParameter
        (stg->lookup ("PhiMin"));
  }

  void KM3PassiveVolumeMgr::ReadBooleanConfiguration(Setting* stg)
  {
    PassiveVolumeParameters *thePTP=(PassiveVolumeParameters*)(fParametersDict[fname]);
    thePTP->IsOk =
      stg->exists("X") &&
      stg->exists("Y") &&
      stg->exists("Z") &&
      stg->exists("RotTheta") &&
      stg->exists("RotPhi") &&
      stg->exists("RefVolume") &&
      stg->exists("OperandVolume") ;

    if (!thePTP->IsOk)
      return;

    thePTP->X = GetFloatingParameter
      (stg->lookup ("X"));
    thePTP->Y = GetFloatingParameter
      (stg->lookup ("Y"));
    thePTP->Z = GetFloatingParameter
      (stg->lookup ("Z"));
    thePTP->RotTheta = GetFloatingParameter
      (stg->lookup ("RotTheta"));
    thePTP->RotPhi = GetFloatingParameter
      (stg->lookup ("RotPhi"));
    stg->lookupValue ("RefVolume",thePTP->RefVolume);
    stg->lookupValue ("OperandVolume",thePTP->OperandVolume);
  }


  KM3PassiveVolumeMgr::~KM3PassiveVolumeMgr ()
  {
  }

  KM3PassiveVolumeMgr::KM3PassiveVolume::KM3PassiveVolume (PassiveVolumeParameters *PTP): KM3ParametrizedVolume::KM3Volume(PTP)
  {
    cout << "creating it " << PTP->Prefix << endl;
  }

  KM3PassiveVolumeMgr::KM3PassiveVolume::~KM3PassiveVolume ()
  {
  }


  G4LogicalVolume* KM3PassiveVolumeMgr::KM3PassiveVolume::Construct(const VolumeParameters* PTP)
  {
    if (fTheLogicalVolume)
      return fTheLogicalVolume;

    if (PTP)
      fTheParameters = new PassiveVolumeParameters (*((PassiveVolumeParameters*)PTP));

    const PassiveVolumeParameters &thePTP=*((PassiveVolumeParameters*)fTheParameters);

    if (thePTP.VolumeKind == "Ellipsoid")        return fTheLogicalVolume = ConstructEllipsoidVolume ();
    else if (thePTP.VolumeKind == "Box")         return fTheLogicalVolume = ConstructBoxVolume ();
    else if (thePTP.VolumeKind == "Con")         return fTheLogicalVolume = ConstructConVolume ();
    else if (thePTP.VolumeKind == "Tub")         return fTheLogicalVolume = ConstructTubVolume ();
    else if (thePTP.VolumeKind == "Sphere")      return fTheLogicalVolume = ConstructSphereVolume ();
    else if (thePTP.VolumeKind == "Torus")       return fTheLogicalVolume = ConstructTorusVolume ();
    else if (thePTP.VolumeKind == "Hexagone")    return fTheLogicalVolume = ConstructHexagoneVolume ();
    else if (thePTP.VolumeKind == "Union"       ||
             thePTP.VolumeKind == "Subtraction" ||
             thePTP.VolumeKind == "Intersection"  ) return fTheLogicalVolume = ConstructBooleanVolume ();


    cerr << "KM3PassiveVolumeMgr:ERROR, wrong specification for the operation: \"" << thePTP.VolumeKind <<"\" should be a geant4 solid volume name" << endl;
    return NULL;
  }

  G4LogicalVolume* KM3PassiveVolumeMgr::KM3PassiveVolume::ConstructEllipsoidVolume ()
  {
    const PassiveVolumeParameters &thePTP=*((PassiveVolumeParameters*)fTheParameters);
    return MakeItLogical (
                          new G4Ellipsoid
                          (thePTP.Name, thePTP.RadiusMax, thePTP.SecondRadiusMax, thePTP.ThirdRadiusMax, thePTP.ZCutMin,thePTP.ZCutMax)
                          );
  }

  G4LogicalVolume* KM3PassiveVolumeMgr::KM3PassiveVolume::ConstructBoxVolume ()
  {
    const PassiveVolumeParameters &thePTP=*((PassiveVolumeParameters*)fTheParameters);
    return MakeItLogical (
                          new G4Box
                          (thePTP.Name, thePTP.X, thePTP.Y, thePTP.Z)
                          );
  }

  G4LogicalVolume* KM3PassiveVolumeMgr::KM3PassiveVolume::ConstructConVolume ()
  {
    const PassiveVolumeParameters &thePTP=*((PassiveVolumeParameters*)fTheParameters);
    return MakeItLogical (
                          new G4Cons
                          (thePTP.Name, thePTP.RadiusMin, thePTP.RadiusMax, thePTP.SecondRadiusMin, thePTP.SecondRadiusMax, thePTP.SecondRadiusMax, thePTP.PhiMin, thePTP.PhiMax)
                          );
  }

  G4LogicalVolume* KM3PassiveVolumeMgr::KM3PassiveVolume::ConstructTubVolume ()
  {
    const PassiveVolumeParameters &thePTP=*((PassiveVolumeParameters*)fTheParameters);
    return MakeItLogical (
                          new G4Tubs
                          (thePTP.Name, thePTP.RadiusMin, thePTP.RadiusMax, thePTP.Z, thePTP.PhiMin, thePTP.PhiMax)
                          );
  }

  G4LogicalVolume* KM3PassiveVolumeMgr::KM3PassiveVolume::ConstructSphereVolume ()
  {
    const PassiveVolumeParameters &thePTP=*((PassiveVolumeParameters*)fTheParameters);
    return MakeItLogical (
                          new G4Sphere
                          (thePTP.Name, thePTP.RadiusMin, thePTP.RadiusMax, thePTP.PhiMin, thePTP.PhiMax, thePTP.ThetaMin, thePTP.ThetaMax)
                          );
  }

  G4LogicalVolume* KM3PassiveVolumeMgr::KM3PassiveVolume::ConstructTorusVolume ()
  {
    const PassiveVolumeParameters &thePTP=*((PassiveVolumeParameters*)fTheParameters);
    return MakeItLogical (
                          new G4Torus
                          (thePTP.Name, thePTP.RadiusMin, thePTP.RadiusMax, thePTP.SecondRadiusMax, thePTP.PhiMin, thePTP.PhiMax)
                          );
  }

  G4LogicalVolume* KM3PassiveVolumeMgr::KM3PassiveVolume::ConstructHexagoneVolume ()
  {
    const PassiveVolumeParameters &thePTP=*((PassiveVolumeParameters*)fTheParameters);
    double z[2] = {0,thePTP.Z};

    double r1=thePTP.RadiusMin;
    double r2=thePTP.RadiusMax;
    double r21=thePTP.SecondRadiusMin;
    double r22=thePTP.SecondRadiusMax;

    double rmax[2] = {r1,r21};
    double rmin[2] = {r2,r22};
    return MakeItLogical (
                          new G4Polyhedra
                          (thePTP.Name.data (), thePTP.PhiMin, 6.3, 6, 2, z, rmax, rmin)
                          );
  }

  G4LogicalVolume* KM3PassiveVolumeMgr::KM3PassiveVolume::ConstructBooleanVolume ()
  {
    //manage with a rotation axis!
    const PassiveVolumeParameters &thePTP=*((PassiveVolumeParameters*)fTheParameters);
    G4LogicalVolume* refLvol = KM3PassiveVolumeMgr::GetIt ()->GetLogicalVolume (thePTP.RefVolume);
    G4LogicalVolume* opLvol = KM3PassiveVolumeMgr::GetIt ()->GetLogicalVolume (thePTP.OperandVolume);
    if (refLvol == NULL)
      {
        clog << "KM3PassiveVolumeMgr::KM3PassiveVolume::ConstructBooleanVolume: warning: trying to get the ref volume among the PMTs and OMs " <<thePTP.RefVolume<< endl;
        refLvol = KM3PhotoTubeMgr::GetIt ()->GetLogicalVolume (thePTP.RefVolume);
            cout << "here the refvol is " << refLvol << endl;
        if (refLvol == NULL)
          {
            refLvol = KM3OpticalModuleMgr::GetIt ()->GetLogicalVolume (thePTP.RefVolume);
          }
      }
    if (opLvol == NULL)
      {
        clog << "KM3PassiveVolumeMgr::KM3PassiveVolume::ConstructBooleanVolume: warning: trying to get the operand volume among the PMTs and OMs " << thePTP.OperandVolume<< endl;
        opLvol = KM3PhotoTubeMgr::GetIt ()->GetLogicalVolume (thePTP.OperandVolume);
        if (opLvol == NULL)
          {
            opLvol = KM3OpticalModuleMgr::GetIt ()->GetLogicalVolume (thePTP.OperandVolume);
          }
      }


    G4VSolid* refVol = refLvol->GetSolid ();
    G4VSolid* opVol  = opLvol->GetSolid ();


    G4RotationMatrix* theRot = new G4RotationMatrix;
    theRot->setTheta(thePTP.RotTheta);
    theRot->setPhi(thePTP.RotPhi);

    return MakeItLogical (
                          [&] () -> G4VSolid*
                          {
                            if (thePTP.VolumeKind == "Union")
                              return new G4UnionSolid
                                (thePTP.Name,
                                 refVol, opVol,
                                 theRot,
                                 G4ThreeVector (thePTP.X,thePTP.Y,thePTP.Z));
                            else if (thePTP.VolumeKind == "Subtraction")
                              return new G4SubtractionSolid
                                (thePTP.Name,
                                 refVol, opVol,
                                 theRot,
                                 G4ThreeVector (thePTP.X,thePTP.Y,thePTP.Z));
                            else if (thePTP.VolumeKind == "Intersection")
                              return new G4IntersectionSolid
                                (thePTP.Name,
                                 refVol, opVol,
                                 theRot,
                                 G4ThreeVector (thePTP.X,thePTP.Y,thePTP.Z));
                            return 0;
                          } ()
                          );

  }

  G4LogicalVolume* KM3PassiveVolumeMgr::KM3PassiveVolume::MakeItLogical (G4VSolid* Solid)
  {
    const PassiveVolumeParameters &thePTP=*((PassiveVolumeParameters*)fTheParameters);
    KM3Material* materialsManager = KM3Material::GetIt ();
    G4LogicalVolume* toreturn = new G4LogicalVolume (Solid, materialsManager->GetMaterial (thePTP.Material), thePTP.Name);

    if (thePTP.SkinMaterial.size () == 0)
      return toreturn;

    cout << "doing the skin" << thePTP.SkinMaterial << endl;
    //    G4OpticalSurface* the_surface =new G4OpticalSurface(thePTP.Prefix + thePTP.SkinMaterial + "_surface",

    //the_surface->SetMaterialPropertiesTable(G4Material::GetMaterial(thePTP.SkinMaterial)->GetMaterialPropertiesTable());
    new G4LogicalSkinSurface(thePTP.Prefix + thePTP.SkinMaterial,
                             toreturn,
                             KM3Material::GetIt()->GetOpticalSurface (thePTP.SkinMaterial));
    return toreturn;

  }
}
