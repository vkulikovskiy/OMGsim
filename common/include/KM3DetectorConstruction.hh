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
#ifndef KM3DetectorConstruction_h
#define KM3DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "G4UnionSolid.hh"

#include <string>
#include <map>

#include <KM3ParametrizedVolume.hh>

//#include "KM3MagneticField.hh"

namespace libconfig
{
  class Config;
}

class G4Tubs;
class G4Box;
class G4Sphere;
class G4Cons;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4Region;
class G4UserLimits;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace km3net
{
  class OpticalModuleMaker;
  class KM3DetectorMessenger;
  class KM3Material;

  class KM3DetectorConstruction: public G4VUserDetectorConstruction, public KM3ParametrizedVolume {
  public:
    static KM3DetectorConstruction* GetIt ();
  private:
    KM3DetectorConstruction();
    ~KM3DetectorConstruction();

    static KM3DetectorConstruction *fInstance;

    //exeption in the KM3ParametrizedVolume familly, no ConstructTheVolume and no construct.
    virtual G4LogicalVolume* Construct (const VolumeParameters* = NULL) {return 0;}
    virtual KM3Volume* ConstructTheVolume (VolumeParameters*) {return 0;}

  public:
    virtual G4VPhysicalVolume* Construct();
    G4LogicalVolume* ConstructOpticalModule1();

    const G4VPhysicalVolume* GetDetector() {
      return physiDetector;
    }
    ;

    G4double GetTargetFullLength() const {
      return fTargetLength;
    }
    ;
    G4double GetWorldFullLength() {
      return fWorldLength;
    }
    ;
    G4double GetSphereRMax() {
      return fSphereRmax;
    }
    ;

    void SetMeanCapRadius (G4double R);
    void DontDraw (G4String DD);
    void SetTargetMaterial(G4String);

    void SetTargetLength(G4double value);
    void SetWorldLength(G4double value) {
      fWorldLength = value;
    }
    ;

    void SetOrientationTheta(G4double value) {
      FOrientationTheta = value;
    }
    void SetOrientationPsi(G4double value) {
      FOrientationPsi = value;
    }
    void SetOrientationPhi(G4double value) {
      FOrientationPhi = value;
    }
    ;

    void SetGenerationVolumeLength(G4double value) {
      fTargetLength=value;
    }
    ;
    void SetGenerationVolumeRadius(G4double value) {
      SetGenerationVolumeLength(value);
    }
    ;


    void SetDetectorTypeName (std::string NAME)
    {
      fDetectorTypeName=NAME;
    }

    static const std::string fOMBasePrefix;

    void SetDataFilesSources (std::string Path);

  private:
    void DefineMaterials();
    void ProcessLogicalVolume(G4LogicalVolume* volume);

    struct PositionParameters: public KM3ParametrizedVolume::VolumeParameters
    {
      bool isOk=true;
      int ID=0;
      double X=0;
      double Y=0;
      double Z=0;
      double Theta=0;
      double Phi=0;
      double Psi=0;
      std::string Name="AntaresOM";
      std::string MotherVolume = "";
      std::string BorderSurfaceMaterial;
      std::string BorderSurfaceWith;
    };

    std::string fDetectorTypeName = "K40cfg";
    std::map<std::string, std::vector<PositionParameters> > fOMList;
    std::map<std::string, std::vector<PositionParameters> > fPMList;
    std::map<std::string, std::vector<PositionParameters> > fVolumeList;

    void ReadConfiguration (std::string);
    void ReadThisConfiguration (libconfig::Setting* cfg,
                                std::vector<PositionParameters>&  ConfList);
    G4VPhysicalVolume* PlaceThem (KM3ParametrizedVolume*, PositionParameters&);

  protected:

    //     G4Tubs*             solidWorld;    // pointer to the solid envelope
    G4Box* solidWorld;    // pointer to the solid envelope
    G4LogicalVolume* logicWorld;    // pointer to the logical envelope
    G4VPhysicalVolume* physiWorld;    // pointer to the physical envelope

    //     G4Tubs*             solidTarget;   // pointer to the solid Target
  public: //!!!!!!!!!!!!!! attention, temporary stuff
    G4VSolid* solidTarget;   // pointer to the solid Target
    G4LogicalVolume* logicTarget;   // pointer to the logical Target
    G4VPhysicalVolume* physiTarget;   // pointer to the physical Target
  protected:

    G4Sphere* solidDetector;  // pointer to the solid Detector
    G4LogicalVolume* logicDetector;  // pointer to the logical Detector
    G4VPhysicalVolume* physiDetector;  // pointer to the physical Detector

    KM3DetectorMessenger* detectorMessenger;  // pointer to the Messenger
    KM3Material* fMaterialsManager;         // material manager

    G4Material* DefaultMater;          // Default material
    G4Material* TargetMater = 0;          // Target material
    G4Material* DetectorMater;         // Detector material

    G4double fWorldLength;            // Full length the world volume
    static G4double fTargetLength;           // Full length of the target
    G4double fWorldRadius;            // Radius of the world volume
    G4double fTargetRadius;           // Radius of the target
    G4double fDetectorThickness;      // Thickness of the Detector
    G4ThreeVector globpos; // vector that shifts storey with PMT to center of the target

    G4Region* targetRegion;
    G4Region* detectorRegion;
    G4UserLimits* stepLimit;             // pointer to user step limits

    G4ThreeVector fOM1Placement;

    G4double fMeanCapRadius=0;

    G4double fSphereRmax;
    G4double FStoreyVolume;

    G4double FOrientationTheta;
    G4double FOrientationPhi=0;
    G4double FOrientationPsi=0;
    G4double fOMCubicVolume;

    G4double fGenerationVolumeLength; // Full length of the generation volume box
    G4double fGenerationVolumeRadius; // Radius of the generation volume sphere

    G4String fDD;

    OpticalModuleMaker* fOMM=NULL;

    std::map <std::string, G4VPhysicalVolume*>  fPlacedVolumeDict;
    std::map <std::string, std::string>         fBorderSurfaceList;

  };
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
