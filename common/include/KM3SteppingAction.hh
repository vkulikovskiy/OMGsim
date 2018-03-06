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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef KM3SteppingAction_h
#define KM3SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include <G4Types.hh>
#include <G4ThreeVector.hh>
#include <G4MaterialPropertyVector.hh>


#include <vector>

class G4MaterialPropertiesTable;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

namespace km3net
{
  class KM3SteppingAction : public G4UserSteppingAction
  {
  public:
    KM3SteppingAction();
    virtual ~KM3SteppingAction();

    virtual void UserSteppingAction(const G4Step*);

  protected:
    unsigned char CheckDetection(const G4ThreeVector& localPoint);
    G4MaterialPropertyVector* _angular_efficiency_err;
    G4MaterialPropertyVector* _angular_efficiency;
    G4MaterialPropertyVector* _efficiency_photocathode;
    G4MaterialPropertyVector* _thickness_photocathode;

    bool LoadMaterialProps (G4MaterialPropertiesTable*);

    double fEnergy=0;
    double fElectronAbsLength = 25e-6;
  };
}
#endif
