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
//
/// \copied from B5DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
#include "DriftChamberSD.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"

#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4PVReplica.hh"
#include "G4UserLimits.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4GenericMessenger.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction(), 
  dcin_wireplane_logical_(nullptr), dcout_wireplane_logical_(nullptr)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::~DetectorConstruction()
{
  //delete fMessenger;

  for (auto visAttributes: fVisAttributes) {
    delete visAttributes;
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Construct materials
  ConstructMaterials();
  auto air = G4Material::GetMaterial("G4_AIR");
  auto vacuum = G4Material::GetMaterial("G4_Galactic");
  auto argonGas = G4Material::GetMaterial("G4_Ar");
  auto scintillator = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  auto csI = G4Material::GetMaterial("G4_CESIUM_IODIDE");
  auto lead = G4Material::GetMaterial("G4_Pb");
  auto carbon = G4Material::GetMaterial("G4_C");
  auto hydrogen = new G4Material("hydrogne", 1., 1.01*g/mole, 1.*g/cm3);

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  // geometries --------------------------------------------------------------

  // experimental hall (world volume)
  auto worldSolid 
    = new G4Box("worldBox",10.*m,3.*m,10.*m);
  auto worldLogical
    = new G4LogicalVolume(worldSolid,vacuum,"worldLogical");
  auto worldPhysical
    = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"worldPhysical",0,
        false,0,checkOverlaps);

  // target 
  auto target_size_x = 50.*mm;
  auto target_size_y = 50.*mm;
  auto target_thickness = 0.5*mm; 
  auto targetSolid 
    = new G4Box("targetBox",target_size_x/2.,target_size_y/2.,target_thickness/2.);
  auto targetLogical
    = new G4LogicalVolume(targetSolid,vacuum,"targetLogical");
  //auto targetPhysical
  //  = new G4PVPlacement(0,G4ThreeVector(),targetLogical,"targetPhysical",
  //      worldLogical,false,0,checkOverlaps);

  // drift chamber (in)
  auto dc_size_r = 10.*cm;
  auto dc_thickness = 1.*mm;
  auto dcin_position = G4ThreeVector(0.,0.,0.);
  auto dcin_solid 
    = new G4Sphere("dcin_sphere",dc_size_r-dc_thickness/2.,dc_size_r+dc_thickness/2.,0.*degree,180.*degree,1.*degree,179.*degree);
  auto dcin_logical
    = new G4LogicalVolume(dcin_solid,air,"dcin_logical");
  auto dcin_physical
    = new G4PVPlacement(0,dcin_position,dcin_logical,"dcin_physical",
        worldLogical,false,0,checkOverlaps);
  // wireplane
  auto dc_wireplane_thickness = 1.*nm;
  auto dcin_wireplane_solid
    = new G4Sphere("dcin_wireplane_sphere", dc_size_r-dc_wireplane_thickness/2., dc_size_r+dc_wireplane_thickness/2., 0.*degree,180.*degree,1.*degree,179.*degree);
  dcin_wireplane_logical_
    = new G4LogicalVolume(dcin_wireplane_solid,argonGas,"dcin_wireplane_logical");
  auto dcin_wireplane_physical
    = new G4PVPlacement(0,G4ThreeVector(),dcin_wireplane_logical_,"dcin_wireplane_physical",
        dcin_logical, false,0,checkOverlaps);

  // drift chamber (out)
  auto dcout_position = G4ThreeVector(0.,0.,0.);
  G4RotationMatrix* dcout_rotation = new G4RotationMatrix();
  dcout_rotation->rotateZ(180.*degree);
  auto dcout_solid 
    = new G4Sphere("dcout_sphere",dc_size_r-dc_thickness/2.,dc_size_r+dc_thickness/2.,0.*degree,180.*degree,1.*degree,179.*degree);
  auto dcout_logical
    = new G4LogicalVolume(dcout_solid,air,"dcout_logical");
  auto dcout_physical
    = new G4PVPlacement(dcout_rotation,dcout_position,dcout_logical,"dcout_physical",
        worldLogical,false,0,checkOverlaps);
  // wireplane
  auto dcout_wireplane_solid
    = new G4Sphere("dcout_wireplane_sphere", dc_size_r-dc_wireplane_thickness/2., dc_size_r+dc_wireplane_thickness/2., 0.*degree,180.*degree,1.*degree,179.*degree);
  dcout_wireplane_logical_
    = new G4LogicalVolume(dcout_wireplane_solid,argonGas,"dcout_wireplane_logical");
  auto dcout_wireplane_physical
    = new G4PVPlacement(0,G4ThreeVector(),dcout_wireplane_logical_,"dcout_wireplane_physical",
        dcout_logical, false,0,checkOverlaps);

  // visualization attributes ------------------------------------------------

  auto visAttributes = new G4VisAttributes(G4Colour::White());
  visAttributes->SetVisibility(false);
  worldLogical->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  visAttributes = new G4VisAttributes(G4Colour::Blue());
  targetLogical->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  visAttributes = new G4VisAttributes(G4Colour::Yellow());
  dcin_logical->SetVisAttributes(visAttributes);
  dcout_logical->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);

  visAttributes = new G4VisAttributes(G4Colour::Green());
  dcin_wireplane_logical_->SetVisAttributes(visAttributes);
  dcout_wireplane_logical_->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);


  // return the world physical volume ----------------------------------------

  return worldPhysical;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  auto sdManager = G4SDManager::GetSDMpointer();
  G4String SDname;

  // sensitive detectors -----------------------------------------------------
  auto dcin = new DriftChamberSD(SDname="/dcin");
  sdManager->AddNewDetector(dcin);
  dcin_wireplane_logical_->SetSensitiveDetector(dcin);

  auto dcout = new DriftChamberSD(SDname="/dcout");
  sdManager->AddNewDetector(dcout);
  dcout_wireplane_logical_->SetSensitiveDetector(dcout);

}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructMaterials()
{
  auto nistManager = G4NistManager::Instance();

  // Air 
  nistManager->FindOrBuildMaterial("G4_AIR");

  // Argon gas
  nistManager->FindOrBuildMaterial("G4_Ar");

  // Scintillator
  // (PolyVinylToluene, C_9H_10)
  nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

  // CsI
  nistManager->FindOrBuildMaterial("G4_CESIUM_IODIDE");

  // Lead
  nistManager->FindOrBuildMaterial("G4_Pb");

  // Carbon
  nistManager->FindOrBuildMaterial("G4_C");

  // Vacuum "Galactic"
  nistManager->FindOrBuildMaterial("G4_Galactic");


  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
