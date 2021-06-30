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
/// \file B4DetectorConstruction.cc
/// \brief Implementation of the B4DetectorConstruction class

#include "B4DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4VSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal
    G4GlobalMagFieldMessenger *B4DetectorConstruction::fMagFieldMessenger = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
B4DetectorConstruction::B4DetectorConstruction()
    : G4VUserDetectorConstruction(),
      fAbsorberPV(nullptr),
      fGapPV(nullptr),
      fCheckOverlaps(true)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4DetectorConstruction::~B4DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *B4DetectorConstruction::Construct()
{
    // Define materials
    DefineMaterials();

    // Define volumes
    return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::DefineMaterials()
{
    // Air defined using NIST Manager
    auto nist = G4NistManager::Instance();
    nist->FindOrBuildMaterial("G4_AIR");

    // Aluminium material
    G4double a; // mass of a mole;
    G4double z; // z=mean number of protons;
    G4double density;
    new G4Material("Aluminium", z = 13., a = 26.981538 * g / mole, density = 2.702 * g / cm3);
    // The aluminium casing

    // Vacuum
    new G4Material("Galactic", z = 1., a = 1.01 * g / mole, density = universe_mean_density,
                   kStateGas, 2.73 * kelvin, 3.e-18 * pascal);

    // Print materials
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *B4DetectorConstruction::DefineVolumes()
{

    // Get materials
    // auto nist = G4NistManager::Instance();
    auto vaccum = G4Material::GetMaterial("Galactic");
    auto world_mat = G4Material::GetMaterial("G4_AIR");
    auto Al = G4Material::GetMaterial("Aluminium");

    if (!vaccum || !world_mat || !Al)
    {
        G4ExceptionDescription msg;
        msg << "Cannot retrieve materials already defined.";
        G4Exception("B4DetectorConstruction::DefineVolumes()",
                    "MyCode0001", FatalException, msg);
    }

    //
    // World
    //

    G4double world_sizeXY = 20 * cm;
    G4double world_sizeZ = 30 * cm;

    G4Box *solidWorld =
        new G4Box("World",                                                    //its name
                  0.5 * world_sizeXY, 0.5 * world_sizeXY, 0.5 * world_sizeZ); //its size

    G4LogicalVolume *worldLog =
        new G4LogicalVolume(solidWorld, //its solid
                            world_mat,  //its material
                            "World");   //its name

    G4VPhysicalVolume *physWorld =
        new G4PVPlacement(0,               //no rotation
                          G4ThreeVector(), //at (0,0,0)
                          worldLog,        //its logical volume
                          "World",         //its name
                          0,               //its mother  volume
                          false,           //no boolean operation
                          0,               //copy number
                          fCheckOverlaps); //overlaps checking

    // //
    // // Envelope
    // //
    // G4Box *solidEnv =
    //     new G4Box("Envelope",                                           //its name
    //               0.5 * env_sizeXY, 0.5 * env_sizeXY, 0.5 * env_sizeZ); //its size

    // G4LogicalVolume *logicEnv =
    //     new G4LogicalVolume(solidEnv,    //its solid
    //                         env_mat,     //its material
    //                         "Envelope"); //its name

    // new G4PVPlacement(0,               //no rotation
    //                   G4ThreeVector(), //at (0,0,0)
    //                   logicEnv,        //its logical volume
    //                   "Envelope",      //its name
    //                   worldLog,        //its mother  volume
    //                   false,           //no boolean operation
    //                   0,               //copy number
    //                   checkOverlaps);  //overlaps checking

    //Aluminium casing======
    G4double caseRi = 3.55 * cm;  // inner radius
    G4double caseRo = 3.65 * cm;  //outer radius
    G4double caseHz = 2.1 * cm;   // half height
    G4double caseAo = 0. * deg;   //start angle
    G4double caseAs = 360. * deg; //spanning angle
    G4Tubs *caseTube = new G4Tubs("Case",
                                  caseRi,
                                  caseRo,
                                  caseHz,
                                  caseAo,
                                  caseAs);
    G4LogicalVolume *caseLog = new G4LogicalVolume(caseTube, Al, "caseLog", 0, 0, 0);
    //position of the cylinder
    G4double cposx = 0. * cm;
    G4double cposy = 0. * cm;
    G4double cposz = 0. * cm;
    new G4PVPlacement(0, G4ThreeVector(cposx, cposy, cposz), caseLog, "case_vol", worldLog, false, 0, fCheckOverlaps);
    /////
    //window and the other side of aluminium case
    G4double winRi = 0. * cm;
    G4double winRo = caseRo;
    G4double winHz = 0.1 * cm;
    G4double winAo = 0. * deg;
    G4double winAs = 360. * deg;
    G4Tubs *window = new G4Tubs("insideCase",
                                winRi,
                                winRo,
                                winHz,
                                winAo,
                                winAs);
    G4LogicalVolume *winLog = new G4LogicalVolume(window, Al, "winLog", 0, 0, 0);
    //position of the window of Al...
    G4double win1posx = 0. * cm;
    G4double win1posy = 0. * cm;
    G4double win1posz = caseHz;
    fGapPV = new G4PVPlacement(0,
                               G4ThreeVector(win1posx, win1posy, win1posz),
                               winLog,
                               "windowPhysical",
                               caseLog,
                               false,
                               0,
                               fCheckOverlaps);
    //position of the other side of the window...
    G4double win2posx = 0. * cm;
    G4double win2posy = 0. * cm;
    G4double win2posz = -caseHz;
    new G4PVPlacement(0,
                      G4ThreeVector(win2posx, win2posy, win2posz),
                      winLog,
                      "othersidePhysical",
                      caseLog,
                      false,
                      0,
                      fCheckOverlaps);
    ///// Aluminium Casing.... DONE../

    //Inner Vaccum inside the aluminium
    G4double ivacRi = 0. * cm;
    G4double ivacRo = caseRi;
    G4double ivacHz = caseHz - winHz;
    G4double ivacAo = 0. * deg;
    G4double ivacAs = 360. * deg;
    G4Tubs *insideCase = new G4Tubs("insideCase",
                                    ivacRi,
                                    ivacRo,
                                    ivacHz,
                                    ivacAo,
                                    ivacAs);
    G4LogicalVolume *insideLog = new G4LogicalVolume(insideCase, vaccum, "insedeCaseLog", 0, 0, 0);
    //position of inner vaccum...
    G4double ivacposx = 0. * cm;
    G4double ivacposy = 0. * cm;
    G4double ivacposz = 0. * cm;
    new G4PVPlacement(0,
                      G4ThreeVector(ivacposx, ivacposy, ivacposz),
                      insideLog,
                      "inside_vol",
                      caseLog,
                      false,
                      0,
                      fCheckOverlaps);

    //// Hpge tube
    // n - type HPGe crystal creation
    G4String name, symbol;
    G4int ncomponents;
    G4double z, a, fractionmass, density;
    // z = 32;
    a = 72.64 * g / mole;
    density = 5.323 * g / cm3;
    G4Element *pureGe = new G4Element(name = "Germanium", symbol = "Ge", z = 32, a);
    // z = 15;
    a = 30.97 * g / mole;
    G4Element *P = new G4Element("Phosphorus", symbol = "P", z = 15, a);
    G4Material *Ge = new G4Material(name = "nHPGe", density, ncomponents = 2);
    Ge->AddElement(pureGe, fractionmass = 99. * perCent);
    Ge->AddElement(P, fractionmass = 1. * perCent);

    G4double hpgeRi = 0. * cm;
    G4double hpgeRo = 3.05 * cm;
    G4double hpgeHz = 1.6 * cm;
    G4double hpgeAo = 0. * deg;
    G4double hpgeAs = 360. * deg;

    G4VSolid *hpgeTube = new G4Tubs("hpgeTube",
                                    hpgeRi,
                                    hpgeRo,
                                    hpgeHz,
                                    hpgeAo,
                                    hpgeAs);

    ///////
    //Bore Hole========
    G4double boreRi = 0. * mm;
    G4double boreRo = 7.072 * mm;
    G4double boreHz = 6.338 * mm;
    G4double boreAo = 0. * deg;
    G4double boreAs = 360. * deg;

    G4VSolid *boreTube = new G4Tubs("boreTube",
                                    boreRi,
                                    boreRo,
                                    boreHz,
                                    boreAo,
                                    boreAs);
    //position of the hole wrt hpge....
    G4double boreposx = 0. * cm;
    G4double boreposy = 0. * cm;
    G4double boreposz = -(hpgeHz - boreHz);
    G4VSolid *hpgeWithBore = new G4SubtractionSolid("hpgeWithBore",                               //name of the solid
                                                    hpgeTube,                                     //parent volume
                                                    boreTube,                                     //substractor volume
                                                    0,                                            // spin of subtractor w.r.t parent
                                                    G4ThreeVector(boreposx, boreposy, boreposz)); //position w.r.t. parent
    G4double taperRi = 4.24264 * cm;                                                              // <-- it is much sensitive towars the tapering
    G4double taperRo = 4.74 * cm;                                                                 // <-- this is taken 10% for more safety
    G4double taperPhiO = 0. * deg;
    G4double taperPhiS = 360. * deg;
    //upper tapering
    G4double taperUThetaO = 0 * deg;    //44.805 * deg;
    G4double taperUThetaS = 48.0 * deg; //2.34 * deg;

    G4VSolid *taperU = new G4Sphere("tapering",                  //Name
                                    taperRi, taperRo,            //inner and outer radius
                                    taperPhiO, taperPhiS,        // starting and spanning angle phi
                                    taperUThetaO, taperUThetaS); //starting and spanning angle Theta
    // position of the upper taper
    G4double taperUposx = 0. * cm;
    G4double taperUposy = 0. * cm;
    G4double taperUposz = (hpgeHz - hpgeRo);
    G4VSolid *hpgeWithUTaper = new G4SubtractionSolid("hpgeWithUTaper",                                   //name of the solid
                                                      hpgeWithBore,                                       //parent volume
                                                      taperU,                                             //substractor volume
                                                      0,                                                  // spin of subtractor w.r.t parent
                                                      G4ThreeVector(taperUposx, taperUposy, taperUposz)); //position w.r.t. parent

    //Lower tapering...
    G4double taperLThetaO = 90.0 * deg;                          //134.805 * deg;
    G4double taperLThetaS = 48.0 * deg;                          //2.34 * deg;
    G4VSolid *taperL = new G4Sphere("tapering",                  //Name
                                    taperRi, taperRo,            //inner and outer radius
                                    taperPhiO, taperPhiS,        // starting and spanning angle phi
                                    taperLThetaO, taperLThetaS); //starting and spanning angle Theta
    G4double taperLposx = 0. * cm;
    G4double taperLposy = 0. * cm;
    G4double taperLposz = -(hpgeHz - hpgeRo);
    G4VSolid *hpgeWithTapers = new G4SubtractionSolid("hpgeWithTapers",                                   //name of the solid
                                                      hpgeWithUTaper,                                     //parent volume
                                                      taperL,                                             //substractor volume
                                                      0,                                                  // spin of subtractor w.r.t parent
                                                      G4ThreeVector(taperLposx, taperLposy, taperLposz)); //position w.r.t. parent
    G4double taperBRi = 10.001 * mm;
    G4double taperBRo = 10.7084 * mm;
    G4double taperBThetaO = 0 * deg;  //40.5516 * deg;
    G4double taperBThetaS = 50 * deg; //8.0883 * deg;
    G4double taperBPhiO = 0. * deg;
    G4double taperBPhiS = 360. * deg;

    G4VSolid *taperB = new G4Sphere("taperingBore",              //Name
                                    taperBRi, taperBRo,          //inner and outer radius
                                    taperBPhiO, taperBPhiS,      // starting and spanning angle phi
                                    taperBThetaO, taperBThetaS); //starting and spanning angle Theta
    //position of bore tapering
    G4double taperBposx = 0. * cm;
    G4double taperBposy = 0. * cm;
    G4double taperBposz = -(hpgeHz + boreRo);

    G4VSolid *hpgeAllTaper = new G4SubtractionSolid("hpgeAllTaper",                                     //name of the solid
                                                    hpgeWithTapers,                                     //parent volume
                                                    taperB,                                             //substractor volume
                                                    0,                                                  // spin of subtractor w.r.t parent
                                                    G4ThreeVector(taperBposx, taperBposy, taperBposz)); //position w.r.t. parent

    //Point contact dimension...
    G4double PCRi = 0. * cm;
    G4double PCRo = boreRo;
    G4double PCHz = 0.05 * mm;
    G4double PCAo = 0 * deg;
    G4double PCAs = 360 * deg;
    //position of point contact crystal
    G4double PCx = 0. * cm;
    G4double PCy = 0. * cm;
    G4double PCz = -(hpgeHz - (2 * boreRo) - PCHz);

    G4VSolid *PCTube = new G4Tubs("poinContactTube",
                                  PCRi,
                                  PCRo,
                                  PCHz,
                                  PCAo,
                                  PCAs);
    //cutting out the main crystal for point contact

    G4VSolid *hpgeActual = new G4SubtractionSolid("hpgeActual",                  //name of the solid
                                                  hpgeAllTaper,                  //parent volume
                                                  PCTube,                        //substractor volume
                                                  0,                             // spin of subtractor w.r.t parent
                                                  G4ThreeVector(PCx, PCy, PCz)); //position w.r.t. parent

    /////////// n -type crystal cutting.... DONE../
    // p - type HPGe crystal formation...
    // z = 5;
    a = 10.81 * g / mole;
    G4Element *B = new G4Element(name = "Boron", symbol = "B", z = 5, a);
    G4Material *pGe = new G4Material(name = "pHPGe", density, ncomponents = 2);
    pGe->AddElement(pureGe, fractionmass = 97. * perCent);
    pGe->AddElement(B, fractionmass = 3. * perCent);

    G4LogicalVolume *pGeLog = new G4LogicalVolume(PCTube, pGe, "pTypeHpgeLog", 0, 0, 0);
    //
    G4LogicalVolume *hpgeLog = new G4LogicalVolume(hpgeActual, Ge, "hpgeLog", 0, 0, 0);
    //position of solid hpge...
    G4double hpgeposx = 0. * m;
    G4double hpgeposy = 0. * m;
    G4double hpgeposz = 0. * m;
    fAbsorberPV = new G4PVPlacement(0, G4ThreeVector(hpgeposx, hpgeposy, hpgeposz), hpgeLog, "hpge", worldLog, false, 0, fCheckOverlaps);

    //position of point contact w.r.t. HPGe...
    new G4PVPlacement(0, G4ThreeVector(PCx, PCy, PCz), pGeLog, "hpge", hpgeLog, false, 0, fCheckOverlaps);

    //
    // Visualization attributes
    //
    worldLog->SetVisAttributes(G4VisAttributes::GetInvisible());

    auto simpleBoxVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    simpleBoxVisAtt->SetVisibility(true);
    hpgeLog->SetVisAttributes(simpleBoxVisAtt);

    /*
  auto worldS = new G4Box("World",                                           // its name
                          worldSizeXY / 2, worldSizeXY / 2, worldSizeZ / 2); // its size

  auto worldLV = new G4LogicalVolume(
      worldS,          // its solid
      defaultMaterial, // its material
      "World");        // its name

  auto worldPV = new G4PVPlacement(
      0,               // no rotation
      G4ThreeVector(), // at (0,0,0)
      worldLV,         // its logical volume
      "World",         // its name
      0,               // its mother  volume
      false,           // no boolean operation
      0,               // copy number
      fCheckOverlaps); // checking overlaps

  //
  // Calorimeter
  //
  auto calorimeterS = new G4Box("Calorimeter",                                         // its name
                                calorSizeXY / 2, calorSizeXY / 2, calorThickness / 2); // its size

  auto calorLV = new G4LogicalVolume(
      calorimeterS,    // its solid
      defaultMaterial, // its material
      "Calorimeter");  // its name

  new G4PVPlacement(
      0,               // no rotation
      G4ThreeVector(), // at (0,0,0)
      calorLV,         // its logical volume
      "Calorimeter",   // its name
      worldLV,         // its mother  volume
      false,           // no boolean operation
      0,               // copy number
      fCheckOverlaps); // checking overlaps

  //
  // Layer
  //
  auto layerS = new G4Box("Layer",                                               // its name
                          calorSizeXY / 2, calorSizeXY / 2, layerThickness / 2); // its size

  auto layerLV = new G4LogicalVolume(
      layerS,          // its solid
      defaultMaterial, // its material
      "Layer");        // its name

  new G4PVReplica(
      "Layer",         // its name
      layerLV,         // its logical volume
      calorLV,         // its mother
      kZAxis,          // axis of replication
      nofLayers,       // number of replica
      layerThickness); // witdth of replica

  //
  // Absorber
  //
  auto absorberS = new G4Box("Abso",                                               // its name
                             calorSizeXY / 2, calorSizeXY / 2, absoThickness / 2); // its size

  auto absorberLV = new G4LogicalVolume(
      absorberS,        // its solid
      absorberMaterial, // its material
      "Abso");          // its name

  fAbsorberPV = new G4PVPlacement(
      0,                                        // no rotation
      G4ThreeVector(0., 0., -gapThickness / 2), // its position
      absorberLV,                               // its logical volume
      "Abso",                                   // its name
      layerLV,                                  // its mother  volume
      false,                                    // no boolean operation
      0,                                        // copy number
      fCheckOverlaps);                          // checking overlaps

  //
  // Gap
  //
  auto gapS = new G4Box("Gap",                                               // its name
                        calorSizeXY / 2, calorSizeXY / 2, gapThickness / 2); // its size

  auto gapLV = new G4LogicalVolume(
      gapS,        // its solid
      gapMaterial, // its material
      "Gap");      // its name

  fGapPV = new G4PVPlacement(
      0,                                        // no rotation
      G4ThreeVector(0., 0., absoThickness / 2), // its position
      gapLV,                                    // its logical volume
      "Gap",                                    // its name
      layerLV,                                  // its mother  volume
      false,                                    // no boolean operation
      0,                                        // copy number
      fCheckOverlaps);                          // checking overlaps

  //
  // print parameters
  //
  G4cout
      << G4endl
      << "------------------------------------------------------------" << G4endl
      << "---> The calorimeter is " << nofLayers << " layers of: [ "
      << absoThickness / mm << "mm of " << absorberMaterial->GetName()
      << " + "
      << gapThickness / mm << "mm of " << gapMaterial->GetName() << " ] " << G4endl
      << "------------------------------------------------------------" << G4endl;

  //
  // Visualization attributes
  //
  worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

  auto simpleBoxVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
  simpleBoxVisAtt->SetVisibility(true);
  calorLV->SetVisAttributes(simpleBoxVisAtt);
*/
    //
    // Always return the physical World
    //
    return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::ConstructSDandField()
{
    // Create global magnetic field messenger.
    // Uniform magnetic field is then created automatically if
    // the field value is not zero.
    G4ThreeVector fieldValue;
    fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
    fMagFieldMessenger->SetVerboseLevel(1);

    // Register the field messenger for deleting
    G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
