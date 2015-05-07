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
// $Id: CMSHFDetectorConstruction.cc,v 1.5 2013/06/18 16:08:49 cowden Exp $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <iostream>

#include "CMSHFDetectorConstruction.hh"
#include "HFStackingAction.hh"
#include "HFPrimaryGeneratorAction.hh"

#include "PMTQE.h"

#include "G4RunManager.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4NistManager.hh"
#include "G4UserLimits.hh"
#include "G4UniformMagField.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"


#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CMSHFDetectorConstruction::CMSHFDetectorConstruction()
:m_isConstructed(false),m_nFib(1.457),m_nClad(1.42)
,m_absFib(28.1*m),m_absClad(28.1*m)
//,m_nSFib(1.59),m_nSClad(1.49)
,m_nSFib(1.457),m_nSClad(1.42)
//,m_absSFib(4.*m),m_absSClad(4.*m)
,m_absSFib(28.1*m),m_absSClad(28.1*m)
,m_nGlass(1.517),m_absGlass(1.*m)
,m_Bfield(NULL)
,m_stacking(NULL)
//,m_gun(NULL)
{

  m_expHall_z = 5.0*m;
  m_expHall_x = 5.0*m;
  m_expHall_y = 5.0*m;

  m_pitch = 0.*deg;
  m_yaw = 0.*deg;

  m_length = 20.*cm;  // half length of the calorimeter
  m_segWidth = 10*cm;
  m_segHeight = 20*cm;
  m_absDepth = 5.*cm;

  // set the position of the detector
  //m_zPos = 7.3*m;
  m_zPos = 0.6*m*cos(m_yaw);
  m_xPos = 0.*m;
  m_yPos = 0.*m;

  // This is the full radius including the cladding
  // since the cladding is a daughter volume of the fiber 
  m_rCFib = 0.400*mm;
  m_rCClad = 0.330*mm;
  m_rCCore = 0.3*mm;

  m_rS1Fib = 0.175*mm;
  m_rS1Clad = 0.115*mm;
  m_rS1Glass = 0.100*mm;
  m_rS1Core = 0.030*mm;

  m_rS2Fib = 0.275*mm;
  m_rS2Clad = 0.215*mm;
  m_rS2Glass = 0.200*mm;
  m_rS2Core = 0.075*mm;

  m_rair = 1.25*mm;

  m_Nseg = 4U;

  m_nQseg = 4;
  m_nSseg = 0;

  // Fiber coupling parameters
  m_rCollar = 1.27*cm;
  m_rOuterCollar = 1.77*cm;
  m_pmtWindowThickness = 5.*mm;
  m_collarLength = 2*cm;


  // dead material around block
  m_deadSeg_bottom = 1U;
  m_deadSeg_top = 1U;
  m_deadSeg_right = 1U;
  m_deadSeg_left = 1U;

  m_fillFibres = false; // by default don't insert fibers
  
  // scintillation properties
  m_scinFastConst = 1.*ns;
  m_scinSlowConst = 10.*ns;
  m_scinYield = 8300./MeV;
  m_scinYieldRatio = 0.8;

  // create null pointers to some physical volumes
  // this is an improtant check when geometry is refreshed as some volumes
  // have not yet been placed.
  m_collar_phys = NULL;
  m_pmtWindow_phys = NULL;
  m_pmtSpace_phys = NULL;

  m_qFibreScin1_log = NULL;
  m_buffScin1_log = NULL;
  m_jacketScin1_log = NULL;
  m_cladScin1_log = NULL;

  m_qFibreScin2_log = NULL;
  m_buffScin2_log = NULL;
  m_jacketScin2_log = NULL;
  m_cladScin2_log = NULL;

  // initialize fiber type
  m_fbtype = FBquartz;

  m_checkOverlaps = false;
  m_messenger = new CMSHFDetectorConstructionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CMSHFDetectorConstruction::~CMSHFDetectorConstruction()
{ 
  delete m_messenger; 
  //delete m_fibLimits;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* CMSHFDetectorConstruction::Construct()
{
  
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  if ( !m_isConstructed ) {
    m_isConstructed = true;
    CalculateConstants();
    DefineMaterials();
    SetupWorld();
    SetupGeometry();
    SetupDetectors();
  }

  // print the calo parameters
  PrintCalorParameters();


//always return the physical World
  return m_expHall_phys;
}


//  Define the materials used in the detector
void CMSHFDetectorConstruction::DefineMaterials()
{ 

  // ------------- Materials -------------

  G4double a, z, density;
  G4int nelements;

  // Air
  // 
  G4Element* N = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);
  G4Element* Si = new G4Element("Silicon", "Si", z=14, a=28.09*g/mole);

  G4NistManager* man = G4NistManager::Instance();
  G4Element* B = man->FindOrBuildElement("B");
  G4Element* Na = man->FindOrBuildElement("Na");
  G4Element* Al = man->FindOrBuildElement("Al");
  G4Element* K = man->FindOrBuildElement("K");
  G4Element* W = man->FindOrBuildElement("W");
  G4Element* Ni = man->FindOrBuildElement("Ni");
  G4Element* Cu = man->FindOrBuildElement("Cu");
  G4Element* C = man->FindOrBuildElement("C");
  G4Element* H = man->FindOrBuildElement("H");
  G4Element* Zn = man->FindOrBuildElement("Zn");

  // air (material)  optical properties are defined below
  m_air = new G4Material("Air", density=1.29*mg/cm3, nelements=2);
  m_air->AddElement(N, 70.*perCent);
  m_air->AddElement(O, 30.*perCent);

  // vacuum (same as air, but less dense) 
  m_vacuum = new G4Material("vacuum", density=0.129*mg/cm3,nelements=2);
  m_vacuum->AddElement(N,70.*perCent);
  m_vacuum->AddElement(O, 30.*perCent);
  

  m_buffer = man->FindOrBuildMaterial("G4_KAPTON");

  // Brass
  m_brass = new G4Material("Brass",8.525*g/cm3,2);
  m_brass->AddElement(Cu,70.*perCent);
  m_brass->AddElement(Zn,30.*perCent);

  // Iron
  m_iron = new G4Material("Iron", z=26, a=55.845*g/mole, density=7.874*g/cm3);

  // Lead
  m_lead = new G4Material("Lead", z=82, a=207.2*g/mole, density=11.35*g/cm3);

  // steel
  m_steel = man->FindOrBuildMaterial("G4_STAINLESS-STEEL");

  // collar material (using lucite since it's in the NIST database)
  m_collarMaterial = man->FindOrBuildMaterial("G4_LUCITE");

  // Quartz
  //
  m_quartz = new G4Material("quartz",2.2*g/cm3,2);
  m_quartz->AddElement(Si,1);
  m_quartz->AddElement(O,2);

  // cladding material
  m_cladCher = new G4Material("claddingCher",2.2*g/cm3,2);
  m_cladCher->AddElement(Si,1);
  m_cladCher->AddElement(O,2);

  // scsf
  m_scsf78 = new G4Material("scsf78",2.2*g/cm3,2);
  m_scsf78->AddElement(Si,1);
  m_scsf78->AddElement(O,2);
  // PS 
  //m_scsf78->AddElement(C,8);
  //m_scsf78->AddElement(H,8);

  // cladding material
  m_cladScin = new G4Material("claddingScin",2.2*g/cm3,2);
  m_cladScin->AddElement(Si,1);
  m_cladScin->AddElement(O,2);
  // PMMA 
  //m_cladScin->AddElement(C,5);
  //m_cladScin->AddElement(H,8);
  //m_cladScin->AddElement(O,2);

  //
  // ------------ Generate & Add Material Properties Table ------------
  //
  const G4int nEntries = 32;

  G4double PhotonEnergy[nEntries] =
            { 2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV,
              2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,
              2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,
              2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,
              2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV,
              3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV,
              3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,
              3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV };


  //
  // Air
  //
  G4double RefractiveIndex2[nEntries] =
            { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00 };

  G4double reflectivity[nEntries] =
            { 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
              0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
              0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
              0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
              0.00, 0.00, 0.00, 0.00 };

  G4double efficiency[nEntries] =
            { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00 };

  G4MaterialPropertiesTable* myMPT2 = new G4MaterialPropertiesTable();
  myMPT2->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex2, nEntries);
  
  m_air->SetMaterialPropertiesTable(myMPT2);

  //
  // vacuum optical properties

  G4double * pmtE = pmtEnergies();
  G4double * pmtEff = pmtEfficiencies();
  G4int pmtInt = pmtEntries;
  m_vacMPT = new G4MaterialPropertiesTable();
  m_vacMPT->AddProperty("REFLECTIVITY",PhotonEnergy, reflectivity, nEntries);
  m_vacMPT->AddProperty("EFFICIENCY", pmtE, pmtEff, pmtInt);
  m_vacuum->SetMaterialPropertiesTable(m_vacMPT);

  G4cout << "PMT Material " << pmtE << " " << pmtEff << " " << pmtEff[0] << G4endl;
  for ( unsigned i=0; i != pmtEntries; i++ ) 
    G4cout << "PMT Eff: " << pmtE[i] << " " << pmtEff[i] << G4endl;

  // Glass
  //
  m_glass = new G4Material("glass",2.23*g/cm3,6);
  m_glass->AddElement(B,0.040064);
  m_glass->AddElement(O,0.539562);
  m_glass->AddElement(Na,0.028191);
  m_glass->AddElement(Al,0.011644);
  m_glass->AddElement(Si,0.377220);
  m_glass->AddElement(K,0.003321);

  G4MaterialPropertiesTable *glassPT = new G4MaterialPropertiesTable();
  G4double glassE[3] = { 1.625*eV, 5.*eV, 12.4*eV };
  G4double glassN[3] = {m_nGlass,m_nGlass,m_nGlass};
  G4double glassAbs[3] = {m_absGlass,m_absGlass,m_absGlass};
  glassPT->AddProperty("RINDEX",glassE,glassN,3);
  glassPT->AddProperty("ABSLENGTH",glassE,glassAbs,3);

  m_glass->SetMaterialPropertiesTable(glassPT);


  // ------------- Surfaces and Optical Properties --------------


  const unsigned nEnergies = 3;
  G4double energies[nEnergies] = { 1.625*eV, 5.*eV, 12.4*eV };
  G4double qRindex[nEnergies] = {m_nFib, m_nFib, m_nFib};   
  G4double qAbsLength[nEnergies] = {m_absFib, m_absFib, m_absFib};

  G4MaterialPropertiesTable *qProps = new G4MaterialPropertiesTable();
  qProps->AddProperty("RINDEX",energies,qRindex,nEnergies);
  qProps->AddProperty("ABSLENGTH",energies,qAbsLength,nEnergies);
  m_quartz->SetMaterialPropertiesTable(qProps);

  G4double cRindex[nEnergies] = {m_nClad,m_nClad,m_nClad};
  G4double cAbsLength[nEnergies] = {m_absClad,m_absClad,m_absClad};
  
  G4MaterialPropertiesTable *cProps = new G4MaterialPropertiesTable();
  cProps->AddProperty("RINDEX",energies,cRindex,nEnergies);
  cProps->AddProperty("ABSLENGTH",energies,cAbsLength,nEnergies);
  m_cladCher->SetMaterialPropertiesTable(cProps);

  // set up the material properties for the scintillating fibers
  DefineScintillator();

}

// define the scintillator properties
void CMSHFDetectorConstruction::DefineScintillator()
{
  const unsigned nScinEnergies = 6;
  G4double scinEnergies[nScinEnergies] = {3.1*eV, 2.9*eV, 2.76*eV, 2.48*eV, 2.25*eV, 2.07*eV };
  G4double scinValues[nScinEnergies] = {0.0, 5.0, 10., 5.0, 2.0, 1.0};
  G4double scinRindex[nScinEnergies] = { m_nSFib, m_nSFib, m_nSFib, m_nSFib, m_nSFib, m_nSFib };
  G4double scinAbsLength[nScinEnergies] = { m_absSFib, m_absSFib, m_absSFib, m_absSFib, m_absSFib, m_absSFib };
  G4double scinCladRindex[nScinEnergies] = { m_nSClad, m_nSClad, m_nSClad, m_nSClad, m_nSClad, m_nSClad };
  G4double scinCladAbsLength[nScinEnergies] = { m_absSClad, m_absSClad, m_absSClad, m_absSClad, m_absSClad, m_absSClad };


  G4MaterialPropertiesTable *scsfProps = new G4MaterialPropertiesTable();
  scsfProps->AddProperty("RINDEX",scinEnergies,scinRindex,nScinEnergies);
  scsfProps->AddProperty("ABSLENGTH",scinEnergies,scinAbsLength,nScinEnergies);

  /*scsfProps->AddProperty("FASTCOMPONENT",scinEnergies,scinValues,nScinEnergies);
  scsfProps->AddProperty("SLOWCOMPONENT",scinEnergies,scinValues,nScinEnergies);

  //
  // see http://infoscience.epfl.ch/record/164027/files/EPFL_TH5033.pdf
  // for scintillation yield
  scsfProps->AddConstProperty("SCINTILLATIONYIELD", m_scinYield);
  scsfProps->AddConstProperty("RESOLUTIONSCALE", 2.0);
  scsfProps->AddConstProperty("FASTTIMECONSTANT", m_scinFastConst);
  scsfProps->AddConstProperty("SLOWTIMECONSTANT", m_scinSlowConst);
  scsfProps->AddConstProperty("YIELDRATIO", m_scinYieldRatio);*/

  G4MaterialPropertiesTable *scsfCladProps = new G4MaterialPropertiesTable();
  scsfCladProps->AddProperty("RINDEX",scinEnergies,scinCladRindex,nScinEnergies);
  scsfCladProps->AddProperty("ABSLENGTH",scinEnergies,scinCladAbsLength,nScinEnergies);

  G4MaterialPropertiesTable *oldProps = m_scsf78->GetMaterialPropertiesTable();

  //m_scsf78->SetMaterialPropertiesTable(scsfProps);
  m_scsf78->SetMaterialPropertiesTable(scsfProps);
  if ( oldProps ) delete oldProps;

  oldProps = m_cladScin->GetMaterialPropertiesTable();
  m_cladScin->SetMaterialPropertiesTable(scsfCladProps);
  if ( oldProps ) delete oldProps;

}


// Setup the world geometry
void CMSHFDetectorConstruction::SetupWorld()
{

  // The experimental Hall
  //
  m_expHall_box = new G4Box("World",m_expHall_x,m_expHall_y,m_expHall_z);

  m_expHall_log
    = new G4LogicalVolume(m_expHall_box,m_air,"World",0,0,0);

  m_expHall_phys
    = new G4PVPlacement(0,G4ThreeVector(0, 0, 0),m_expHall_log,"World",0,false,0);

}



// Setup the logical volumes
void CMSHFDetectorConstruction::SetupGeometry()
{ 

  //
  // ------------- primitives -----------
  G4cout << "CMSHF Constructing block: " << 2.*m_Wdx << "x" << 2.*m_Wdy << "x" << 2.*m_length << G4endl;
  m_absBlock = new G4Box("block",m_segWidth/2.,m_segHeight/2.,m_absDepth/2.);

  m_qFibreCher = new G4Tubs("quartzFibre",0.,m_rCFib,m_length,0.,2.*pi);
  m_cladCher_tube = new G4Tubs("cladding",m_rCCore,m_rCClad,m_length,0,2.*pi);
  m_buffCher_tube = new G4Tubs("buffer",m_rCClad,m_rCFib,m_length,0,2.*pi);

  m_qFibreScin1 = new G4Tubs("scsfFibre",0.,m_rS1Fib,m_length,0.,2.*pi);
  m_cladScin1_tube = new G4Tubs("cladding",m_rS1Glass,m_rS1Clad,m_length,0,2.*pi);
  m_buffScin1_tube = new G4Tubs("buffer",m_rS1Clad,m_rS1Fib,m_length,0,2.*pi);
  m_jacketScin1_tube = new G4Tubs("jacket",m_rS1Core,m_rS1Glass,m_length,0,2.*pi);

  m_qFibreScin2 = new G4Tubs("scsfFibre",0.,m_rS2Fib,m_length,0.,2.*pi);
  m_cladScin2_tube = new G4Tubs("cladding",m_rS2Glass,m_rS2Clad,m_length,0,2.*pi);
  m_buffScin2_tube = new G4Tubs("buffer",m_rS2Clad,m_rS2Fib,m_length,0,2.*pi);
  m_jacketScin2_tube = new G4Tubs("jacket",m_rS2Core,m_rS2Glass,m_length,0,2.*pi);

  m_glass_box = new G4Box("Glass",m_segWidth/2.,m_segHeight/2.,1.*cm);

  m_collar_tube = new G4Tubs("collar",m_rCollar,m_rOuterCollar,m_collarLength/2.,0.,2.*pi);
  m_pmtWindow_plate = new G4Tubs("pmtWindow",0.,m_rCollar,m_pmtWindowThickness/2.,0.,2.*pi);

  //
  // ------------- Volumes --------------
  const unsigned segTot = m_Nseg;

  m_absBlock_log = new G4LogicalVolume(m_absBlock,m_iron,"block",0,0,0);
  m_pbBlock_log = new G4LogicalVolume(m_absBlock,m_lead,"lead",0,0,0);
  // air gap 

  
  // Quartz fibre
  m_qFibreCher_log = new G4LogicalVolume(m_qFibreCher,m_quartz,"quartzFibreLog",0,0,0);
  
  // Cladding on fibre
  m_cladCher_log = new G4LogicalVolume(m_cladCher_tube,m_cladCher,"claddingCher",0,0,0);

  // buffer on Cherenkov Fibre
  m_buffCher_log = new G4LogicalVolume(m_buffCher_tube,m_buffer,"bufferCherLog",0,0,0);
  
  // Scintillating scsf-78
  m_qFibreScin1_log = new G4LogicalVolume(m_qFibreScin1,m_scsf78,"scsf78FibreLog",0,0,0);

  // jacket
  m_jacketScin1_log = new G4LogicalVolume(m_jacketScin1_tube,m_quartz,"jacket",0,0,0);
  
  // Cladding on fibre
  m_cladScin1_log = new G4LogicalVolume(m_cladScin1_tube,m_cladScin,"claddingScin",0,0,0);

  // buffer on Scintillating Fibre
  m_buffScin1_log = new G4LogicalVolume(m_buffScin1_tube,m_buffer,"bufferScinLog",0,0,0);

  // Scintillating scsf-78
  m_qFibreScin2_log = new G4LogicalVolume(m_qFibreScin2,m_scsf78,"scsf78FibreLog",0,0,0);

  // jacket
  m_jacketScin2_log = new G4LogicalVolume(m_jacketScin2_tube,m_quartz,"jacket",0,0,0);
  
  // Cladding on fibre
  m_cladScin2_log = new G4LogicalVolume(m_cladScin2_tube,m_cladScin,"claddingScin",0,0,0);

  // buffer on Scintillating Fibre
  m_buffScin2_log = new G4LogicalVolume(m_buffScin2_tube,m_buffer,"bufferScinLog",0,0,0);

  // Glass 
  //
  m_glass_log = new G4LogicalVolume(m_glass_box,m_glass,"glassLog",0,0,0); 

  // collar and pmt window
  m_collar_log = new G4LogicalVolume(m_collar_tube,m_collarMaterial,"collar",0,0,0);
  m_pmtWindow_log = new G4LogicalVolume(m_pmtWindow_plate,m_glass,"pmtWindow",0,0,0);
  m_pmtSpace_log = new G4LogicalVolume(m_pmtWindow_plate,m_vacuum,"pmtSpace",0,0,0);

  // set limits on the time to propagate photons
  // needs a G4StepLimiter needs to be added to OpticalPhoton
  // process for the following to do anything.
  /*m_fibLimits  = new G4UserLimits();
  m_fibLimits->SetUserMaxTime(0.0*ns);
  m_qFibre_log->SetUserLimits(m_fibLimits);
  m_clad_log->SetUserLimits(m_fibLimits);
  m_glass_log->SetUserLimits(m_fibLimits);
  m_expHall_log->SetUserLimits(m_fibLimits);*/

}

// Setup the detectors physical volumes
void CMSHFDetectorConstruction::SetupDetectors()
{ 

  int count=0,scsfCount=0;

  const unsigned segTot = m_Nseg;
  m_pbBlock_phys.resize(segTot);

  // setup block and wedge transformation
  G4RotationMatrix rot;
  rot.rotateX(m_pitch);
  rot.rotateY(m_yaw);

  G4ThreeVector glassOffset = rot(G4ThreeVector(0.,0.,m_length+1.*cm));

  // set the fibre orientation in stacking action
  if ( m_stacking ) {
    m_stacking->SetFibreDirection(rot(G4ThreeVector(0.,0.,1.)));
    m_stacking->SetFibLength(m_length*2.);
  }

  // place first absorber segment
  m_absBlock_phys = new G4PVPlacement(G4Transform3D(rot,G4ThreeVector(0.,0.,m_absDepth/2.)),m_absBlock_log,"absorber",m_expHall_log,false,0,m_checkOverlaps);
  
  m_glass_phys = new G4PVPlacement(G4Transform3D(rot,m_segPositions[0]+glassOffset),m_glass_log,"glass",m_expHall_log,false,0,m_checkOverlaps); 

  // ------------------------------------------------------
  // insert fibers 
  char name[50];
  unsigned segRow=0,segCol=0;
  unsigned fibCountRow=0; // count of fibres (total) in a row across all segments
  unsigned rowCount=0; // count of the number of rows across all segments
  const unsigned nFibSide = m_Nseg;

  if ( m_fillFibres ) PlaceFibres();


  // place the lead blocks
  for ( unsigned i=0; i != segTot; i++ ) {
    m_pbBlock_phys[i] = new G4PVPlacement(G4Transform3D(rot,m_segPositions[i]),m_pbBlock_log,"absorber",m_expHall_log,false,i,m_checkOverlaps);
  }


}

void CMSHFDetectorConstruction::PlaceFibres()
{

  G4RotationMatrix rot;
  rot.rotateY(90*deg);

  // place two quartz fibres
  //m_fibresCher.push_back(new G4PVPlacement(G4Transform3D(rot,G4ThreeVector(0.,0.,10.*cm)),m_qFibreCher_log,"Cfib0",m_expHall_log,false,0,m_checkOverlaps));

  if ( m_fbtype == FBquartz ) {
    // place quartz fibers in 10x10 grid
    const unsigned nSide = 10U;
    for ( unsigned i=0; i != nSide; i++ ) {
      const double yPos = 9.*m_rCFib - i*2.*m_rCFib;
      for ( unsigned j=0; j != nSide; j++ ) {
        const double zPos = 10.*cm + 9.*m_rCFib - j*2.*m_rCFib;
        char name[50];
        sprintf(name,"Cfib%d",i*nSide+j);
        m_fibresCher.push_back(new G4PVPlacement(G4Transform3D(rot,G4ThreeVector(0.,yPos,zPos)),m_qFibreCher_log,name,m_expHall_log,false,0,m_checkOverlaps));
      }
    }
  
    // place clad and buffer in the first fibre
    m_claddingCher = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),m_cladCher_log,"Qclad",m_qFibreCher_log,false,0,m_checkOverlaps);
    m_bufferCher = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),m_buffCher_log,"buffer",m_qFibreCher_log,false,0,m_checkOverlaps);

  }
  else if ( m_fbtype == FBphase1 ) { 
    // place phase-I fibers in a 15x15 grid
    const unsigned nSide = 15U;
    for ( unsigned i=0; i != nSide; i++ ) {
      const double yPos = 14.*m_rS1Fib - i*2.*m_rS1Fib;
      for ( unsigned j=0; j != nSide; j++ ) {
        const double zPos = 10.*cm + 14.*m_rS1Fib - j*2.*m_rS1Fib;
        char name[50];
        sprintf(name,"Sfib%d",i*nSide+j);
        m_fibresScin.push_back(new G4PVPlacement(G4Transform3D(rot,G4ThreeVector(0.,yPos,zPos)),m_qFibreScin1_log,name,m_expHall_log,false,0,m_checkOverlaps));
      }
    }
    // place clad and buffer in the first fibre
    m_jacketScin = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),m_jacketScin1_log,"Sjacket",m_qFibreScin1_log,false,0,m_checkOverlaps);
    m_claddingScin = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),m_cladScin1_log,"Sclad",m_qFibreScin1_log,false,0,m_checkOverlaps);
    m_bufferScin = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),m_buffScin1_log,"buffer",m_qFibreScin1_log,false,0,m_checkOverlaps); 
  }

  else if ( m_fbtype == FBphase2 ) {
    // place phase-II fibers in a 15x15 grid
    const unsigned nSide = 15U;
    for ( unsigned i=0; i != nSide; i++ ) {
      const double yPos = 14.*m_rS2Fib - i*2.*m_rS2Fib;
      for ( unsigned j=0; j != nSide; j++ ) {
        const double zPos = 10.*cm + 14.*m_rS2Fib - j*2.*m_rS2Fib;
        char name[50];
        sprintf(name,"Sfib%d",i*nSide+j);
        m_fibresScin.push_back(new G4PVPlacement(G4Transform3D(rot,G4ThreeVector(0.,yPos,zPos)),m_qFibreScin2_log,name,m_expHall_log,false,0,m_checkOverlaps));
      }
    }
    // place clad and buffer in the first fibre
    m_jacketScin = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),m_jacketScin2_log,"Sjacket",m_qFibreScin2_log,false,0,m_checkOverlaps);
    m_claddingScin = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),m_cladScin2_log,"Sclad",m_qFibreScin2_log,false,0,m_checkOverlaps);
    m_bufferScin = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),m_buffScin2_log,"buffer",m_qFibreScin2_log,false,0,m_checkOverlaps);
  }


  // place two scintillating fibers
  /*const double x = sqrt((m_rCFib+m_rSFib)*(m_rCFib+m_rSFib)-m_rCFib*m_rCFib);
  m_fibresScin.push_back(new G4PVPlacement(0,G4ThreeVector(x,0.,0.),m_qFibreScin_log,"Sfib0",m_expHall_log,false,0,m_checkOverlaps));
  m_fibresScin.push_back(new G4PVPlacement(0,G4ThreeVector(-x,0.,0.),m_qFibreScin_log,"Sfib1",m_expHall_log,false,1,m_checkOverlaps));*/
  // place clad and biffer and jacket on scintillating fiber
  /*m_jacketScin = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),m_jacketScin_log,"jacket",m_qFibreScin_log,false,0,m_checkOverlaps);
  m_claddingScin = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),m_cladScin_log,"Sclad",m_qFibreScin_log,false,0,m_checkOverlaps);
  m_bufferScin = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),m_buffScin_log,"buffer",m_qFibreScin_log,false,0,m_checkOverlaps);*/


  // place light guide (collar + pmt window)
  m_collar_phys = new G4PVPlacement(G4Transform3D(rot,G4ThreeVector(-m_length-0.9*m_collarLength/2.,0.,10.*cm)),m_collar_log,"collar",m_expHall_log,false,0,m_checkOverlaps);
  m_pmtWindow_phys = new G4PVPlacement(G4Transform3D(rot,G4ThreeVector(-m_length-1.*mm-m_pmtWindowThickness/2.,0.,10.*cm)),m_pmtWindow_log,"pmtWindow",m_expHall_log,false,0,m_checkOverlaps);
  m_pmtSpace_phys = new G4PVPlacement(G4Transform3D(rot,G4ThreeVector(-m_length-1.*mm-1.5*m_pmtWindowThickness,0.,10.*cm)),m_pmtSpace_log,"pmtSpace",m_expHall_log,false,0,m_checkOverlaps);


  G4OpticalSurface *photoCathode = new G4OpticalSurface("photoCathode");
  new G4LogicalBorderSurface("PhotoCathode",m_pmtWindow_phys,m_pmtSpace_phys,photoCathode);

  photoCathode->SetType(dielectric_metal);
  photoCathode->SetFinish(polished);
  photoCathode->SetModel(glisur);

  photoCathode->SetMaterialPropertiesTable(m_vacMPT);

  // skin surface
  //new G4LogicalSkinSurface("PhotoCathode",m_pmtSpace_log,photoCathode);

}


void CMSHFDetectorConstruction::CalculateConstants()
{
  
  if ( m_stacking ) {
     m_stacking->SetFiberRadius(m_rCClad);
     m_stacking->SetFibLength(m_length);
  }

  // calculate segment positions
  const double xOffset = 0.;
  const double yOffset = 0.*cm;
  const double zOffset = 20.*cm;

  m_segPositions.resize(m_Nseg);

 
  G4cout << " == Updating segment positions " << m_Nseg << " segments == " << G4endl; 
  for ( unsigned i=0; i != m_Nseg; i++ ) {
    m_segPositions[i] = G4ThreeVector( xOffset, yOffset, zOffset+i*m_absDepth);
    G4cout << "  [i] " << m_segPositions[i] << G4endl;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void CMSHFDetectorConstruction::SetPositionXYZ(const G4ThreeVector &detPos)
{
  m_xPos = detPos.x();
  m_yPos = detPos.y();
  m_zPos = detPos.z();

  /*if ( !m_isConstructed ) return;

  ClearPhysicalVolumes();
  //ClearLogicalVolumes();

  CalculateConstants();
  //SetupGeometry();
  SetupDetectors();

  G4RunManager::GetRunManager()->GeometryHasBeenModified();*/

}


void CMSHFDetectorConstruction::SetLength(G4double l)
{

  m_length = l;
  if ( m_stacking ) {
    m_stacking->SetFibLength(m_length);
  }

  /*if ( !m_isConstructed ) return;

  ClearPhysicalVolumes();
  ClearLogicalVolumes();

  CalculateConstants();
  SetupGeometry();
  SetupDetectors(); 

  //if ( m_gun ) {
  //  m_gun->SetInitDist(m_length);
  //}

  G4RunManager::GetRunManager()->GeometryHasBeenModified(); */

}


void CMSHFDetectorConstruction::SetNSeg(unsigned N)
{

  /*if ( m_isConstructed ) {
    ClearPhysicalVolumes();
    ClearLogicalVolumes();
  
    m_Nseg = N;

    CalculateConstants();
    SetupGeometry();
    SetupDetectors();
  
    G4RunManager::GetRunManager()->GeometryHasBeenModified();

  }*/

  m_Nseg = N;
  
}

void CMSHFDetectorConstruction::SetDeadTop(unsigned N)
{

  m_deadSeg_top = N;

  /*if ( m_isConstructed ) {
    ClearPhysicalVolumes();
    ClearLogicalVolumes();
  

    CalculateConstants();
    SetupGeometry();
    SetupDetectors();
  
    G4RunManager::GetRunManager()->GeometryHasBeenModified();

  }*/
  
}

void CMSHFDetectorConstruction::SetDeadBottom(unsigned N)
{
  m_deadSeg_bottom = N;

  /*if ( m_isConstructed ) {
    ClearPhysicalVolumes();
    ClearLogicalVolumes();
  

    CalculateConstants();
    SetupGeometry();
    SetupDetectors();
  
    G4RunManager::GetRunManager()->GeometryHasBeenModified();

  }*/
  
}

void CMSHFDetectorConstruction::SetDeadRight(unsigned N)
{
  m_deadSeg_right = N;

  /*if ( m_isConstructed ) {
    ClearPhysicalVolumes();
    ClearLogicalVolumes();
  

    CalculateConstants();
    SetupGeometry();
    SetupDetectors();
  
    G4RunManager::GetRunManager()->GeometryHasBeenModified();

  }*/

}

void CMSHFDetectorConstruction::SetDeadLeft(unsigned N)
{
  m_deadSeg_left = N;

  /*if ( m_isConstructed ) {
    ClearPhysicalVolumes();
    ClearLogicalVolumes();

    CalculateConstants();
    SetupGeometry();
    SetupDetectors();
  
    G4RunManager::GetRunManager()->GeometryHasBeenModified();

  }*/

}


void CMSHFDetectorConstruction::SetPitchAndYaw(G4double pitch, G4double yaw) {

  m_pitch = pitch;
  m_yaw = yaw;

  /*if ( m_isConstructed ){
    ClearPhysicalVolumes();
    ClearLogicalVolumes();
  
    CalculateConstants();
    SetupGeometry();
    SetupDetectors();
   
    G4RunManager::GetRunManager()->GeometryHasBeenModified(); 
    
  }*/
}


void CMSHFDetectorConstruction::SetFibres(const bool b) {
  m_fillFibres = b;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}



void CMSHFDetectorConstruction::SetScinYield(const G4double yield )
{
  m_scinYield = yield;
  
  if ( !m_isConstructed ) return;

  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void CMSHFDetectorConstruction::SetScinFastConst(const G4double tau)
{
  m_scinFastConst = tau;

  if ( !m_isConstructed ) return;

  G4RunManager::GetRunManager()->GeometryHasBeenModified();

}

void CMSHFDetectorConstruction::SetScinSlowConst(const G4double tau)
{
  m_scinSlowConst = tau;

  if ( !m_isConstructed ) return;

  G4RunManager::GetRunManager()->GeometryHasBeenModified();

}

void CMSHFDetectorConstruction::SetScinYieldRatio(const G4double r)
{
  m_scinYieldRatio = r;

  if ( !m_isConstructed ) return;

  G4RunManager::GetRunManager()->GeometryHasBeenModified();

} 

void CMSHFDetectorConstruction::SetFiberType( const FiberType fb)
{
  m_fbtype = fb;

  if ( !m_isConstructed ) return;
 
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}


void CMSHFDetectorConstruction::SetOverlapCheck(G4bool check)
{
  m_checkOverlaps = check;

  if ( !m_isConstructed ) return;

  G4RunManager::GetRunManager()->GeometryHasBeenModified(); 
  
}


void CMSHFDetectorConstruction::ClearPhysicalVolumes()
{

  assert(m_isConstructed);

  const unsigned segTot = m_Nseg;

  m_absBlock_log->RemoveDaughter(m_absBlock_phys);
  delete m_absBlock_phys;

  m_glass_log->RemoveDaughter( m_glass_phys );
  delete m_glass_phys;

  const unsigned nLead = m_pbBlock_phys.size();
  for ( unsigned i=0; i != nLead; i++ ) {
    m_pbBlock_log->RemoveDaughter( m_pbBlock_phys[i] );
    delete m_pbBlock_phys[i];
  }
  m_pbBlock_phys.clear();


  m_buffCher_log->ClearDaughters();
  m_cladCher_log->ClearDaughters();
  m_qFibreCher_log->ClearDaughters();
  m_absBlock_log->ClearDaughters();

  // clear daughters of logical volumes
  if ( m_jacketScin1_log ) {
    m_jacketScin1_log->ClearDaughters();
    m_buffScin1_log->ClearDaughters();
    m_cladScin1_log->ClearDaughters();  
    m_qFibreScin1_log->ClearDaughters();
  }


  if ( m_jacketScin2_log ) {
     m_jacketScin2_log->ClearDaughters();
    m_buffScin2_log->ClearDaughters(); 
    m_cladScin2_log->ClearDaughters();
    m_qFibreScin2_log->ClearDaughters();
  }

  const unsigned nFibs = m_fibresCher.size();
  for ( unsigned i=0; i != nFibs; i++ ) {
    delete m_fibresCher[i];
  }
  m_fibresCher.clear();

  const unsigned nS1Fibs = m_fibresScin.size();
  for ( unsigned i=0; i != nS1Fibs; i++ ) {
    delete m_fibresScin[i];
  }
  m_fibresScin.clear();

  if ( m_collar_phys ) { 
    m_collar_log->RemoveDaughter( m_collar_phys );
    delete m_collar_phys;
    m_collar_phys = NULL;
  }

  if ( m_pmtWindow_phys ) {
    m_pmtWindow_log->RemoveDaughter( m_pmtWindow_phys ) ;
    delete m_pmtWindow_phys;
    m_pmtWindow_phys = NULL;
  }

  if ( m_pmtSpace_phys ) {
    m_pmtSpace_log->RemoveDaughter( m_pmtSpace_phys );
    delete m_pmtSpace_phys;
    m_pmtSpace_phys = NULL;
  }


}

void CMSHFDetectorConstruction::ClearLogicalVolumes()
{

  assert(m_isConstructed);

  delete m_absBlock_log;
  delete m_qFibreCher_log;
  delete m_cladCher_log;
  delete m_buffCher_log;

  if ( m_qFibreScin1_log ) {
    delete m_qFibreScin1_log;
    delete m_cladScin1_log;
    delete m_buffScin1_log;
    delete m_jacketScin1_log;
    m_qFibreScin1_log = 0;
    m_cladScin1_log = 0;
    m_buffScin1_log = 0;
    m_jacketScin1_log = 0;
  }

  if ( m_qFibreScin2_log ) {
    delete m_qFibreScin2_log;
    delete m_cladScin2_log;
    delete m_buffScin2_log;
    delete m_jacketScin2_log;
    m_qFibreScin2_log = 0;
    m_cladScin2_log = 0;
    m_buffScin2_log = 0;
    m_jacketScin2_log = 0;
  }

  delete m_glass_log;

  delete m_collar_log;
  delete m_pmtWindow_log;
  delete m_pmtSpace_log;

}

void CMSHFDetectorConstruction::RefreshGeometry()
{
    ClearPhysicalVolumes();
    ClearLogicalVolumes();

    DefineScintillator();
    CalculateConstants();
    SetupGeometry();
    SetupDetectors();
  
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void CMSHFDetectorConstruction::SetMagneticField(const G4ThreeVector &vec)
{

  if ( m_Bfield ) delete m_Bfield;

  m_Bfield = new G4UniformMagField(vec);

  G4FieldManager *fldMngr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fldMngr->SetDetectorField(m_Bfield);
  fldMngr->CreateChordFinder(m_Bfield);

  G4cout << "Creating Magnetic Field: " << vec << G4endl;

}



void CMSHFDetectorConstruction::SetStackingAction( HFStackingAction *sa )
{
  m_stacking = sa;
}



