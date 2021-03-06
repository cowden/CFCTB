
#include "CMSHFDetectorConstructionMessenger.hh"

#include "CMSHFDetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"

#include <iostream>


CMSHFDetectorConstructionMessenger::CMSHFDetectorConstructionMessenger(CMSHFDetectorConstruction *det)
:m_detector(det)
{

  std::cout << "CMSHF Constructing messenger" << std::endl;
  m_dir = new G4UIdirectory("/testBeam/");
  m_dir->SetGuidance("CMSHF W calorimeter control");

  m_posCmd = new G4UIcmdWith3VectorAndUnit("/testBeam/detPosition", this);
  m_posCmd->SetGuidance("Set the position of the center of the detector element.");
  m_posCmd->SetParameterName("X","Y","Z",true);
  m_posCmd->SetDefaultValue(G4ThreeVector(0.,0.,0.7));
  m_posCmd->SetDefaultUnit("m");  
  m_posCmd->AvailableForStates(G4State_Idle);

  m_lengthCmd = new G4UIcmdWithADoubleAndUnit("/testBeam/length",this);
  m_lengthCmd->SetGuidance("Set the 1/2 Length of the calorimeter");
  m_lengthCmd->SetParameterName("length",false);
  m_lengthCmd->SetDefaultValue(60);
  m_lengthCmd->SetDefaultUnit("cm");
  m_lengthCmd->AvailableForStates(G4State_Idle);

  m_NfibCmd = new G4UIcmdWithAnInteger("/testBeam/Nfib",this);
  m_NfibCmd->SetGuidance("Set the number of fibers along side of segment");
  m_NfibCmd->SetParameterName("Nfib",false);
  m_NfibCmd->SetDefaultValue(32);
  m_NfibCmd->AvailableForStates(G4State_Idle);

  m_NsegCmd = new G4UIcmdWithAnInteger("/testBeam/Nseg",this);
  m_NsegCmd->SetGuidance("Set the number of segments along side of detector");
  m_NsegCmd->SetParameterName("Nseg",false);
  m_NsegCmd->SetDefaultValue(13);
  m_NsegCmd->AvailableForStates(G4State_Idle);

  m_fibreIndexCmd =  new G4UIcmdWithADouble("/testBeam/fibreIndex",this);
  m_fibreIndexCmd->SetGuidance("Set the fibre index of refraction");
  m_fibreIndexCmd->SetParameterName("fibreIndex",false);
  m_fibreIndexCmd->SetDefaultValue(1.457);
  m_fibreIndexCmd->AvailableForStates(G4State_Idle);
  
  m_cladIndexCmd =  new G4UIcmdWithADouble("/testBeam/cladIndex",this);
  m_cladIndexCmd->SetGuidance("Set the clad index of refraction");
  m_cladIndexCmd->SetParameterName("cladIndex",false);
  m_cladIndexCmd->SetDefaultValue(1.42);
  m_cladIndexCmd->AvailableForStates(G4State_Idle);

  m_overlapCheckCmd = new G4UIcmdWithABool("/testBeam/checkOverlaps",this);
  m_overlapCheckCmd->SetGuidance("Set flag to check geometry for overlapping volumes");
  m_overlapCheckCmd->SetDefaultValue(false);
  m_overlapCheckCmd->AvailableForStates(G4State_Idle);

  m_fieldCmd = new G4UIcmdWith3VectorAndUnit("/testBeam/magField", this);
  m_fieldCmd->SetGuidance("Define a uniform magnetic field");
  m_fieldCmd->SetParameterName("X","Y","Z",true);
  m_fieldCmd->SetDefaultValue(G4ThreeVector(0.,0.,0.7));
  m_fieldCmd->SetDefaultUnit("tesla");  
  m_fieldCmd->AvailableForStates(G4State_Idle);

  m_rotCmd = new G4UIcmdWith3VectorAndUnit("/testBeam/rotate", this);
  m_rotCmd->SetGuidance("Define the pitch and yaw of the detector.  Requires a 3 vector (the third element is a dummy inpur)");
  m_rotCmd->SetDefaultValue(G4ThreeVector(0.,0.,0.));
  m_rotCmd->SetDefaultUnit("deg");
  m_rotCmd->AvailableForStates(G4State_Idle);

  m_fillFibres = new G4UIcmdWithABool("/testBeam/fibres",this);
  m_fillFibres->SetGuidance("Turn on fibres, this is important");
  m_fillFibres->SetDefaultValue(false);
  m_fillFibres->AvailableForStates(G4State_Idle);

  m_NsegTopCmd = new G4UIcmdWithAnInteger("/testBeam/NdeadTop",this);
  m_NsegTopCmd->SetGuidance("Set number of dead segments on top of block");
  m_NsegTopCmd->SetDefaultValue(1);
  m_NsegTopCmd->AvailableForStates(G4State_Idle);

  m_NsegBottomCmd = new G4UIcmdWithAnInteger("/testBeam/NdeadBottom",this);
  m_NsegBottomCmd->SetGuidance("Set number of dead segments on bottom of block");
  m_NsegBottomCmd->SetDefaultValue(1);
  m_NsegBottomCmd->AvailableForStates(G4State_Idle);

  m_NsegLeftCmd = new G4UIcmdWithAnInteger("/testBeam/NdeadLeft",this);
  m_NsegLeftCmd->SetGuidance("Set number of dead segments on left of block");
  m_NsegLeftCmd->SetDefaultValue(1);
  m_NsegLeftCmd->AvailableForStates(G4State_Idle);

  m_NsegRightCmd = new G4UIcmdWithAnInteger("/testBeam/NdeadRight",this);
  m_NsegRightCmd->SetGuidance("Set number of dead segments on right of block");
  m_NsegRightCmd->SetDefaultValue(1);
  m_NsegRightCmd->AvailableForStates(G4State_Idle);

  m_refreshCmd = new G4UIcmdWithABool("/testBeam/refreshGeometry",this);
  m_refreshCmd->SetGuidance("Force refresh of geometry and rebuild");
  m_refreshCmd->SetDefaultValue(true);
  m_refreshCmd->AvailableForStates(G4State_Idle);

  m_materialCmd = new G4UIcmdWithABool("/testBeam/buildWorBrass",this);
  m_materialCmd->SetGuidance("Build detector with W (true) or Brass (false) absorber");
  m_materialCmd->SetDefaultValue(true);
  m_materialCmd->AvailableForStates(G4State_Idle);

  m_fibSpaceCmd = new G4UIcmdWithADoubleAndUnit("/testBeam/fibreSpacing",this);
  m_fibSpaceCmd->SetGuidance("Set the center to center separation of fibres");
  m_fibSpaceCmd->SetDefaultValue(1.2);
  m_fibSpaceCmd->SetDefaultUnit("mm");
  m_fibSpaceCmd->AvailableForStates(G4State_Idle);

  m_scinYieldCmd = new G4UIcmdWithADouble("/testBeam/scinYield",this);
  m_scinYieldCmd->SetGuidance("Set the scintillation yield in photons/MeV");
  m_scinYieldCmd->SetDefaultValue(8300);
  m_scinYieldCmd->AvailableForStates(G4State_Idle);

  m_scinYieldRatioCmd = new G4UIcmdWithADouble("/testBeam/scinYieldRatio",this);
  m_scinYieldRatioCmd->SetGuidance("Set the ratio of scintillation light between fast and slow components");
  m_scinYieldRatioCmd->SetDefaultValue(0.8);
  m_scinYieldRatioCmd->AvailableForStates(G4State_Idle);

  m_scinFastConstCmd = new G4UIcmdWithADoubleAndUnit("/testBeam/scinFastTimeConst",this);
  m_scinFastConstCmd->SetGuidance("Set the fast time constant for the scintillation response. ");
  m_scinFastConstCmd->SetDefaultValue(1.);
  m_scinFastConstCmd->SetDefaultUnit("ns");
  m_scinFastConstCmd->AvailableForStates(G4State_Idle);

  m_scinSlowConstCmd = new G4UIcmdWithADoubleAndUnit("/testBeam/scinSlowTimeConst",this);
  m_scinSlowConstCmd->SetGuidance("Set the slow time constant for the scintillation response.");
  m_scinSlowConstCmd->SetDefaultValue(10.);
  m_scinSlowConstCmd->SetDefaultUnit("ns");
  m_scinSlowConstCmd->AvailableForStates(G4State_Idle);

}

CMSHFDetectorConstructionMessenger::~CMSHFDetectorConstructionMessenger()
{

  delete m_posCmd;
  delete m_lengthCmd;
  delete m_NfibCmd;
  delete m_NsegCmd;
  delete m_fibreIndexCmd;
  delete m_cladIndexCmd;
  delete m_overlapCheckCmd;
  delete m_rotCmd;
  delete m_fillFibres;
  delete m_NsegTopCmd;
  delete m_NsegBottomCmd;
  delete m_NsegRightCmd;
  delete m_NsegLeftCmd;
  delete m_refreshCmd;
  delete m_materialCmd;
  delete m_fibSpaceCmd;
  delete m_scinYieldCmd;
  delete m_scinYieldRatioCmd;
  delete m_scinFastConstCmd;
  delete m_scinSlowConstCmd;

}

void CMSHFDetectorConstructionMessenger::SetNewValue(G4UIcommand* cmd,G4String newValue)
{

  std::cout << "CMSHF parsing commands " << cmd << " " << newValue << std::endl;

  if ( cmd == m_lengthCmd ) {
    m_detector->SetLength(m_lengthCmd->GetNewDoubleValue(newValue));
  } 
  else if ( cmd == m_posCmd ) {
    m_detector->SetPositionXYZ(m_posCmd->GetNew3VectorValue(newValue));
  }
  else if ( cmd == m_NfibCmd ) {
  }
  else if ( cmd == m_NsegCmd ) {
    m_detector->SetNSeg(m_NsegCmd->GetNewIntValue(newValue));
  }
  else if ( cmd == m_fibreIndexCmd ) {
  } 
  else if ( cmd == m_cladIndexCmd ) {
  }
  else if ( cmd == m_overlapCheckCmd ) {
    m_detector->SetOverlapCheck(m_overlapCheckCmd->GetNewBoolValue(newValue));
  }
  else if ( cmd == m_fieldCmd ) {
    m_detector->SetMagneticField(m_fieldCmd->GetNew3VectorValue(newValue));
  }
  else if ( cmd == m_rotCmd ) {
    const G4ThreeVector & vec = m_rotCmd->GetNew3VectorValue(newValue);
    m_detector->SetPitchAndYaw(vec.x(),vec.y());
  }
  else if ( cmd == m_fillFibres ) {
    m_detector->SetFibres( m_fillFibres->GetNewBoolValue(newValue));
  }
  else if ( cmd == m_NsegTopCmd ) {
    m_detector->SetDeadTop( m_NsegTopCmd->GetNewIntValue(newValue));
  }
  else if ( cmd == m_NsegBottomCmd ) {
    m_detector->SetDeadBottom( m_NsegBottomCmd->GetNewIntValue(newValue));
  }
  else if ( cmd == m_NsegLeftCmd ) {
    m_detector->SetDeadLeft( m_NsegLeftCmd->GetNewIntValue(newValue));
  }
  else if ( cmd == m_NsegRightCmd ) {
    m_detector->SetDeadRight( m_NsegRightCmd->GetNewIntValue(newValue));
  }
  else if ( cmd == m_refreshCmd ) {
    m_detector->RefreshGeometry();
  }
  else if ( cmd == m_materialCmd ) {
  }
  else if ( cmd == m_fibSpaceCmd ) {
  }
  else if ( cmd == m_scinYieldCmd ) {
    m_detector->SetScinYield( m_scinYieldCmd->GetNewDoubleValue(newValue)/MeV );
  }
  else if ( cmd == m_scinYieldRatioCmd ) {
    m_detector->SetScinYieldRatio( m_scinYieldRatioCmd->GetNewDoubleValue(newValue) );
  }
  else if ( cmd == m_scinFastConstCmd ) {
    m_detector->SetScinFastConst( m_scinFastConstCmd->GetNewDoubleValue(newValue) );
  }
  else if ( cmd == m_scinSlowConstCmd ) {
    m_detector->SetScinSlowConst( m_scinSlowConstCmd->GetNewDoubleValue(newValue) );
  }

}
