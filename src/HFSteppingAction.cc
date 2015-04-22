
#include "HFSteppingAction.hh"

#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4SteppingManager.hh"
#include "G4ProcessManager.hh"

#include "HFDataFormat.hh"

HFSteppingAction::HFSteppingAction(HFDataFormat *df)
:m_df(df)
,m_optDef( G4OpticalPhoton::OpticalPhotonDefinition() )
,m_fibLength(2.*m)
{ }

HFSteppingAction::~HFSteppingAction()
{ }

void HFSteppingAction::UserSteppingAction(const G4Step * step)
{

  // get volume of the current step
  G4VPhysicalVolume* preVolume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  const G4String & preName = preVolume->GetName();
  
  // get next volume
  G4VPhysicalVolume * postVolume = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume();


  // get the particle
  G4Track * theTrack = step->GetTrack();
  const G4ThreeVector & pos = theTrack->GetPosition();

  const double time = theTrack->GetGlobalTime();
  if ( time > 100*ns ) theTrack->SetTrackStatus(fStopAndKill);

  const G4String & postName = postVolume->GetName();


  // store ionization energy loss in the scintillating fiber
  if ( theTrack->GetDefinition() != m_optDef && preVolume == postVolume && preName.contains("Sfib") ) {
    // collect ionization energy loss in the scintillating fiber core
    const double E = step->GetTotalEnergyDeposit();

    const G4ThreeVector & touchTrans = theTrack->GetTouchable()->GetTranslation();
    const double dist = (touchTrans.z()+m_fibLength/2.-pos.z())/m_fibLength;
    const double depth = m_fibLength*(1.-dist);

    IoniStruct is(E,pos,depth,time);
    m_df->fillIonization(is); 
  } 

  // record photons collected at the photo-cathod
  // by analyzing the boundary process
  const  G4StepPoint * thePostPoint = step->GetPostStepPoint();
  const G4VPhysicalVolume *thePostPV = thePostPoint->GetPhysicalVolume();
  if ( theTrack->GetDefinition() == m_optDef 
    && thePostPoint->GetStepStatus() == fGeomBoundary ) {
    
    // get the boundary process
    static G4ThreadLocal G4OpBoundaryProcess *boundary=NULL;
    if ( !boundary ) {
      G4ProcessManager *pm = step->GetTrack()->GetDefinition()->GetProcessManager();
      G4int nprocesses = pm->GetProcessListLength();
      G4ProcessVector *pv = pm->GetProcessList();
      for ( G4int i=0; i != nprocesses; i++ ){
	if ((*pv)[i]->GetProcessName() == "OpBoundary") {
	  boundary = (G4OpBoundaryProcess*)(*pv)[i];
	  break;
	}
      }
    }

    const G4OpBoundaryProcessStatus boundaryStatus = boundary->GetStatus(); 

    //G4cout << "Boundary: " << preName << " => " << postName << " status: " << boundaryStatus << G4endl;

    // check for detection.
    // since the photo-cathode is the only material which is given an efficiency
    // it is the only physcal volume considered
    if ( boundaryStatus == Detection ) {
      const G4String & origName = step->GetTrack()->GetOriginTouchable()->GetVolume()->GetName();
      const bool isFibre = origName.contains("Cfib") || origName.contains("Qclad");
      const bool isScinFibre = origName.contains("Sfib") || origName.contains("Sclad") || origName.contains("jacket");
      G4cout << "Found boundary detection!! " << origName << G4endl;
      
      const G4DynamicParticle * theParticle = theTrack->GetDynamicParticle();
      const double wavelength = hbarc*twopi/theParticle->GetTotalEnergy()*1.e+6;

      const G4ThreeVector & pol = theParticle->GetPolarization();
      const G4ThreeVector & vertpos = theTrack->GetVertexPosition();

      SteppingStruct st(pos,time,theTrack->GetLocalTime(),theTrack->GetTrackLength(),wavelength,pol.x(),pol.y(),vertpos);
     
      if ( isFibre ) m_df->fillSteppingAction( st, fCherenkov );
      else if ( isScinFibre ) m_df->fillSteppingAction( st, fScintillation );
    }


  } 



}

