
#include "HFSteppingAction.hh"

#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

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
  
  // get next volume
  G4VPhysicalVolume * postVolume = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume();

  const G4String & preName = preVolume->GetName();
  const bool isFibre = preName.contains("Cfib") || preName.contains("Qclad");
  const bool isScinFibre = preName.contains("Sfib") || preName.contains("Sclad") || preName.contains("jacket");

  // get the particle
  G4Track * theTrack = step->GetTrack();
  const G4ThreeVector & pos = theTrack->GetPosition();

  const double time = theTrack->GetGlobalTime();
  if ( time > 100*ns ) theTrack->SetTrackStatus(fStopAndKill);

  // record photons tracked to PMT face
  if ( theTrack->GetDefinition() == m_optDef && preVolume != postVolume 
	&&  ( isFibre || isScinFibre ) ) {

    const G4String & postName = postVolume->GetName();
    if ( postName.contains("glass") ) {

      const G4DynamicParticle * theParticle = theTrack->GetDynamicParticle();
      const double wavelength = hbarc*twopi/theParticle->GetTotalEnergy()*1.e+6;

      const G4ThreeVector & pol = theParticle->GetPolarization();
      const G4ThreeVector & vertpos = theTrack->GetVertexPosition();

      SteppingStruct st(pos,time,theTrack->GetLocalTime(),theTrack->GetTrackLength(),wavelength,pol.x(),pol.y(),vertpos);
      if ( wavelength > 350. )  {
  	if ( isFibre ) m_df->fillSteppingAction( st, fCherenkov );
  	else if ( isScinFibre ) m_df->fillSteppingAction( st, fScintillation );
      }

      // kill the track after readout
      theTrack->SetTrackStatus(fStopAndKill);

    } else if ( postName.contains("World") || postName.contains("absorber") ) {
      // kill the track since it left the fiber
      theTrack->SetTrackStatus(fStopAndKill);
    }

  } else if ( preVolume == postVolume && preName.contains("Sfib") ) {
    const double E = step->GetTotalEnergyDeposit();

    const G4ThreeVector & touchTrans = theTrack->GetTouchable()->GetTranslation();
    const double dist = (touchTrans.z()+m_fibLength/2.-pos.z())/m_fibLength;
    const double depth = m_fibLength*(1.-dist);

    IoniStruct is(E,pos,depth,time);
    m_df->fillIonization(is); 
  } 

}

