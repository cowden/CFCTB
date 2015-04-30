#ifndef HFDataFormat_hh
#define HFDataFormat_hh

///////////////////////////////////////////////////
// HFDataFormat
// Root ntuple dumping class and helper functions
// to coordinate the output data from the various 
// geant4 user action classes.
//////////////////////////////////////////////////

#include <vector>
#include <string>
#include <iostream>


#include "TFile.h"
#include "TTree.h"

#include "globals.hh"

#include "G4ThreeVector.hh"
#include "HFDataFormatMessenger.hh"

///
///
/// enumeration of readout type (Cherenkov or scintillation)
enum ROType {
 fCherenkov=0,
 fScintillation,
 fFiber
};

///
///
/// enumeration of detector faces (sides, front, back)
enum Face {
 fSide=0,
 fFront,
 fBack
};

///
/// The StackingStruct helps to pass data from the HFStackingAction
/// to the ntuple filling procedure.
struct StackingStruct {
  double wavelength;
  double energy;
  double na;
  double x;
  double y;
  double depth;
  double t;
  double tprop;


  inline StackingStruct(double w, double e, double n,double gx, double gy,double gd, double gt, double pt)
    :wavelength(w),energy(e),na(n)
    ,x(gx),y(gy),depth(gd),t(gt),tprop(pt)
    { }
};


///
/// The ParticleStruct helps to pass data about particles
/// As tracks are created, the HFStackingAction class can record
/// these in an ntuple for studies of the shower development in
/// simulation.
struct ParticleStruct {
  int pdgID;
  double px;
  double py;
  double pz;
  double x;
  double y;
  double z;
  double e;

  inline ParticleStruct(int ID, const G4ThreeVector & mom, const G4ThreeVector & pos,double E)
    :pdgID(ID),px(mom.x()),py(mom.y()),pz(mom.z())
    ,x(pos.x()),y(pos.y()),z(pos.z())
    ,e(E)
    { }

};

///
/// The GeneratorStruct helps to pass data about the intial particle
/// to the ntuple
struct GeneratorStruct {
  double x;
  double y;
  double e;

  inline GeneratorStruct(double gx, double gy, double ge)
    :x(gx),y(gy),e(ge)
    { }
  
};

///
/// The SteppingStruct helps to pass data from the stepping action
/// mostly information for PMT simulation
struct SteppingStruct {
  double x;
  double y;
  double z;
  double t;
  double lt; // local time (i.e. time track has been alive)
  double tl; // track length
  double lambda;
  double polX;
  double polY;
  double vx;
  double vy;
  double vz;

  inline SteppingStruct(const G4ThreeVector &mom, double pt,
    double localTime, double tLength, double l, double ppx, double ppy, 
    const G4ThreeVector &vert)
    :x(mom.x()),y(mom.y()),z(mom.z()),t(pt),lt(localTime),tl(tLength),lambda(l),polX(ppx),polY(ppy)
    ,vx(vert.x()),vy(vert.y()),vz(vert.z())
    { }
};

///
/// The ionization struct helps pass energy loss in a step from stepping action to
/// the ntuple
struct IoniStruct {
  double E;
  double x;
  double y;
  double depth;
  double t;

  inline IoniStruct(const double e, const G4ThreeVector &pos, double d, double T)
    :E(e),x(pos.x()),y(pos.y()),depth(d),t(T)
    { }

};


class HFDataFormat {

public:

  // constructor, provide output ROOT file name
  HFDataFormat(const std::string &fileName);

  virtual ~HFDataFormat();


  // fill from StackingAction
  void fillStackingAction(const StackingStruct &);
  void fillStackingAction(const StackingStruct &, const ROType);
  // fill from EventAction 

  // fill shower particles when new track is created
  void fillParticle(const ParticleStruct &);
 
  // accummulate leakage
  void accLeakage(Face f, const double E);

  // fill from primary generator action
  void fillGenerator(const GeneratorStruct &);

  // fill from stepping action (PMT);
  void fillSteppingAction(const SteppingStruct &);
  void fillSteppingAction(const SteppingStruct &, const ROType);
  void fillIonizationCore(const IoniStruct &);
  void fillIonizationQuartz(const IoniStruct &);
  void fillIonizationClad(const IoniStruct &);

  

  // store event in tree, and clear vectors
  void store();

  // dump file and close   
  void fileDump();

  // set the file name and create tree
  void SetFileName(const G4String &fileName);

  // generate tree branches
  void generateTrees();

  // set and get store options
  void SetStoreOpticalInfo(G4bool store)    { _storeOpticalInfo = store; }
  void SetStoreParticleInfo(G4bool store)    { _storeParticleInfo = store; }
  void SetStoreGeneratorInfo(G4bool store)    { _storeGeneratorInfo = store; }
  void SetStorePMTInfo(G4bool store) { _storePMTInfo = store; }
 
  G4bool GetStoreOpticalInfo()          { return _storeOpticalInfo; }
  G4bool GetStoreParticleInfo()          { return _storeParticleInfo; }
  G4bool GetStoreGeneratorInfo()          { return _storeGeneratorInfo; }
  G4bool GetStorePMTInfo()  { return _storePMTInfo; }
 

private:


  // messenger
  HFDataFormatMessenger* _messenger;
  // options in storage
  G4bool _storeOpticalInfo;
  G4bool _storeParticleInfo;
  G4bool _storeGeneratorInfo;
  G4bool _storePMTInfo;

  // clear vectors in trees
  // clear stacking vectors
  void clearStacking();
  // clear pmt vectors
  void clearPMT();
  // clear particle vectors
  void clearParticle();
  // clear generator vectors
  void clearGenerator();
  // clear leakage variables
  void clearLeakage();

  // ------------------- member data -------------------
  TFile * m_file;

  TTree * m_event;

  // tree branches
  // event branches
  std::vector<double>  m_opt_wavelength;
  std::vector<double>  m_opt_energy;
  std::vector<double>  m_opt_na;
  std::vector<double>  m_opt_fx;
  std::vector<double>  m_opt_fy;
  std::vector<double>  m_opt_depth;
  std::vector<double>  m_opt_t;
  std::vector<double>  m_opt_tprop;

  std::vector<double>  m_scin_wavelength;
  std::vector<double>  m_scin_energy;
  std::vector<double>  m_scin_na;
  std::vector<double>  m_scin_fx;
  std::vector<double>  m_scin_fy;
  std::vector<double>  m_scin_depth;
  std::vector<double>  m_scin_t;
  std::vector<double>  m_scin_tprop;

  // energy loss in fibres
  std::vector<double> m_scinIon_E;
  std::vector<double> m_scinIon_t;
  std::vector<double> m_scinIon_x;
  std::vector<double> m_scinIon_y;
  std::vector<double> m_scinIon_depth;

  std::vector<double> m_quartzIon_E;
  std::vector<double> m_quartzIon_t;
  std::vector<double> m_quartzIon_x;
  std::vector<double> m_quartzIon_y;
  std::vector<double> m_quartzIon_depth;

  std::vector<double> m_cladIon_E;
  std::vector<double> m_cladIon_t;
  std::vector<double> m_cladIon_x;
  std::vector<double> m_cladIon_y;
  std::vector<double> m_cladIon_depth;


  // photons absorbed by photocathode 
  std::vector<double> m_pmt_x;
  std::vector<double> m_pmt_y;
  std::vector<double> m_pmt_z;
  std::vector<double> m_pmt_vx;
  std::vector<double> m_pmt_vy;
  std::vector<double> m_pmt_vz;
  std::vector<double> m_pmt_t;
  std::vector<double> m_pmt_lt;
  std::vector<double> m_pmt_tl;
  std::vector<double> m_pmt_wavelength;
  std::vector<double> m_pmt_polX;
  std::vector<double> m_pmt_polY;

  // photons absorbed by the photocathode
  std::vector<double> m_pmtScin_x;
  std::vector<double> m_pmtScin_y;
  std::vector<double> m_pmtScin_z;
  std::vector<double> m_pmtScin_vx;
  std::vector<double> m_pmtScin_vy;
  std::vector<double> m_pmtScin_vz;
  std::vector<double> m_pmtScin_t;
  std::vector<double> m_pmtScin_lt;
  std::vector<double> m_pmtScin_tl;
  std::vector<double> m_pmtScin_wavelength;
  std::vector<double> m_pmtScin_polX;
  std::vector<double> m_pmtScin_polY;


  std::vector<double> m_fib_x;
  std::vector<double> m_fib_y;
  std::vector<double> m_fib_z;
  std::vector<double> m_fib_vx;
  std::vector<double> m_fib_vy;
  std::vector<double> m_fib_vz;
  std::vector<double> m_fib_t;
  std::vector<double> m_fib_lt;
  std::vector<double> m_fib_tl;
  std::vector<double> m_fib_wavelength;
  std::vector<double> m_fib_polX;
  std::vector<double> m_fib_polY;

  // shower particle branches
  std::vector<int>  m_part_pdgId;
  std::vector<double>  m_part_px;
  std::vector<double>  m_part_py;
  std::vector<double>  m_part_pz;
  std::vector<double>  m_part_x;
  std::vector<double>  m_part_y;
  std::vector<double>  m_part_z;
  std::vector<double>  m_part_e;

  // generator (beam) parameters
  std::vector<double>  m_gen_x;
  std::vector<double>  m_gen_y;
  std::vector<double>  m_gen_e;

  // leakage
  double m_totLeak;
  double m_latLeak;
  double m_frontLeak;
  double m_backLeak;


};


#endif
