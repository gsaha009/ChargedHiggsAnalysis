#include "configana.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <functional>
#include <numeric>
#include <string>
#include <climits>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <sstream>

#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TH1K.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"

#include "AnaBase.h"

using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::string;
using std::vector;
using std::map;
using std::pair;
using std::setprecision;
using std::setw;
using std::setiosflags;
using std::resetiosflags;

// -----------
// Constructor
// -----------
AnaBase::AnaBase()
  : chain_(new TChain("Events")),
    chainRun_(new TChain("Runs"))
{
  cout << setiosflags(ios::fixed); 
  fileList_.clear();
  singleMuonHltPathList_.clear();
  singleElectronHltPathList_.clear();
  doubleMuonHltPathList_.clear();
  doubleEgHltPathList_.clear();
  muonEgHltPathList_.clear();
  singleMuonHltPtrList_.clear();
  singleElectronHltPtrList_.clear();
  doubleMuonHltPtrList_.clear();
  doubleEgHltPtrList_.clear();
  muonEgHltPtrList_.clear();
  brList_.clear();
}
// ----------
// Destructor
// ----------
AnaBase::~AnaBase() 
{
  clearEvent();
}
// ------------------------
// Clear the clones arrays
// ------------------------
void AnaBase::clearEvent() {
  //previously HLT paths, prescales and results were cleaned here
}

// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool AnaBase::beginJob() 
{ 
  //  if (isMC_ && usePUWt_ && !readPileUpHist()) return false;
  if (isMC() && !openScaleFactorRootFiles()) {
    cerr << "Scale Factors are not accessible !!!\n";
    return false;
  }
  
  // Open the output ROOT file
  TFile* f = TFile::Open(histFile_.c_str(), "RECREATE");
  histf_.reset(std::move(f));
  // Open the output FakeExtrapolated ROOT file
  //TFile* fakef = TFile::Open(fakehistFile_.c_str(), "RECREATE");
  //fakehistf_.reset(std::move(fakef));

  nEvents_ = static_cast<int>(chain_->GetEntries()); 
  if (nEvents_ <= 0) {
    cerr << "******* nEvents = " << nEvents_ << ", returning!" << endl;
    return false;
  }
  if (maxEvt_ > 0) nEvents_ = std::min(nEvents_, maxEvt_);
  cout << " >>> nEvents = " << nEvents_ << endl;

  std::cout<<"Chain Name: "<<chain()->GetName()<<"\t"<<"Chain Entries: "<<chain()->GetEntries()<<std::endl;
  treeReader_ = new TTreeReader(chain_);
  std::cout<<"treeReader with TChain::Events\n";
  std::cout<<"Chain Name: "<<chainRun()->GetName()<<"\t"<<"Chain Entries: "<<chainRun()->GetEntries()<<std::endl;
  treeReaderRun_ = new TTreeReader(chainRun_);
  std::cout<<"treeReader with TChain::Runs\n";

  if(!treeReader_) {
    std::cout << "Event Tree no found!!!\n";
    std::exit(EXIT_FAILURE);
  }
  else if(!treeReaderRun_) {
    std::cout << "Run Tree no found!!!\n";
    std::exit(EXIT_FAILURE);
  }

  treeReader_->ls();
  treeReaderRun_->ls();

  openFiles();

  return true;
}

// ---------------------------
// Initialize all the branches
// ---------------------------
bool AnaBase::init(){
  if (isMC_){
    //from TTree::Runs
    if (branchFound(evtWtSum_.c_str())) genEventSumw_                  = new TTreeReaderValue< double >(*treeReaderRun_, evtWtSum_.c_str());
    //from TTree::Events 
    if (branchFound("genWeight")) genEvWt                              = new TTreeReaderValue< float >(*treeReader_, "genWeight");
    if (branchFound("puWeight")) PU_Weight                             = new TTreeReaderValue< float >(*treeReader_, "puWeight");
    if (branchFound("puWeightUp")) PU_WeightUp                         = new TTreeReaderValue< float >(*treeReader_, "puWeightUp");
    if (branchFound("puWeightDown")) PU_WeightDown                     = new TTreeReaderValue< float >(*treeReader_, "puWeightDown");
    
    if (branchFound("btagWeight_CSVV2")) btagWeight_CSVV2              = new TTreeReaderValue< float >(*treeReader_, "btagWeight_CSVV2");
  }
  //TTree::Events
  if (branchFound("run")) run_                                         = new TTreeReaderValue< unsigned int >(*treeReader_, "run");
  if (branchFound("event")) event_                                     = new TTreeReaderValue< unsigned long long >(*treeReader_, "event");
  if (branchFound("luminosityBlock")) lumis_                           = new TTreeReaderValue< unsigned int >(*treeReader_, "luminosityBlock");

  //Primary Vertices
  if (branchFound("PV_ndof")) PV_ndf                                   = new TTreeReaderValue< float >(*treeReader_, "PV_ndof");
  if (branchFound("PV_x")) PV_xPos                                     = new TTreeReaderValue< float >(*treeReader_, "PV_x");
  if (branchFound("PV_y")) PV_yPos                                     = new TTreeReaderValue< float >(*treeReader_, "PV_y");
  if (branchFound("PV_z")) PV_zPos                                     = new TTreeReaderValue< float >(*treeReader_, "PV_z");
  if (branchFound("PV_chi2")) PV_chi2                                  = new TTreeReaderValue< float >(*treeReader_, "PV_chi2");
  if (branchFound("PV_score")) PV_score                                = new TTreeReaderValue< float >(*treeReader_, "PV_score");
  if (branchFound("PV_npvs")) nPV                                      = new TTreeReaderValue< int >(*treeReader_, "PV_npvs");
  if (branchFound("PV_npvsGood"))nGoodPV                               = new TTreeReaderValue< int >(*treeReader_, "PV_npvsGood");

  //Muon
  if (branchFound("nMuon")) nMuon                                      = new TTreeReaderValue< unsigned int >(*treeReader_, "nMuon");
  if (branchFound("Muon_pt")) Muon_pt                                  = new TTreeReaderArray< float >(*treeReader_, "Muon_pt");
  if (branchFound("Muon_corrected_pt")) Muon_corrpt                    = new TTreeReaderArray< float >(*treeReader_, "Muon_corrected_pt");
  if (branchFound("Muon_eta")) Muon_eta                                = new TTreeReaderArray< float >(*treeReader_, "Muon_eta");
  if (branchFound("Muon_phi")) Muon_phi                                = new TTreeReaderArray< float >(*treeReader_, "Muon_phi");
  if (branchFound("Muon_charge")) Muon_charge                          = new TTreeReaderArray< int >(*treeReader_, "Muon_charge");
  if (branchFound("Muon_mass")) Muon_mass                              = new TTreeReaderArray< float >(*treeReader_, "Muon_mass");
  if (branchFound("Muon_sip3d")) Muon_sip3d                            = new TTreeReaderArray< float >(*treeReader_, "Muon_sip3d");
  if (branchFound("Muon_jetIdx")) Muon_jetIdx                          = new TTreeReaderArray< int >(*treeReader_, "Muon_jetIdx");
  if (branchFound("Muon_softId")) Muon_LooseId                         = new TTreeReaderArray< bool >(*treeReader_, "Muon_softId");
  if (branchFound("Muon_mediumId")) Muon_MediumId                      = new TTreeReaderArray< bool >(*treeReader_, "Muon_mediumId");
  if (branchFound("Muon_tightId")) Muon_TightId                        = new TTreeReaderArray< bool >(*treeReader_, "Muon_tightId");
  if (branchFound("Muon_highPtId")) Muon_highPtId                      = new TTreeReaderArray< unsigned char >(*treeReader_, "Muon_highPtId");
  if (branchFound("Muon_mvaId")) Muon_mvaId                            = new TTreeReaderArray< unsigned char >(*treeReader_, "Muon_mvaId");
  if (branchFound("Muon_pfRelIso03_all")) Muon_pfRelIso03_all          = new TTreeReaderArray< float >(*treeReader_, "Muon_pfRelIso03_all");
  if (branchFound("Muon_pfRelIso03_chg")) Muon_pfRelIso03_chg          = new TTreeReaderArray< float >(*treeReader_, "Muon_pfRelIso03_chg");
  if (branchFound("Muon_pfRelIso04_all")) Muon_pfRelIso04_all          = new TTreeReaderArray< float >(*treeReader_, "Muon_pfRelIso04_all");
  if (isMC_) {
   if (branchFound("Muon_genPartIdx")) Muon_genPartIdx                 = new TTreeReaderArray< int >(*treeReader_, "Muon_genPartIdx");
   if (branchFound("Muon_genPartFlav")) Muon_genPartFlv                = new TTreeReaderArray< unsigned char >(*treeReader_, "Muon_genPartFlav");
  }
  if (branchFound("Muon_isGlobal")) Muon_isGlobal                      = new TTreeReaderArray< bool >(*treeReader_, "Muon_isGlobal");
  if (branchFound("Muon_isPFcand")) Muon_isPFcand                      = new TTreeReaderArray< bool >(*treeReader_, "Muon_isPFcand");
  if (branchFound("Muon_isTracker")) Muon_isTracker                    = new TTreeReaderArray< bool >(*treeReader_, "Muon_isTracker");
  if (branchFound("Muon_dxy")) Muon_dxy                                = new TTreeReaderArray< float >(*treeReader_, "Muon_dxy");
  if (branchFound("Muon_dxy")) Muon_dz                                 = new TTreeReaderArray< float >(*treeReader_, "Muon_dz");
  if (branchFound("Muon_tightCharge")) Muon_tightCharge                = new TTreeReaderArray< int >(*treeReader_, "Muon_tightCharge");
  if (branchFound("Muon_miniPFRelIso_all")) Muon_miniPFRelIso_all      = new TTreeReaderArray< float >(*treeReader_, "Muon_miniPFRelIso_all");

  //Electron
  if (branchFound("nElectron")) nElectron                              = new TTreeReaderValue< unsigned int >(*treeReader_, "nElectron");
  if (branchFound("Electron_pt")) Electron_pt                          = new TTreeReaderArray< float >(*treeReader_, "Electron_pt");
  if (branchFound("Electron_eta")) Electron_eta                        = new TTreeReaderArray< float >(*treeReader_, "Electron_eta");
  if (branchFound("Electron_phi")) Electron_phi                        = new TTreeReaderArray< float >(*treeReader_, "Electron_phi");
  if (branchFound("Electron_charge")) Electron_charge                  = new TTreeReaderArray< int >(*treeReader_, "Electron_charge");
  if (branchFound("Electron_mass")) Electron_mass                      = new TTreeReaderArray< float >(*treeReader_, "Electron_mass");
  if (branchFound("Electron_sip3d")) Electron_sip3d                    = new TTreeReaderArray< float >(*treeReader_, "Electron_sip3d");
  if (branchFound("Electron_jetIdx")) Electron_jetIdx                  = new TTreeReaderArray< int >(*treeReader_, "Electron_jetIdx");
  if (branchFound("Electron_photonIdx")) Electron_phoIdx               = new TTreeReaderArray< int >(*treeReader_, "Electron_photonIdx");
  if (branchFound("Electron_pfRelIso03_all")) Electron_pfRelIso03_all  = new TTreeReaderArray< float >(*treeReader_, "Electron_pfRelIso03_all");
  if (branchFound("Electron_pfRelIso03_chg")) Electron_pfRelIso03_chg  = new TTreeReaderArray< float >(*treeReader_, "Electron_pfRelIso03_chg");
  if (branchFound("Electron_mvaFall17V2Iso_WP80")) Electron_mvaFall17V2Iso_WP80 = new TTreeReaderArray< bool >(*treeReader_, "Electron_mvaFall17V2Iso_WP80");
  if (branchFound("Electron_mvaFall17V2Iso_WP90")) Electron_mvaFall17V2Iso_WP90 = new TTreeReaderArray< bool >(*treeReader_, "Electron_mvaFall17V2Iso_WP90");
  if (branchFound("Electron_mvaFall17V1Iso_WP90")) Electron_mvaFall17V1Iso_WP90 = new TTreeReaderArray< bool >(*treeReader_, "Electron_mvaFall17V1Iso_WP90");
  if (isMC_) {
    if (branchFound("Electron_genPartIdx")) Electron_genPartIdx        = new TTreeReaderArray< int >(*treeReader_, "Electron_genPartIdx");
    if (branchFound("Electron_genPartFlav")) Electron_genPartFlv       = new TTreeReaderArray< unsigned char >(*treeReader_, "Electron_genPartFlav");
  }
  if (branchFound("Electron_dxy")) Electron_dxy                        = new TTreeReaderArray< float >(*treeReader_, "Electron_dxy");
  if (branchFound("Electron_dz")) Electron_dz                          = new TTreeReaderArray< float >(*treeReader_, "Electron_dz");
  if (branchFound("Electron_mvaFall17V2noIso_WPL")) Electron_mvaFall17V2noIso_WPL = new TTreeReaderArray< bool >(*treeReader_, "Electron_mvaFall17V2noIso_WPL");
  if (branchFound("Electron_mvaFall17V1noIso_WPL")) Electron_mvaFall17V1noIso_WPL = new TTreeReaderArray< bool >(*treeReader_, "Electron_mvaFall17V1noIso_WPL");
  if (branchFound("Electron_lostHits")) Electron_lostHits              = new TTreeReaderArray< unsigned char >(*treeReader_, "Electron_lostHits");
  if (branchFound("Electron_convVeto")) Electron_convVeto              = new TTreeReaderArray< bool >(*treeReader_, "Electron_convVeto");
  if (branchFound("Electron_deltaEtaSC")) Electron_deltaEtaSC          = new TTreeReaderArray< float >(*treeReader_, "Electron_deltaEtaSC");

  //MET
  if (branchFound("MET_pt")) Met_pt                                    = new TTreeReaderValue< float >(*treeReader_, "MET_pt");
  if (branchFound("MET_pt_nom")) Met_nomPt                             = new TTreeReaderValue< float >(*treeReader_, "MET_pt_nom");
  if (branchFound("MET_phi")) Met_phi                                  = new TTreeReaderValue< float >(*treeReader_, "MET_phi");
  if (branchFound("MET_phi_nom")) Met_nomPhi                           = new TTreeReaderValue< float >(*treeReader_, "MET_phi_nom");
  if (branchFound("MET_significance")) Met_significance                = new TTreeReaderValue< float >(*treeReader_, "MET_significance");
  if (branchFound("MET_sumEt")) Met_sumEt                              = new TTreeReaderValue< float >(*treeReader_, "MET_sumEt");

  //AK4Jet
  if (branchFound("nJet")) nJet                                        = new TTreeReaderValue< unsigned int >(*treeReader_, "nJet");
  if (branchFound("Jet_pt")) Jet_pt                                    = new TTreeReaderArray< float >(*treeReader_, "Jet_pt");
  if (branchFound("Jet_pt_nom")) Jet_nomPt                             = new TTreeReaderArray< float >(*treeReader_, "Jet_pt_nom");
  if (branchFound("Jet_eta")) Jet_eta                                  = new TTreeReaderArray< float >(*treeReader_, "Jet_eta");
  if (branchFound("Jet_phi")) Jet_phi                                  = new TTreeReaderArray< float >(*treeReader_, "Jet_phi");
  if (branchFound("Jet_nConstituents")) Jet_nConstituents              = new TTreeReaderArray< int >(*treeReader_, "Jet_nConstituents");
  if (branchFound("Jet_nMuons")) Jet_nMuons                            = new TTreeReaderArray< int >(*treeReader_, "Jet_nMuons");
  if (branchFound("Jet_nElectrons")) Jet_nElectrons                    = new TTreeReaderArray< int >(*treeReader_, "Jet_nElectrons");
  if (branchFound("Jet_puId")) Jet_puId                                = new TTreeReaderArray< int >(*treeReader_, "Jet_puId");
  if (branchFound("Jet_muonIdx1")) Jet_muIdx1                          = new TTreeReaderArray< int >(*treeReader_, "Jet_muonIdx1");
  if (branchFound("Jet_muonIdx2")) Jet_muIdx2                          = new TTreeReaderArray< int >(*treeReader_, "Jet_muonIdx2");
  if (branchFound("Jet_electronIdx1")) Jet_elIdx1                      = new TTreeReaderArray< int >(*treeReader_, "Jet_electronIdx1");
  if (branchFound("Jet_electronIdx2")) Jet_elIdx2                      = new TTreeReaderArray< int >(*treeReader_, "Jet_electronIdx2");
  if (branchFound("Jet_jetId")) Jet_jetId                              = new TTreeReaderArray< int >(*treeReader_, "Jet_jetId");
  if (branchFound("Jet_qgl")) Jet_qgl                                  = new TTreeReaderArray< float >(*treeReader_, "Jet_qgl");
  if (branchFound("Jet_neHEF")) Jet_neHEF                              = new TTreeReaderArray< float >(*treeReader_, "Jet_neHEF");
  if (branchFound("Jet_neEmEF")) Jet_neEmEF                            = new TTreeReaderArray< float >(*treeReader_, "Jet_neEmEF");
  if (branchFound("Jet_chHEF")) Jet_chHEF                              = new TTreeReaderArray< float >(*treeReader_, "Jet_chHEF");
  if (branchFound("Jet_chEmEF")) Jet_chEmEF                            = new TTreeReaderArray< float >(*treeReader_, "Jet_chEmEF");
  if (branchFound("Jet_area")) Jet_area                                = new TTreeReaderArray< float >(*treeReader_, "Jet_area");
  if (branchFound("Jet_mass")) Jet_mass                                = new TTreeReaderArray< float >(*treeReader_, "Jet_mass");
  if (branchFound("Jet_mass_nom")) Jet_nomMass                         = new TTreeReaderArray< float >(*treeReader_, "Jet_mass_nom");
  if (branchFound("Jet_btagCSVV2")) Jet_btagCSVV2                      = new TTreeReaderArray< float >(*treeReader_, "Jet_btagCSVV2");
  if (branchFound("Jet_btagDeepFlavB")) Jet_btagDeepFlavB              = new TTreeReaderArray< float >(*treeReader_, "Jet_btagDeepFlavB");

  //FatJet
  if (branchFound("nFatJet")) nFatJet                                  = new TTreeReaderValue< unsigned int >(*treeReader_, "nFatJet");
  if (branchFound("FatJet_pt")) FatJet_pt                              = new TTreeReaderArray< float >(*treeReader_, "FatJet_pt");
  if (branchFound("FatJet_eta")) FatJet_eta                            = new TTreeReaderArray< float >(*treeReader_, "FatJet_eta");
  if (branchFound("FatJet_phi")) FatJet_phi                            = new TTreeReaderArray< float >(*treeReader_, "FatJet_phi");
  if (branchFound("FatJet_mass")) FatJet_mass                          = new TTreeReaderArray< float >(*treeReader_, "FatJet_mass");
  if (branchFound("FatJet_msoftdrop")) FatJet_msoftdrop                = new TTreeReaderArray< float >(*treeReader_, "FatJet_msoftdrop");
  if (branchFound("FatJet_jetId")) FatJet_jetId                        = new TTreeReaderArray< int >(*treeReader_, "FatJet_jetId"); 
  if (branchFound("FatJet_btagDeepB")) FatJet_btagDeepB                = new TTreeReaderArray< float >(*treeReader_, "FatJet_btagDeepB");
  if (branchFound("FatJet_btagCSVV2")) FatJet_btagCSVV2                = new TTreeReaderArray< float >(*treeReader_, "FatJet_btagCSVV2");
  if (branchFound("FatJet_n2b1")) FatJet_n2b1                          = new TTreeReaderArray< float >(*treeReader_, "FatJet_n2b1");
  if (branchFound("FatJet_n3b1")) FatJet_n3b1                          = new TTreeReaderArray< float >(*treeReader_, "FatJet_n3b1");
  if (branchFound("FatJet_hadronFlavour")) FatJet_hadronFlavour        = new TTreeReaderArray< int >(*treeReader_, "FatJet_hadronFlavour");
  if (branchFound("FatJet_nBHadrons")) FatJet_nBHadrons                = new TTreeReaderArray< unsigned char >(*treeReader_, "FatJet_nBHadrons");
  if (branchFound("FatJet_nCHadrons")) FatJet_nCHadrons                = new TTreeReaderArray< unsigned char >(*treeReader_, "FatJet_nCHadrons");
  if (branchFound("FatJet_rawFactor")) FatJet_rawFactor                = new TTreeReaderArray< float >(*treeReader_, "FatJet_rawFactor");
  if (branchFound("FatJet_subJetIdx1")) FatJet_subJetIdx1              = new TTreeReaderArray< int >(*treeReader_, "FatJet_subJetIdx1");
  if (branchFound("FatJet_subJetIdx2")) FatJet_subJetIdx2              = new TTreeReaderArray< int >(*treeReader_, "FatJet_subJetIdx2");
  if (branchFound("FatJet_tau1")) FatJet_tau1                          = new TTreeReaderArray< float >(*treeReader_, "FatJet_tau1");
  if (branchFound("FatJet_tau2")) FatJet_tau2                          = new TTreeReaderArray< float >(*treeReader_, "FatJet_tau2");
  if (branchFound("FatJet_tau3")) FatJet_tau3                          = new TTreeReaderArray< float >(*treeReader_, "FatJet_tau3");
  if (branchFound("FatJet_tau3")) FatJet_tau4                          = new TTreeReaderArray< float >(*treeReader_, "FatJet_tau4");
  if (branchFound("FatJet_deepTag_WvsQCD")) FatJet_deepTag_WvsQCD      = new TTreeReaderArray< float >(*treeReader_, "FatJet_deepTag_WvsQCD");
  if (branchFound("FatJet_deepTag_ZvsQCD")) FatJet_deepTag_ZvsQCD      = new TTreeReaderArray< float >(*treeReader_, "FatJet_deepTag_ZvsQCD");
  if (branchFound("FatJet_deepTag_TvsQCD")) FatJet_deepTag_TvsQCD      = new TTreeReaderArray< float >(*treeReader_, "FatJet_deepTag_TvsQCD");
  if (branchFound("FatJet_deepTagMD_WvsQCD")) FatJet_deepTagMD_WvsQCD  = new TTreeReaderArray< float >(*treeReader_, "FatJet_deepTagMD_WvsQCD");
  if (branchFound("FatJet_deepTagMD_ZvsQCD")) FatJet_deepTagMD_ZvsQCD  = new TTreeReaderArray< float >(*treeReader_, "FatJet_deepTagMD_ZvsQCD");
  if (branchFound("FatJet_deepTagMD_TvsQCD")) FatJet_deepTagMD_TvsQCD  = new TTreeReaderArray< float >(*treeReader_, "FatJet_deepTagMD_TvsQCD");
  if (branchFound("FatJet_electronIdx3SJ")) FatJet_electronIdx3SJ      = new TTreeReaderArray< int >(*treeReader_, "FatJet_electronIdx3SJ");
  if (branchFound("FatJet_muonIdx3SJ")) FatJet_muonIdx3SJ              = new TTreeReaderArray< int >(*treeReader_, "FatJet_muonIdx3SJ");

  // SubJet
  if (branchFound("nSubJet")) nSubJet                                  = new TTreeReaderValue< unsigned int >(*treeReader_, "nSubJet");
  if (branchFound("SubJet_pt")) SubJet_pt                              = new TTreeReaderArray< float >(*treeReader_, "SubJet_pt");
  if (branchFound("SubJet_eta")) SubJet_eta                            = new TTreeReaderArray< float >(*treeReader_, "SubJet_eta");
  if (branchFound("SubJet_phi")) SubJet_phi                            = new TTreeReaderArray< float >(*treeReader_, "SubJet_phi");
  if (branchFound("SubJet_mass")) SubJet_mass                          = new TTreeReaderArray< float >(*treeReader_, "SubJet_mass");
  if (branchFound("SubJet_btagDeepB")) SubJet_btagDeepB                = new TTreeReaderArray< float >(*treeReader_, "SubJet_btagDeepB");
  if (branchFound("SubJet_rawFactor")) SubJet_rawFactor                = new TTreeReaderArray< float >(*treeReader_, "SubJet_rawFactor");


  //Tau
  if (branchFound("nTau")) nTau                                        = new TTreeReaderValue< unsigned int >(*treeReader_, "nTau");
  if (branchFound("Tau_pt")) Tau_pt                                    = new TTreeReaderArray< float >(*treeReader_, "Tau_pt");
  if (branchFound("Tau_eta")) Tau_eta                                  = new TTreeReaderArray< float >(*treeReader_, "Tau_eta");
  if (branchFound("Tau_phi")) Tau_phi                                  = new TTreeReaderArray< float >(*treeReader_, "Tau_phi");
  if (branchFound("Tau_mass")) Tau_mass                                = new TTreeReaderArray< float >(*treeReader_, "Tau_mass");
  if (branchFound("Tau_dxy")) Tau_dxy                                  = new TTreeReaderArray< float >(*treeReader_, "Tau_dxy");
  if (branchFound("Tau_dz")) Tau_dz                                    = new TTreeReaderArray< float >(*treeReader_, "Tau_dz"); 
  if (branchFound("Tau_charge")) Tau_charge                            = new TTreeReaderArray< int >(*treeReader_, "Tau_charge");
  if (branchFound("Tau_decayMode")) Tau_decayMode                      = new TTreeReaderArray< int >(*treeReader_, "Tau_decayMode");
  if (branchFound("Tau_idDecayMode")) Tau_idDecayMode                  = new TTreeReaderArray< bool >(*treeReader_, "Tau_idDecayMode");
  if (branchFound("Tau_idDecayModeNewDMs")) Tau_idDecayModeNewDMs      = new TTreeReaderArray< bool >(*treeReader_, "Tau_idDecayModeNewDMs");
  if (branchFound("Tau_jetIdx")) Tau_jetIdx                            = new TTreeReaderArray< int >(*treeReader_, "Tau_jetIdx");
  if (branchFound("Tau_idAntiEle")) Tau_idAntiEle                      = new TTreeReaderArray< unsigned char >(*treeReader_, "Tau_idAntiEle");
  if (branchFound("Tau_idAntiMu")) Tau_idAntiMu                        = new TTreeReaderArray< unsigned char >(*treeReader_, "Tau_idAntiMu");
  if (branchFound("Tau_idMVAoldDM")) Tau_idMVAoldDM                    = new TTreeReaderArray< unsigned char >(*treeReader_, "Tau_idMVAoldDM");
  if (branchFound("Tau_idDeepTau2017v2VSjet")) Tau_idDeepTau2017v2VSjet= new TTreeReaderArray< unsigned char >(*treeReader_, "Tau_idDeepTau2017v2VSjet");
  if (branchFound("Tau_idDeepTau2017v2VSmu")) Tau_idDeepTau2017v2VSmu  = new TTreeReaderArray< unsigned char >(*treeReader_, "Tau_idDeepTau2017v2VSmu");
  if (branchFound("Tau_idDeepTau2017v2VSe")) Tau_idDeepTau2017v2VSe    = new TTreeReaderArray< unsigned char >(*treeReader_, "Tau_idDeepTau2017v2VSe");
  if (branchFound("Tau_idDeepTau2017v2p1VSjet")) Tau_idDeepTau2017v2p1VSjet= new TTreeReaderArray< unsigned char >(*treeReader_, "Tau_idDeepTau2017v2p1VSjet");
  if (branchFound("Tau_idDeepTau2017v2p1VSmu")) Tau_idDeepTau2017v2p1VSmu  = new TTreeReaderArray< unsigned char >(*treeReader_, "Tau_idDeepTau2017v2p1VSmu");
  if (branchFound("Tau_idDeepTau2017v2p1VSe")) Tau_idDeepTau2017v2p1VSe    = new TTreeReaderArray< unsigned char >(*treeReader_, "Tau_idDeepTau2017v2p1VSe");

  if (isMC_){
    //LHE
    if (branchFound("LHE_Njets")) LHEnJets                             = new TTreeReaderValue< unsigned char >(*treeReader_, "LHE_Njets");
    if (readGenInfo_){
      if (branchFound("nLHEPart")) nLHEPart                              = new TTreeReaderValue< unsigned int >(*treeReader_, "nLHEPart");
      if (branchFound("LHEPart_pt")) LHEPart_pt                          = new TTreeReaderArray< float >(*treeReader_, "LHEPart_pt");
      if (branchFound("LHEPart_eta")) LHEPart_eta                        = new TTreeReaderArray< float >(*treeReader_, "LHEPart_eta");
      if (branchFound("LHEPart_phi")) LHEPart_phi                        = new TTreeReaderArray< float >(*treeReader_, "LHEPart_phi");
      if (branchFound("LHEPart_mass")) LHEPart_mass                      = new TTreeReaderArray< float >(*treeReader_, "LHEPart_mass");
      if (branchFound("LHEPart_pdgId")) LHEPart_pdgId                    = new TTreeReaderArray< int >(*treeReader_, "LHEPart_pdgId");
      //GenParticle
      if (branchFound("nGenPart")) nGenPart                              = new TTreeReaderValue< unsigned int >(*treeReader_, "nGenPart");
      if (branchFound("GenPart_pt")) GenPart_pt                          = new TTreeReaderArray< float >(*treeReader_, "GenPart_pt");
      if (branchFound("GenPart_eta")) GenPart_eta                        = new TTreeReaderArray< float >(*treeReader_, "GenPart_eta");
      if (branchFound("GenPart_phi")) GenPart_phi                        = new TTreeReaderArray< float >(*treeReader_, "GenPart_phi");
      if (branchFound("GenPart_mass")) GenPart_mass                      = new TTreeReaderArray< float >(*treeReader_, "GenPart_mass");
      if (branchFound("GenPart_genPartIdxMother")) GenPart_motherIdx     = new TTreeReaderArray< int >(*treeReader_, "GenPart_genPartIdxMother");
      if (branchFound("GenPart_pdgId")) GenPart_pdgId                    = new TTreeReaderArray< int >(*treeReader_, "GenPart_pdgId");
      if (branchFound("GenPart_status")) GenPart_status                  = new TTreeReaderArray< int >(*treeReader_, "GenPart_status");
      if (branchFound("GenPart_statusFlags")) GenPart_statusFlags        = new TTreeReaderArray< int >(*treeReader_, "GenPart_statusFlags");
    }
  }
  // Getting the HLT branch pointers
  for (auto& hlt : doubleMuonHltPathList_) {
    bool found = branchFound(hlt.c_str());
    if (found) {
      TTreeReaderValue< bool >* HLT = new TTreeReaderValue< bool >(*treeReader_, hlt.c_str());
      doubleMuonHltPtrList_.push_back(HLT);
    }
  }
  for (auto& hlt : singleMuonHltPathList_) {
    bool found = branchFound(hlt.c_str());
    if (found) {
      TTreeReaderValue< bool >* HLT = new TTreeReaderValue< bool >(*treeReader_, hlt.c_str());
      singleMuonHltPtrList_.push_back(HLT);
    }
  }
  for (auto& hlt : doubleEgHltPathList_) {
    bool found = branchFound(hlt.c_str());
    if (found) {
      TTreeReaderValue< bool >* HLT = new TTreeReaderValue< bool >(*treeReader_, hlt.c_str());
      doubleEgHltPtrList_.push_back(HLT);
    }
  }
  for (auto& hlt : singleElectronHltPathList_) {
    bool found = branchFound(hlt.c_str());
    if (found) {
      TTreeReaderValue< bool >* HLT = new TTreeReaderValue< bool >(*treeReader_, hlt.c_str());
      singleElectronHltPtrList_.push_back(HLT);
    }
  }
  for (auto& hlt : muonEgHltPathList_) {
    bool found = branchFound(hlt.c_str());
    if (found) {
      TTreeReaderValue< bool >* HLT = new TTreeReaderValue< bool >(*treeReader_, hlt.c_str());
      muonEgHltPtrList_.push_back(HLT);
    }
  }
  
  return true;
}



void AnaBase::endJob() 
{
  //closeFiles();
}
void AnaBase::closeHistFile() //Called by closeFiles
{ 
  histf_->cd();
  histf_->Write();
  histf_->Close();
  //fakehistf_->cd();
  //fakehistf_->Write();
  //fakehistf_->Close();
}
double AnaBase::lumiWt(double evtWeightSum, bool verbose) const 
{
  double nevt = (evtWeightSum > -1) ? evtWeightSum : AnaUtil::cutValue(lumiWtMap(), "nevents");
  if (verbose) 
    std::cout << "-- intLumi: " << AnaUtil::cutValue(lumiWtMap(), "intLumi") 
	      << " xsec: " << AnaUtil::cutValue(lumiWtMap(), "xsec") 
	      << " nevt: " << nevt << std::endl; 
  return (AnaUtil::cutValue(lumiWtMap(), "intLumi") * AnaUtil::cutValue(lumiWtMap(), "xsec") / nevt); 
}

// ---------------------------------
// Add input Root files to the chain
// ---------------------------------
int AnaBase::setInputFile(const string& fname) //Called by readJob
{
  auto found = fname.find("root:");
  if (found == string::npos && gSystem->AccessPathName(fname.c_str())) {
    cerr << ">>> Warning: File <<" << fname << ">> was not found!!" << endl;
    return static_cast<int>(chain_->GetEntries()); 
  }
  chain_->AddFile(fname.c_str(), -1);
  chainRun_ ->AddFile(fname.c_str(), -1);
  return static_cast<int>(chain_->GetEntries()); 
}
// ---------------------------------------
// Get total number of events in the chain
// --------------------------------------
int AnaBase::getEntries() const 
{
  return static_cast<int>(chain_->GetEntriesFast());
}

// ------------------------------------------------------
// Open the output file with a global filehandle, C++ way
// ------------------------------------------------------
bool AnaBase::openFiles() //Called by beginJob  
{
  // not using the dump funtion right now.
  // so commented out
  /*
  fLog_.open(logFile_.c_str(), ios::out);
  if (!fLog_) {
    cerr << "File: " << logFile_ << " could not be opened!" << endl;
    return false;
  }
  fLog_ << setiosflags(ios::fixed);
  */
  evLog_.open(evFile_.c_str(), ios::out);
  if (!evLog_) {
    cerr << "File: " << evFile_ << " could not be opened!" << endl;
    return false;
  }
  evLog_ << setiosflags(ios::fixed);

  selEvLog_.open(selEvFile_.c_str(), ios::out);
  if (!selEvLog_) {
    cerr << "File: " << selEvFile_ << " could not be opened!" << endl;
    return false;
  }
  selEvLog_ << setiosflags(ios::fixed);

  return true;
}
// ------------------------
// Close the output file
// ------------------------
void AnaBase::closeFiles() 
{
  /*
  if (fLog_) {
    fLog_ << resetiosflags(ios::fixed); 
    fLog_.close();
  }
  */
  if (evLog_) {
    evLog_ << resetiosflags(ios::fixed); 
    evLog_.close();
  }
  if (selEvLog_) {
    selEvLog_ << resetiosflags(ios::fixed); 
    selEvLog_.close();
  }
  closeHistFile();
  if (chain_) delete chain_;
  if (chainRun_) delete chainRun_;
  if (treeReader_) delete treeReader_;
  if (treeReaderRun_) delete treeReaderRun_;
}

bool AnaBase::branchFound(const string& b)
{
  TBranch* branch1 = chain_->GetBranch(b.c_str());  // Get branch pointer
  TBranch* branch2 = chainRun_->GetBranch(b.c_str());  // Get branch pointer
  if (branch1 == nullptr && branch2 == nullptr) {
    cout << ">>> SetBranchAddress: <<<" << b << ">>> N O T   F O U N D !!!" << endl;
    return false;
  }
  cout << ">>> SetBranchAddress: <" << b << "> found!" << endl;
  brList_.push_back(b);
  return true;
}
int AnaBase::getEntry(int lflag) const
{
  int nbytes = 0;
  for (const auto& v: brList_) {
    TBranch* branch = chain_->GetBranch(v.c_str());
    if (branch == nullptr) {
      cout << ">>> Branch: " << v << " not found!" << endl;
      continue;
    }
    nbytes += branch->GetEntry(lflag);
  }
  return nbytes;
}
// not used yet
void AnaBase::enableBranches() 
{
  chain_->SetBranchStatus("*", kFALSE); // Disable all branches
  for (const auto& v: brList_)
    chain_->SetBranchStatus(v.c_str(), kTRUE);
}
bool AnaBase::readJob(const string& jobFile, int& nFiles)
{
  // Open the file containing the datacards
  ifstream fin(jobFile.c_str(), ios::in);    
  if (!fin) {
    cerr << "==> Input File: <<" << jobFile << ">> could not be opened!" << endl;
    return false;
  }

  static constexpr int BUF_SIZE = 256;
  char buf[BUF_SIZE];
  while (fin.getline(buf, BUF_SIZE, '\n')) {  // Pops off the newline character
    string line(buf);
    if (line.empty() || line == "START") continue;   

    // enable '#' and '//' style comments
    if (line.substr(0,1) == "#" || line.substr(0,2) == "//") continue;
    if (line == "END") break;

    // Split the line into words
    vector<string> tokens;
    AnaUtil::tokenize(line, tokens);
    int vsize = tokens.size();
    assert(vsize > 1);

    const string& key   = tokens.at(0);
    const string& value = tokens.at(1);
    if (key == "dataType") {
      string vtmp(value);
      std::transform(vtmp.begin(), vtmp.end(), vtmp.begin(), ::toupper);
      vector<string> dt;
      AnaUtil::tokenize(vtmp, dt, "#");
      if (dt.size()) {
        isMC_ = (dt.at(0) == "MC") ? true : false;
        if (isMC_ && dt.size() > 1) {
          isSignal_ = (dt.at(1) == "SIGNAL") ? true : false;
        }
      }
    }
    else if (key == "era") 
      era_ = std::stoi(value.c_str());    
    else if (key == "evtWtSum")
      evtWtSum_ = value;
    else if (key == "dataset")
      dataset_ = value;
    else if (key == "readTrigObject") 
      readTrigObject_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "readGenInfo") 
      readGenInfo_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "readPFObject") 
      readPFObject_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "useTrigger") 
      useTrigger_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "useLumiWt") 
      useLumiWt_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "usePUWt") 
      usePUWt_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "logFile")
      logFile_ = value;
    else if (key == "eventFile")
      evFile_  = value;
    else if (key == "selEventFile")
      selEvFile_  = value;
    else if (key == "logOption") 
      logOption_ = strtol(value.c_str(), NULL, 2);
    else if (key == "maxEvent") 
      maxEvt_ = std::stoi(value.c_str());
    else if (key == "nFiles") 
      nFiles_ = std::stoi(value.c_str());
    else if (key == "startEvent") 
      firstEvt_ = std::stoi(value.c_str());
    else if (key == "endEvent") 
      lastEvt_ = std::stoi(value.c_str());
    else if (key == "bunchX") 
      bunchCrossing_ = std::stoi(value.c_str());
    else if (key == "histFile") 
      histFile_ = value;
    //else if (key == "fakehistFile") 
    //  fakehistFile_ = value;
    else if (key == "puHistFile") 
      puHistFile_ = value;
    else if (key == "puHistogram") 
      puHistogram_ = value;
    else if (key == "useTrueNInt") 
      useTrueNInt_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "inputFile") 
      AnaUtil::buildList(tokens, fileList_);
    else if (key == "muonIdSFRootFile")
      muonIdSFRootFile_ = value;
    else if (key == "looseMuonIdSFhistName")
      looseMuonIdSFhistName_ = value;
    else if (key == "medMuonIdSFhistName")
      medMuonIdSFhistName_ = value;
    else if (key == "tightMuonIdSFhistName")
      tightMuonIdSFhistName_ = value;
    else if (key == "electronLooseIdSFRootFile")
      electronLooseIdSFRootFile_ = value;
    else if (key == "looseEleIdSFhistName")
      looseEleIdSFhistName_ = value;
    else if (key == "electronTightIdSFRootFile")
      electronTightIdSFRootFile_ = value;
    else if (key == "tightEleIdSFhistName")
      tightEleIdSFhistName_ = value;
    else if (key == "muonTightIsoSFRootFile")
      muonTightIsoSFRootFile_ = value;
    else if (key == "tightMuIsoSFhistName")
      tightMuIsoSFhistName_ = value;
    else if (key == "DoubleMuon") 
      AnaUtil::buildList(tokens, doubleMuonHltPathList_);
    else if (key == "SingleMuon") 
      AnaUtil::buildList(tokens, singleMuonHltPathList_);
    else if (key == "DoubleEG") 
      AnaUtil::buildList(tokens, doubleEgHltPathList_);
    else if (key == "SingleElectron") 
      AnaUtil::buildList(tokens, singleElectronHltPathList_);
    else if (key == "MuonEG") 
      AnaUtil::buildList(tokens, muonEgHltPathList_);
    else if (key == "eventId" && tokens.size() == 4) 
      AnaUtil::buildMap(tokens, eventIdMap_);
    else {
      if (0) cout << "==> " << line << endl;
      AnaUtil::storeCuts(tokens, hmap_);
    }
  }
  // Close the file
  fin.close();

  if (!isMC_) usePUWt_ = false;
  if (!isSignal_) readGenInfo_ = false;

  // Build the chain of root files
  for (const auto& fname: fileList_) {
    cout << ">>> INFO. Adding input file " << fname << " to TChain " << endl;
    ++nFiles;
    int nevt = setInputFile(fname);
    std::cout<<">>>Right now the chain contains :: "<<nevt<<" events\n";
    if (nFiles_ > 0 && nFiles > nFiles_) break;
    //    if (maxEvt_ > 0 && nevt >= maxEvt_) break;
  }
  if (!nFiles) {
    cerr << ">>> WARN. Input Root file list is empty! exiting ..." << endl;
    return false;
  }

  return true;
}
void AnaBase::printJob(ostream& os) const
{
  os << "       datatype: " << ((isMC_) ? "mc" : "data") << endl
     << "    readGenInfo: " << readGenInfo_ << endl
     << "        logFile: " << logFile_ << endl 
     << "      eventFile: " << evFile_ << endl
     << "       histFile: " << histFile_ << endl
     << "      selEvFile: " << selEvFile_ << endl
     << "      useLumiWt: " << std::boolalpha << useLumiWt_ << endl
     << "        usePUWt: " << std::boolalpha << usePUWt_ << endl
     << "     puHistFile: " << puHistFile_ << endl
     << "    puHistogram: " << puHistogram_ << endl
     << "    useTrueNInt: " << std::boolalpha << useTrueNInt_ << endl
     << "     useTrigger: " << std::boolalpha << useTrigger_ << endl
     << "      logOption: " << logOption_ << endl
     << "       maxEvent: " << maxEvt_ 
     << endl;

  // InputFiles
  if (chain_ != nullptr) {
    TObjArray *fileElements = chain_->GetListOfFiles();
    os << ">>> INFO. nFiles: " << fileElements->GetEntries() 
       << ", Files to analyse:" 
       << endl;
    TIter next(fileElements);
    TChainElement *chEl = 0;
    while (( chEl = dynamic_cast<TChainElement*>(next()) ))
      os << chEl->GetTitle() 
         << endl;
  }
  else
    AnaUtil::showList(fileList_, ">>> INFO. inputFiles:", os);

  // EventID 
  AnaUtil::showMap<string, int>(eventIdMap_, "Event List:", os);

  // Cuts
  AnaUtil::showCuts(hmap_, os);
}

bool AnaBase::openScaleFactorRootFiles() {
  // MUON POG ID SCALE FACTOR
  const char* fname_MuonIdSF = gSystem->ExpandPathName(muonIdSFRootFile_.c_str());
  if (gSystem->AccessPathName(fname_MuonIdSF)) {
    cerr << ">>> Warning: File <<" << muonIdSFRootFile_ << ">> not found!!" << endl;
    return false;
  }
  TFile* file_MuonIdSF = TFile::Open(fname_MuonIdSF);

  // --- Loose Muon ID SF
  file_MuonIdSF -> GetObject(looseMuonIdSFhistName_.c_str(), looseMuonIdSFhist_);
  if (!looseMuonIdSFhist_) {
    cerr << ">>> Warning: Histogram <<" << looseMuonIdSFhistName_ << ">> not found!!" << endl;
    return false;
  }
  looseMuonIdSFhist_->SetDirectory(0);

  // --- Medium Muon ID SF
  file_MuonIdSF ->GetObject(medMuonIdSFhistName_.c_str(), medMuonIdSFhist_);
  if (!medMuonIdSFhist_) {
    cerr << ">>> Warning: Histogram <<" << medMuonIdSFhistName_ << ">> not found!!" << endl;
    return false;
  }
  medMuonIdSFhist_->SetDirectory(0);

  // --- Tight Muon ID SF
  file_MuonIdSF ->GetObject(tightMuonIdSFhistName_.c_str(), tightMuonIdSFhist_);
  if (!tightMuonIdSFhist_) {
    cerr << ">>> Warning: Histogram <<" << tightMuonIdSFhist_ << ">> not found!!" << endl;
    return false;
  }
  tightMuonIdSFhist_->SetDirectory(0);

  file_MuonIdSF ->Close();
  delete file_MuonIdSF;

  // Electron POG ID SCALE FACTOR
  const char* fname_EleLooseIdSF = gSystem->ExpandPathName(electronLooseIdSFRootFile_.c_str());
  if (gSystem->AccessPathName(fname_EleLooseIdSF)) {
    cerr << ">>> Warning: File <<" << electronLooseIdSFRootFile_ << ">> not found!!" << endl;
    return false;
  }
  TFile* file_EleLooseIdSF = TFile::Open(fname_EleLooseIdSF);

  file_EleLooseIdSF -> GetObject(looseEleIdSFhistName_.c_str(), looseEleIdSFhist_);
  if (!looseEleIdSFhist_) {
    cerr << ">>> Warning: Histogram <<" << looseEleIdSFhist_ << ">> not found!!" << endl;
    return false;
  }
  looseEleIdSFhist_->SetDirectory(0);

  file_EleLooseIdSF ->Close();
  delete file_EleLooseIdSF;


  const char* fname_EleTightIdSF = gSystem->ExpandPathName(electronTightIdSFRootFile_.c_str());
  if (gSystem->AccessPathName(fname_EleTightIdSF)) {
    cerr << ">>> Warning: File <<" << electronTightIdSFRootFile_ << ">> not found!!" << endl;
    return false;
  }
  TFile* file_EleTightIdSF = TFile::Open(fname_EleTightIdSF);

  file_EleTightIdSF -> GetObject(tightEleIdSFhistName_.c_str(), tightEleIdSFhist_);
  if (!tightEleIdSFhist_) {
    cerr << ">>> Warning: Histogram <<" << tightEleIdSFhist_ << ">> not found!!" << endl;
    return false;
  }
  tightEleIdSFhist_->SetDirectory(0);

  file_EleTightIdSF ->Close();
  delete file_EleTightIdSF;

  // Iso Scale factor
  const char* fname_MuonTightIsoSF = gSystem->ExpandPathName(muonTightIsoSFRootFile_.c_str());
  if (gSystem->AccessPathName(fname_MuonTightIsoSF)) {
    cerr << ">>> Warning: File <<" << muonTightIsoSFRootFile_ << ">> not found!!" << endl;
    return false;
  }
  TFile* file_MuTightIsoSF = TFile::Open(fname_MuonTightIsoSF);

  file_MuTightIsoSF -> GetObject(tightMuIsoSFhistName_.c_str(), tightMuIsoSFhist_);
  if (!tightMuIsoSFhist_) {
    cerr << ">>> Warning: Histogram <<" << tightMuIsoSFhist_ << ">> not found!!" << endl;
    return false;
  }
  tightMuIsoSFhist_->SetDirectory(0);

  file_MuTightIsoSF ->Close();
  delete file_MuTightIsoSF;

  return true;
}
/*
double AnaBase::wtPileUp(float nPU, bool verbose) const {
  return puHist_->GetBinContent(puHist_->GetXaxis()->FindBin(nPU));
}
*/
double AnaBase::getIdSF(std::string IdType, float pt, float eta, std::string Flav) const {
  double SF = 0.0;
  if (Flav == "Muon") {
    if (IdType == "Loose") {
      if (pt < looseMuonIdSFhist_->GetXaxis()->GetXmax()){
	Int_t binX = looseMuonIdSFhist_->GetXaxis()->FindBin(pt);
	Int_t binY = looseMuonIdSFhist_->GetYaxis()->FindBin(std::abs(eta));
	SF = looseMuonIdSFhist_->GetBinContent(binX, binY);
      }
      else SF = 1.0;
    }
    else if (IdType == "Medium") {
      if (pt < medMuonIdSFhist_->GetXaxis()->GetXmax()){
	Int_t binX = medMuonIdSFhist_->GetXaxis()->FindBin(pt);
	Int_t binY = medMuonIdSFhist_->GetYaxis()->FindBin(std::abs(eta));
	SF = medMuonIdSFhist_->GetBinContent(binX, binY);
      }
      else SF = 1.0;
    }
    else if (IdType == "Tight") {
      if (pt < tightMuonIdSFhist_->GetXaxis()->GetXmax()){
	Int_t binX = tightMuonIdSFhist_->GetXaxis()->FindBin(pt);
	Int_t binY = tightMuonIdSFhist_->GetYaxis()->FindBin(std::abs(eta));
	SF = tightMuonIdSFhist_->GetBinContent(binX, binY);
      }
    }
  }
  else if (Flav == "Electron") {
    if (IdType == "Loose") {
      if (pt < looseEleIdSFhist_->GetYaxis()->GetXmax()){
	Int_t binX = looseEleIdSFhist_->GetXaxis()->FindBin(std::abs(eta));
	Int_t binY = looseEleIdSFhist_->GetYaxis()->FindBin(pt);
	SF = looseEleIdSFhist_->GetBinContent(binX, binY);
      }
      else SF = 1.00;
    }
    else if (IdType == "Tight") {
      if (pt < tightEleIdSFhist_->GetYaxis()->GetXmax()){
	Int_t binX = tightEleIdSFhist_->GetXaxis()->FindBin(std::abs(eta));
	Int_t binY = tightEleIdSFhist_->GetYaxis()->FindBin(pt);
	SF = tightEleIdSFhist_->GetBinContent(binX, binY);
      }
      else SF = 1.0;
    }
  }

  return SF;
}

double AnaBase::getIsoSF(std::string IsoType, float pt, float eta, std::string Flav) const {
  double SF = 0.0;
  if (Flav == "Muon") {
    if (IsoType == "Tight") {
      if (pt < tightMuIsoSFhist_->GetXaxis()->GetXmax()){
        Int_t binX = tightMuIsoSFhist_->GetXaxis()->FindBin(pt);
        Int_t binY = tightMuIsoSFhist_->GetYaxis()->FindBin(std::abs(eta));
        SF = tightMuIsoSFhist_->GetBinContent(binX, binY);
      }
      else SF = 1.0;
    }
  }
  return SF;
}
