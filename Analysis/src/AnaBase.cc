#include "configana.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <sstream>

#include "AnaBase.h"

#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TFile.h"

using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::stoi;

using std::ostream;
using std::string;
using std::vector;
using std::map;
using std::pair;
using std::setprecision;
using std::setw;
using std::setiosflags;
using std::resetiosflags;
using std::make_unique;

// -----------
// Constructor
// -----------
AnaBase::AnaBase() :
  chain_(make_unique<TChain>("Events")),
  chainRun_(make_unique<TChain>("Runs"))
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
}
// ----------
// Destructor
// ----------
AnaBase::~AnaBase() 
{
}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool AnaBase::beginJob() 
{ 
  if (isMC()) {
    if (!SFHandler_.openRootFiles()) {
      cerr << "==> ERROR. Scale Factors are not accessible!!!" << endl;
      return false;
    }
  }
  
  // Open the output ROOT file
  histf_ = make_unique<TFile>(histFile_.c_str(), "RECREATE");

  // Open the output FakeExtrapolated ROOT file
  //fakehistf_ = make_unique<TFile>(fakehistFile_.c_str(), "RECREATE");

  nEvents_ = static_cast<int>(chain_->GetEntries()); 
  if (nEvents_ <= 0) {
    cerr << "==> ERROR. No events found in the chain!" << endl;
    return false;
  }
  if (maxEvt_ > 0) nEvents_ = std::min(nEvents_, maxEvt_);
  cout << "==> INFO. nEvents = " << nEvents_ << endl;

  cout << "==> INFO. Create TTreeReader from TChain - Events" 
       << endl
       << "Chain Name: " << chain_->GetName()
       << ", Chain Entries: " << chain_->GetEntries() 
       << endl;
  treeReader_ = make_unique<TTreeReader>(chain_.get());
  assert(treeReader_);
  treeReader_->ls();

  cout << "==> INFO. Create TTreeReader with TChain - Runs" 
       << endl
       << "Chain Name: " << chainRun_->GetName()
       << ", Chain Entries: " << chainRun_->GetEntries() 
       << endl;
  treeReaderRun_ = make_unique<TTreeReader>(chainRun_.get());
  assert(treeReaderRun_);
  treeReaderRun_->ls();

  // Set Branch Addresses
  if (!init()) {
    cout << "==> ERROR. SetBranchAddress failed!!!" << endl;
    return false;
  }

  openFiles();

  return true;
}
// ---------------------------
// Initialize all the branches
// ---------------------------
bool AnaBase::init() {
  if (isMC_) {
    // from TTree::Runs
    if (branchFound(evtWtSum_.c_str(), "Runs")) genEventSumw_ = make_unique<TTreeReaderValue<double>>(*treeReaderRun_, evtWtSum_.c_str());

    // from TTree::Events (default in branchFound)
    if (branchFound("genWeight"))    genEvWt       = make_unique<TTreeReaderValue<float>>(*treeReader_, "genWeight");
    if (branchFound("puWeight"))     PU_Weight     = make_unique<TTreeReaderValue<float>>(*treeReader_, "puWeight");
    if (branchFound("puWeightUp"))   PU_WeightUp   = make_unique<TTreeReaderValue<float>>(*treeReader_, "puWeightUp");
    if (branchFound("puWeightDown")) PU_WeightDown = make_unique<TTreeReaderValue<float>>(*treeReader_, "puWeightDown");
  
    if (branchFound("btagWeight_CSVV2")) btagWeight_CSVV2 = make_unique<TTreeReaderValue<float>>(*treeReader_, "btagWeight_CSVV2");
  }

  // TTree::Events
  if (branchFound("run"))             run_   = make_unique<TTreeReaderValue<unsigned int>>(*treeReader_, "run");
  if (branchFound("event"))           event_ = make_unique<TTreeReaderValue<unsigned long long>>(*treeReader_, "event");
  if (branchFound("luminosityBlock")) lumis_ = make_unique<TTreeReaderValue<unsigned int>>(*treeReader_, "luminosityBlock");

  // Primary Vertices
  if (branchFound("PV_ndof"))     PV_ndf   = make_unique<TTreeReaderValue<float>>(*treeReader_, "PV_ndof");
  if (branchFound("PV_x"))        PV_xPos  = make_unique<TTreeReaderValue<float>>(*treeReader_, "PV_x");
  if (branchFound("PV_y"))        PV_yPos  = make_unique<TTreeReaderValue<float>>(*treeReader_, "PV_y");
  if (branchFound("PV_z"))        PV_zPos  = make_unique<TTreeReaderValue<float>>(*treeReader_, "PV_z");
  if (branchFound("PV_chi2"))     PV_chi2  = make_unique<TTreeReaderValue<float>>(*treeReader_, "PV_chi2");
  if (branchFound("PV_score"))    PV_score = make_unique<TTreeReaderValue<float>>(*treeReader_, "PV_score");
  if (branchFound("PV_npvs"))     nPV      = make_unique<TTreeReaderValue<int>>(*treeReader_, "PV_npvs");
  if (branchFound("PV_npvsGood")) nGoodPV  = make_unique<TTreeReaderValue<int>>(*treeReader_, "PV_npvsGood");

  // Muon
  if (branchFound("nMuon"))               nMuon               = make_unique<TTreeReaderValue<unsigned int>>(*treeReader_, "nMuon");
  if (branchFound("Muon_pt"))             Muon_pt             = make_unique<TTreeReaderArray<float>>(*treeReader_, "Muon_pt");
  if (branchFound("Muon_corrected_pt"))   Muon_corrpt         = make_unique<TTreeReaderArray<float>>(*treeReader_, "Muon_corrected_pt");
  if (branchFound("Muon_eta"))            Muon_eta            = make_unique<TTreeReaderArray<float>>(*treeReader_, "Muon_eta");
  if (branchFound("Muon_phi"))            Muon_phi            = make_unique<TTreeReaderArray<float>>(*treeReader_, "Muon_phi");
  if (branchFound("Muon_charge"))         Muon_charge         = make_unique<TTreeReaderArray<int>>(*treeReader_, "Muon_charge");
  if (branchFound("Muon_mass"))           Muon_mass           = make_unique<TTreeReaderArray<float>>(*treeReader_, "Muon_mass");
  if (branchFound("Muon_sip3d"))          Muon_sip3d          = make_unique<TTreeReaderArray<float>>(*treeReader_, "Muon_sip3d");
  if (branchFound("Muon_jetIdx"))         Muon_jetIdx         = make_unique<TTreeReaderArray<int>>(*treeReader_, "Muon_jetIdx");
  if (branchFound("Muon_softId"))         Muon_LooseId        = make_unique<TTreeReaderArray<bool>>(*treeReader_, "Muon_softId");
  if (branchFound("Muon_mediumId"))       Muon_MediumId       = make_unique<TTreeReaderArray<bool>>(*treeReader_, "Muon_mediumId");
  if (branchFound("Muon_tightId"))        Muon_TightId        = make_unique<TTreeReaderArray<bool>>(*treeReader_, "Muon_tightId");
  if (branchFound("Muon_highPtId"))       Muon_highPtId       = make_unique<TTreeReaderArray<unsigned char>>(*treeReader_, "Muon_highPtId");
  if (branchFound("Muon_mvaId"))          Muon_mvaId          = make_unique<TTreeReaderArray<unsigned char>>(*treeReader_, "Muon_mvaId");
  if (branchFound("Muon_pfRelIso03_all")) Muon_pfRelIso03_all = make_unique<TTreeReaderArray<float>>(*treeReader_, "Muon_pfRelIso03_all");
  if (branchFound("Muon_pfRelIso03_chg")) Muon_pfRelIso03_chg = make_unique<TTreeReaderArray<float>>(*treeReader_, "Muon_pfRelIso03_chg");
  if (branchFound("Muon_pfRelIso04_all")) Muon_pfRelIso04_all = make_unique<TTreeReaderArray<float>>(*treeReader_, "Muon_pfRelIso04_all");
  if (isMC_ && readGenInfo_) {
    if (branchFound("Muon_genPartIdx"))   Muon_genPartIdx     = make_unique<TTreeReaderArray<int>>(*treeReader_, "Muon_genPartIdx");
    if (branchFound("Muon_genPartFlav"))  Muon_genPartFlv     = make_unique<TTreeReaderArray<unsigned char>>(*treeReader_, "Muon_genPartFlav");
  }
  if (branchFound("Muon_isGlobal"))       Muon_isGlobal       = make_unique<TTreeReaderArray<bool>>(*treeReader_, "Muon_isGlobal");
  if (branchFound("Muon_isPFcand"))       Muon_isPFcand       = make_unique<TTreeReaderArray<bool>>(*treeReader_, "Muon_isPFcand");
  if (branchFound("Muon_isTracker"))      Muon_isTracker      = make_unique<TTreeReaderArray<bool>>(*treeReader_, "Muon_isTracker");
  if (branchFound("Muon_dxy"))            Muon_dxy            = make_unique<TTreeReaderArray<float>>(*treeReader_, "Muon_dxy");
  if (branchFound("Muon_dz"))             Muon_dz             = make_unique<TTreeReaderArray<float>>(*treeReader_, "Muon_dz");
  if (branchFound("Muon_tightCharge"))    Muon_tightCharge    = make_unique<TTreeReaderArray<int>>(*treeReader_, "Muon_tightCharge");
  if (branchFound("Muon_miniPFRelIso_all")) Muon_miniPFRelIso_all = make_unique<TTreeReaderArray<float>>(*treeReader_, "Muon_miniPFRelIso_all");

  // Electron
  if (branchFound("nElectron"))                   nElectron               = make_unique<TTreeReaderValue<unsigned int>>(*treeReader_, "nElectron");
  if (branchFound("Electron_pt"))                 Electron_pt             = make_unique<TTreeReaderArray<float>>(*treeReader_, "Electron_pt");
  if (branchFound("Electron_eta"))                Electron_eta            = make_unique<TTreeReaderArray<float>>(*treeReader_, "Electron_eta");
  if (branchFound("Electron_phi"))                Electron_phi            = make_unique<TTreeReaderArray<float>>(*treeReader_, "Electron_phi");
  if (branchFound("Electron_charge"))             Electron_charge         = make_unique<TTreeReaderArray<int>>(*treeReader_, "Electron_charge");
  if (branchFound("Electron_mass"))               Electron_mass           = make_unique<TTreeReaderArray<float>>(*treeReader_, "Electron_mass");
  if (branchFound("Electron_sip3d"))              Electron_sip3d          = make_unique<TTreeReaderArray<float>>(*treeReader_, "Electron_sip3d");
  if (branchFound("Electron_jetIdx"))             Electron_jetIdx         = make_unique<TTreeReaderArray<int>>(*treeReader_, "Electron_jetIdx");
  if (branchFound("Electron_photonIdx"))          Electron_phoIdx         = make_unique<TTreeReaderArray<int>>(*treeReader_, "Electron_photonIdx");
  if (branchFound("Electron_pfRelIso03_all"))     Electron_pfRelIso03_all = make_unique<TTreeReaderArray<float>>(*treeReader_, "Electron_pfRelIso03_all");
  if (branchFound("Electron_pfRelIso03_chg"))     Electron_pfRelIso03_chg = make_unique<TTreeReaderArray<float>>(*treeReader_, "Electron_pfRelIso03_chg");
  if (branchFound("Electron_mvaFall17V2Iso_WP80"))
                   Electron_mvaFall17V2Iso_WP80 = make_unique<TTreeReaderArray<bool>>(*treeReader_, "Electron_mvaFall17V2Iso_WP80");
  if (branchFound("Electron_mvaFall17V2Iso_WP90"))
                   Electron_mvaFall17V2Iso_WP90 = make_unique<TTreeReaderArray<bool>>(*treeReader_, "Electron_mvaFall17V2Iso_WP90");  
  if (branchFound("Electron_mvaFall17V1Iso_WP90"))
                   Electron_mvaFall17V1Iso_WP90 = make_unique<TTreeReaderArray<bool>>(*treeReader_, "Electron_mvaFall17V1Iso_WP90");

  if (isMC_ && readGenInfo_) {
    if (branchFound("Electron_genPartIdx"))  Electron_genPartIdx = make_unique<TTreeReaderArray<int>>(*treeReader_, "Electron_genPartIdx");
    if (branchFound("Electron_genPartFlav")) Electron_genPartFlv = make_unique<TTreeReaderArray<unsigned char>>(*treeReader_, "Electron_genPartFlav");
  }
  if (branchFound("Electron_dxy")) Electron_dxy = make_unique<TTreeReaderArray<float>>(*treeReader_, "Electron_dxy");
  if (branchFound("Electron_dz"))  Electron_dz  = make_unique<TTreeReaderArray<float>>(*treeReader_, "Electron_dz");
  if (branchFound("Electron_mvaFall17V2noIso_WPL"))
                   Electron_mvaFall17V2noIso_WPL = make_unique<TTreeReaderArray<bool>>(*treeReader_, "Electron_mvaFall17V2noIso_WPL");
  if (branchFound("Electron_mvaFall17V1noIso_WPL"))
                   Electron_mvaFall17V1noIso_WPL = make_unique<TTreeReaderArray<bool>>(*treeReader_, "Electron_mvaFall17V1noIso_WPL");
  if (branchFound("Electron_lostHits")) Electron_lostHits = make_unique<TTreeReaderArray<unsigned char>>(*treeReader_, "Electron_lostHits");
  if (branchFound("Electron_convVeto")) Electron_convVeto = make_unique<TTreeReaderArray<bool>>(*treeReader_, "Electron_convVeto");
  if (branchFound("Electron_deltaEtaSC")) Electron_deltaEtaSC = make_unique<TTreeReaderArray<float>>(*treeReader_, "Electron_deltaEtaSC");

  // MET
  if (branchFound("MET_pt"))           Met_pt           = make_unique<TTreeReaderValue<float>>(*treeReader_, "MET_pt");
  if (branchFound("MET_pt_nom"))       Met_nomPt        = make_unique<TTreeReaderValue<float>>(*treeReader_, "MET_pt_nom");
  if (branchFound("MET_phi"))          Met_phi          = make_unique<TTreeReaderValue<float>>(*treeReader_, "MET_phi");
  if (branchFound("MET_phi_nom"))      Met_nomPhi       = make_unique<TTreeReaderValue<float>>(*treeReader_, "MET_phi_nom");
  if (branchFound("MET_significance")) Met_significance = make_unique<TTreeReaderValue<float>>(*treeReader_, "MET_significance");
  if (branchFound("MET_sumEt"))        Met_sumEt        = make_unique<TTreeReaderValue<float>>(*treeReader_, "MET_sumEt");

  // AK4Jet
  if (branchFound("nJet"))              nJet              = make_unique<TTreeReaderValue<unsigned int>>(*treeReader_, "nJet");
  if (branchFound("Jet_pt"))            Jet_pt            = make_unique<TTreeReaderArray<float>>(*treeReader_, "Jet_pt");
  if (branchFound("Jet_pt_nom"))        Jet_nomPt         = make_unique<TTreeReaderArray<float>>(*treeReader_, "Jet_pt_nom");
  if (branchFound("Jet_eta"))           Jet_eta           = make_unique<TTreeReaderArray<float>>(*treeReader_, "Jet_eta");
  if (branchFound("Jet_phi"))           Jet_phi           = make_unique<TTreeReaderArray<float>>(*treeReader_, "Jet_phi");
  if (branchFound("Jet_nConstituents")) Jet_nConstituents = make_unique<TTreeReaderArray<int>>(*treeReader_, "Jet_nConstituents");
  if (branchFound("Jet_nMuons"))        Jet_nMuons        = make_unique<TTreeReaderArray<int>>(*treeReader_, "Jet_nMuons");
  if (branchFound("Jet_nElectrons"))    Jet_nElectrons    = make_unique<TTreeReaderArray<int>>(*treeReader_, "Jet_nElectrons");
  if (branchFound("Jet_puId"))          Jet_puId          = make_unique<TTreeReaderArray<int>>(*treeReader_, "Jet_puId");
  if (branchFound("Jet_muonIdx1"))      Jet_muIdx1        = make_unique<TTreeReaderArray<int>>(*treeReader_, "Jet_muonIdx1");
  if (branchFound("Jet_muonIdx2"))      Jet_muIdx2        = make_unique<TTreeReaderArray<int>>(*treeReader_, "Jet_muonIdx2");
  if (branchFound("Jet_electronIdx1"))  Jet_elIdx1        = make_unique<TTreeReaderArray<int>>(*treeReader_, "Jet_electronIdx1");
  if (branchFound("Jet_electronIdx2"))  Jet_elIdx2        = make_unique<TTreeReaderArray<int>>(*treeReader_, "Jet_electronIdx2");
  if (branchFound("Jet_jetId"))         Jet_jetId         = make_unique<TTreeReaderArray<int>>(*treeReader_, "Jet_jetId");
  if (branchFound("Jet_qgl"))           Jet_qgl           = make_unique<TTreeReaderArray<float>>(*treeReader_, "Jet_qgl");
  if (branchFound("Jet_neHEF"))         Jet_neHEF         = make_unique<TTreeReaderArray<float>>(*treeReader_, "Jet_neHEF");
  if (branchFound("Jet_neEmEF"))        Jet_neEmEF        = make_unique<TTreeReaderArray<float>>(*treeReader_, "Jet_neEmEF");
  if (branchFound("Jet_chHEF"))         Jet_chHEF         = make_unique<TTreeReaderArray<float>>(*treeReader_, "Jet_chHEF");
  if (branchFound("Jet_chEmEF"))        Jet_chEmEF        = make_unique<TTreeReaderArray<float>>(*treeReader_, "Jet_chEmEF");
  if (branchFound("Jet_area"))          Jet_area          = make_unique<TTreeReaderArray<float>>(*treeReader_, "Jet_area");
  if (branchFound("Jet_mass"))          Jet_mass          = make_unique<TTreeReaderArray<float>>(*treeReader_, "Jet_mass");
  if (branchFound("Jet_mass_nom"))      Jet_nomMass       = make_unique<TTreeReaderArray<float>>(*treeReader_, "Jet_mass_nom");
  if (branchFound("Jet_btagCSVV2"))     Jet_btagCSVV2     = make_unique<TTreeReaderArray<float>>(*treeReader_, "Jet_btagCSVV2");
  if (branchFound("Jet_btagDeepFlavB")) Jet_btagDeepFlavB = make_unique<TTreeReaderArray<float>>(*treeReader_, "Jet_btagDeepFlavB");

  // FatJet
  if (branchFound("nFatJet"))                 nFatJet                 = make_unique<TTreeReaderValue<unsigned int>>(*treeReader_, "nFatJet");
  if (branchFound("FatJet_pt"))               FatJet_pt               = make_unique<TTreeReaderArray<float>>(*treeReader_, "FatJet_pt");
  if (branchFound("FatJet_eta"))              FatJet_eta              = make_unique<TTreeReaderArray<float>>(*treeReader_, "FatJet_eta");
  if (branchFound("FatJet_phi"))              FatJet_phi              = make_unique<TTreeReaderArray<float>>(*treeReader_, "FatJet_phi");
  if (branchFound("FatJet_mass"))             FatJet_mass             = make_unique<TTreeReaderArray<float>>(*treeReader_, "FatJet_mass");
  if (branchFound("FatJet_msoftdrop"))        FatJet_msoftdrop        = make_unique<TTreeReaderArray<float>>(*treeReader_, "FatJet_msoftdrop");
  if (branchFound("FatJet_jetId"))            FatJet_jetId            = make_unique<TTreeReaderArray<int>>(*treeReader_, "FatJet_jetId"); 
  if (branchFound("FatJet_btagDeepB"))        FatJet_btagDeepB        = make_unique<TTreeReaderArray<float>>(*treeReader_, "FatJet_btagDeepB");
  if (branchFound("FatJet_btagCSVV2"))        FatJet_btagCSVV2        = make_unique<TTreeReaderArray<float>>(*treeReader_, "FatJet_btagCSVV2");
  if (branchFound("FatJet_n2b1"))             FatJet_n2b1             = make_unique<TTreeReaderArray<float>>(*treeReader_, "FatJet_n2b1");
  if (branchFound("FatJet_n3b1"))             FatJet_n3b1             = make_unique<TTreeReaderArray<float>>(*treeReader_, "FatJet_n3b1");
  if (branchFound("FatJet_hadronFlavour"))    FatJet_hadronFlavour    = make_unique<TTreeReaderArray<int>>(*treeReader_, "FatJet_hadronFlavour");
  if (branchFound("FatJet_nBHadrons"))        FatJet_nBHadrons        = make_unique<TTreeReaderArray<unsigned char>>(*treeReader_, "FatJet_nBHadrons");
  if (branchFound("FatJet_nCHadrons"))        FatJet_nCHadrons        = make_unique<TTreeReaderArray<unsigned char>>(*treeReader_, "FatJet_nCHadrons");
  if (branchFound("FatJet_rawFactor"))        FatJet_rawFactor        = make_unique<TTreeReaderArray<float>>(*treeReader_, "FatJet_rawFactor");
  if (branchFound("FatJet_subJetIdx1"))       FatJet_subJetIdx1       = make_unique<TTreeReaderArray<int>>(*treeReader_, "FatJet_subJetIdx1");
  if (branchFound("FatJet_subJetIdx2"))       FatJet_subJetIdx2       = make_unique<TTreeReaderArray<int>>(*treeReader_, "FatJet_subJetIdx2");
  if (branchFound("FatJet_tau1"))             FatJet_tau1             = make_unique<TTreeReaderArray<float>>(*treeReader_, "FatJet_tau1");
  if (branchFound("FatJet_tau2"))             FatJet_tau2             = make_unique<TTreeReaderArray<float>>(*treeReader_, "FatJet_tau2");
  if (branchFound("FatJet_tau3"))             FatJet_tau3             = make_unique<TTreeReaderArray<float>>(*treeReader_, "FatJet_tau3");
  if (branchFound("FatJet_tau4"))             FatJet_tau4             = make_unique<TTreeReaderArray<float>>(*treeReader_, "FatJet_tau4");
  if (branchFound("FatJet_deepTag_WvsQCD"))   FatJet_deepTag_WvsQCD   = make_unique<TTreeReaderArray<float>>(*treeReader_, "FatJet_deepTag_WvsQCD");
  if (branchFound("FatJet_deepTag_ZvsQCD"))   FatJet_deepTag_ZvsQCD   = make_unique<TTreeReaderArray<float>>(*treeReader_, "FatJet_deepTag_ZvsQCD");
  if (branchFound("FatJet_deepTag_TvsQCD"))   FatJet_deepTag_TvsQCD   = make_unique<TTreeReaderArray<float>>(*treeReader_, "FatJet_deepTag_TvsQCD");
  if (branchFound("FatJet_deepTagMD_WvsQCD")) FatJet_deepTagMD_WvsQCD = make_unique<TTreeReaderArray<float>>(*treeReader_, "FatJet_deepTagMD_WvsQCD");
  if (branchFound("FatJet_deepTagMD_ZvsQCD")) FatJet_deepTagMD_ZvsQCD = make_unique<TTreeReaderArray<float>>(*treeReader_, "FatJet_deepTagMD_ZvsQCD");
  if (branchFound("FatJet_deepTagMD_TvsQCD")) FatJet_deepTagMD_TvsQCD = make_unique<TTreeReaderArray<float>>(*treeReader_, "FatJet_deepTagMD_TvsQCD");

  if (branchFound("FatJet_electronIdx3SJ"))   FatJet_electronIdx3SJ   = make_unique<TTreeReaderArray<int>>(*treeReader_, "FatJet_electronIdx3SJ");
  if (branchFound("FatJet_muonIdx3SJ"))       FatJet_muonIdx3SJ       = make_unique<TTreeReaderArray<int>>(*treeReader_, "FatJet_muonIdx3SJ");

  // SubJet
  if (branchFound("nSubJet")) nSubJet                     = make_unique<TTreeReaderValue<unsigned int>>(*treeReader_, "nSubJet");
  if (branchFound("SubJet_pt")) SubJet_pt                 = make_unique<TTreeReaderArray<float>>(*treeReader_, "SubJet_pt");
  if (branchFound("SubJet_eta")) SubJet_eta               = make_unique<TTreeReaderArray<float>>(*treeReader_, "SubJet_eta");
  if (branchFound("SubJet_phi")) SubJet_phi               = make_unique<TTreeReaderArray<float>>(*treeReader_, "SubJet_phi");
  if (branchFound("SubJet_mass")) SubJet_mass             = make_unique<TTreeReaderArray<float>>(*treeReader_, "SubJet_mass");
  if (branchFound("SubJet_btagDeepB")) SubJet_btagDeepB   = make_unique<TTreeReaderArray<float>>(*treeReader_, "SubJet_btagDeepB");
  if (branchFound("SubJet_rawFactor")) SubJet_rawFactor   = make_unique<TTreeReaderArray<float>>(*treeReader_, "SubJet_rawFactor");

  // Tau
  if (branchFound("nTau"))            nTau            = make_unique<TTreeReaderValue<unsigned int>>(*treeReader_, "nTau");
  if (branchFound("Tau_pt"))          Tau_pt          = make_unique<TTreeReaderArray<float>>(*treeReader_, "Tau_pt");
  if (branchFound("Tau_eta"))         Tau_eta         = make_unique<TTreeReaderArray<float>>(*treeReader_, "Tau_eta");
  if (branchFound("Tau_phi"))         Tau_phi         = make_unique<TTreeReaderArray<float>>(*treeReader_, "Tau_phi");
  if (branchFound("Tau_mass"))        Tau_mass        = make_unique<TTreeReaderArray<float>>(*treeReader_, "Tau_mass");
  if (branchFound("Tau_dxy"))         Tau_dxy         = make_unique<TTreeReaderArray<float>>(*treeReader_, "Tau_dxy");
  if (branchFound("Tau_dz"))          Tau_dz          = make_unique<TTreeReaderArray<float>>(*treeReader_, "Tau_dz"); 
  if (branchFound("Tau_charge"))      Tau_charge      = make_unique<TTreeReaderArray<int>>(*treeReader_, "Tau_charge");
  if (branchFound("Tau_decayMode"))   Tau_decayMode   = make_unique<TTreeReaderArray<int>>(*treeReader_, "Tau_decayMode");
  if (branchFound("Tau_idDecayMode")) Tau_idDecayMode = make_unique<TTreeReaderArray<bool>>(*treeReader_, "Tau_idDecayMode");
  if (branchFound("Tau_idDecayModeNewDMs")) 
                   Tau_idDecayModeNewDMs = make_unique<TTreeReaderArray<bool>>(*treeReader_, "Tau_idDecayModeNewDMs");
  if (branchFound("Tau_jetIdx"))      Tau_jetIdx      = make_unique<TTreeReaderArray<int>>(*treeReader_, "Tau_jetIdx");
  if (branchFound("Tau_idAntiEle"))   Tau_idAntiEle   = make_unique<TTreeReaderArray<unsigned char>>(*treeReader_, "Tau_idAntiEle");
  if (branchFound("Tau_idAntiMu"))    Tau_idAntiMu    = make_unique<TTreeReaderArray<unsigned char>>(*treeReader_, "Tau_idAntiMu");
  if (branchFound("Tau_idMVAoldDM"))  Tau_idMVAoldDM  = make_unique<TTreeReaderArray<unsigned char>>(*treeReader_, "Tau_idMVAoldDM");
  if (branchFound("Tau_idDeepTau2017v2VSjet")) 
                   Tau_idDeepTau2017v2VSjet = make_unique<TTreeReaderArray<unsigned char>>(*treeReader_, "Tau_idDeepTau2017v2VSjet");
  if (branchFound("Tau_idDeepTau2017v2VSmu")) 
                   Tau_idDeepTau2017v2VSmu = make_unique<TTreeReaderArray<unsigned char>>(*treeReader_, "Tau_idDeepTau2017v2VSmu");
  if (branchFound("Tau_idDeepTau2017v2VSe")) 
                   Tau_idDeepTau2017v2VSe = make_unique<TTreeReaderArray<unsigned char>>(*treeReader_, "Tau_idDeepTau2017v2VSe");
  if (branchFound("Tau_idDeepTau2017v2p1VSjet")) 
                   Tau_idDeepTau2017v2p1VSjet = make_unique<TTreeReaderArray<unsigned char>>(*treeReader_, "Tau_idDeepTau2017v2p1VSjet");
  if (branchFound("Tau_idDeepTau2017v2p1VSmu")) 
                   Tau_idDeepTau2017v2p1VSmu = make_unique<TTreeReaderArray<unsigned char>>(*treeReader_, "Tau_idDeepTau2017v2p1VSmu");
  if (branchFound("Tau_idDeepTau2017v2p1VSe")) 
                   Tau_idDeepTau2017v2p1VSe = make_unique<TTreeReaderArray<unsigned char>>(*treeReader_, "Tau_idDeepTau2017v2p1VSe");

  if (isMC_) {
    // LHE
    if (branchFound("LHE_Njets")) LHEnJets  = make_unique<TTreeReaderValue<unsigned char>>(*treeReader_, "LHE_Njets");
    if (readGenInfo_) {
      if (branchFound("nLHEPart"))      nLHEPart      = make_unique<TTreeReaderValue<unsigned int>>(*treeReader_, "nLHEPart");
      if (branchFound("LHEPart_pt"))    LHEPart_pt    = make_unique<TTreeReaderArray<float>>(*treeReader_, "LHEPart_pt");
      if (branchFound("LHEPart_eta"))   LHEPart_eta   = make_unique<TTreeReaderArray<float>>(*treeReader_, "LHEPart_eta");
      if (branchFound("LHEPart_phi"))   LHEPart_phi   = make_unique<TTreeReaderArray<float>>(*treeReader_, "LHEPart_phi");
      if (branchFound("LHEPart_mass"))  LHEPart_mass  = make_unique<TTreeReaderArray<float>>(*treeReader_, "LHEPart_mass");      
      if (branchFound("LHEPart_pdgId")) LHEPart_pdgId = make_unique<TTreeReaderArray<int>>(*treeReader_, "LHEPart_pdgId");

      // GenParticle
      if (branchFound("nGenPart"))                 nGenPart            = make_unique<TTreeReaderValue<unsigned int>>(*treeReader_, "nGenPart");
      if (branchFound("GenPart_pt"))               GenPart_pt          = make_unique<TTreeReaderArray<float>>(*treeReader_, "GenPart_pt");
      if (branchFound("GenPart_eta"))              GenPart_eta         = make_unique<TTreeReaderArray<float>>(*treeReader_, "GenPart_eta");
      if (branchFound("GenPart_phi"))              GenPart_phi         = make_unique<TTreeReaderArray<float>>(*treeReader_, "GenPart_phi");
      if (branchFound("GenPart_mass"))             GenPart_mass        = make_unique<TTreeReaderArray<float>>(*treeReader_, "GenPart_mass");
      if (branchFound("GenPart_genPartIdxMother")) GenPart_motherIdx   = make_unique<TTreeReaderArray<int>>(*treeReader_, "GenPart_genPartIdxMother");
      if (branchFound("GenPart_pdgId"))            GenPart_pdgId       = make_unique<TTreeReaderArray<int>>(*treeReader_, "GenPart_pdgId");
      if (branchFound("GenPart_status"))           GenPart_status      = make_unique<TTreeReaderArray<int>>(*treeReader_, "GenPart_status");
      if (branchFound("GenPart_statusFlags"))      GenPart_statusFlags = make_unique<TTreeReaderArray<int>>(*treeReader_, "GenPart_statusFlags");
    }
  }
  // Setting the HLT branch pointers
  setHltPtrList(doubleMuonHltPathList_, doubleMuonHltPtrList_);
  setHltPtrList(singleMuonHltPathList_, singleMuonHltPtrList_);
  setHltPtrList(doubleEgHltPathList_, doubleEgHltPtrList_);
  setHltPtrList(singleElectronHltPathList_, singleElectronHltPtrList_);
  setHltPtrList(muonEgHltPathList_, muonEgHltPtrList_);
  setHltPtrList(doubleMuonHltPathList_, doubleMuonHltPtrList_);

  setHltPtrList(singleMuonHltForFakePathList_, singleMuonHltForFakePtrList_);
  setHltPtrList(singleElectronHltForFakePathList_, singleElectronHltForFakePtrList_);

  return true;
}
void AnaBase::setHltPtrList(const vector<string>& hltPathList, std::vector<std::unique_ptr<TTreeReaderValue<bool>>>& hltPtrList) {
  for (const auto& hlt: hltPathList) {
    if (branchFound(hlt)) 
      hltPtrList.push_back(std::move(make_unique<TTreeReaderValue<bool>>(*treeReader_, hlt.c_str())));
  }
}
vector<bool> AnaBase::getHLTscores(const vector<std::unique_ptr<TTreeReaderValue<bool>>>& ptrs) {
  vector<bool> hltScores;
  for (const auto& hlt: ptrs)
    hltScores.push_back(*hlt->Get());

  return hltScores;
}
bool AnaBase::isDuplicate(bool passDoubleMuonHLT, bool passDoubleEgHLT, bool passMuonEgHLT,
			  bool passSingleMuonHLT, bool passSingleEleHLT, const string& dataset) {
  map<string, bool> options {
    {"DoubleMuon", passDoubleMuonHLT},
    {"DoubleEG", passDoubleEgHLT && !passDoubleMuonHLT},
    {"MuonEG", passMuonEgHLT && !(passDoubleMuonHLT || passDoubleEgHLT)},
    {"SingleMuon", passSingleMuonHLT && !(passDoubleMuonHLT || passDoubleEgHLT || passMuonEgHLT)},
    {"SingleElectron", passSingleEleHLT && !(passDoubleMuonHLT || passDoubleEgHLT || passMuonEgHLT || passSingleMuonHLT)}
  };
  map<string, bool>::iterator it = options.find(dataset.c_str());
  if (it != options.end() && it->second) return true;

  return false;
}
// if an event has paased HLT condition or not
bool AnaBase::isTriggered(const vector<string>& paths, const vector<bool>& scores) {
  for (size_t i = 0; i < paths.size(); ++i)
    if (scores[i]) return true;

  return false;
}
void AnaBase::endJob() 
{
}
void AnaBase::closeHistFile()
{ 
  histf_->cd();
  histf_->Write();
  histf_->Close();
}
double AnaBase::lumiWt(double evtWeightSum, bool verbose) const 
{
  double nevt = (evtWeightSum > -1) ? evtWeightSum : AnaUtil::cutValue(lumiWtMap(), "nevents");
  if (verbose) 
    cout << "==> INFO. intLumi: " << AnaUtil::cutValue(lumiWtMap(), "intLumi") 
	 << " xsec: " << AnaUtil::cutValue(lumiWtMap(), "xsec") 
	 << " nevt: " << nevt 
	 << endl; 
  return (AnaUtil::cutValue(lumiWtMap(), "intLumi") * AnaUtil::cutValue(lumiWtMap(), "xsec") / nevt); 
}
// ---------------------------------
// Add input Root files to the chain
// ---------------------------------
int AnaBase::setInputFile(const string& fname) //Called by readJob
{
  auto found = fname.find("root:");
  if (found == string::npos && gSystem->AccessPathName(fname.c_str())) {
    cerr << "==> WARN: File <<" << fname << ">> was not found!!" << endl;
    return static_cast<int>(chain_->GetEntries()); 
  }
  chain_->AddFile(fname.c_str(), -1);
  chainRun_->AddFile(fname.c_str(), -1);

  return static_cast<int>(chain_->GetEntries()); 
}
// ---------------------------------------
// Get total number of events in the chain
// --------------------------------------
int AnaBase::getEntries() const 
{
  return static_cast<int>(chain_->GetEntriesFast());
}
// ------------------------------------------
// Open output files with global filehandles
// ------------------------------------------
bool AnaBase::openFiles()
{
  fLog_.open(logFile_.c_str(), ios::out);
  if (!fLog_) {
    cerr << "==> ERROR. File: " << logFile_ << " could not be opened!" << endl;
    return false;
  }
  fLog_ << setiosflags(ios::fixed);

  evLog_.open(evFile_.c_str(), ios::out);
  if (!evLog_) {
    cerr << "==> ERROR. File: " << evFile_ << " could not be opened!" << endl;
    return false;
  }
  evLog_ << setiosflags(ios::fixed);

  selEvLog_.open(selEvFile_.c_str(), ios::out);
  if (!selEvLog_) {
    cerr << "==> ERROR. File: " << selEvFile_ << " could not be opened!" << endl;
    return false;
  }
  selEvLog_ << setiosflags(ios::fixed);

  return true;
}
// ------------------------
// Close output files
// ------------------------
void AnaBase::closeFiles() 
{
  if (fLog_) {
    fLog_ << resetiosflags(ios::fixed); 
    fLog_.close();
  }
  if (evLog_) {
    evLog_ << resetiosflags(ios::fixed); 
    evLog_.close();
  }
  if (selEvLog_) {
    selEvLog_ << resetiosflags(ios::fixed); 
    selEvLog_.close();
  }
  closeHistFile();
}
bool AnaBase::branchFound(const string& b, const string& tname)
{
  TBranch* branch = (tname == "Events") ? chain_->GetBranch(b.c_str())
                                        : chainRun_->GetBranch(b.c_str());  // Get branch pointer
  if (branch == nullptr) {
    cout << ">>> Branch: <" << b << "> N O T   F O U N D !!!" << endl;
    return false;
  }
  cout << ">>> Branch: <" << b << "> found!" << endl;

  return true;
}
bool AnaBase::readJob(const string& jobFile, int& nFiles)
{
  // Open the file containing the datacards
  std::ifstream fin(jobFile.c_str(), ios::in);    
  if (!fin) {
    cerr << "==> ERROR. Input File: <<" << jobFile << ">> could not be opened!" << endl;
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
        if (isMC_ && dt.size() > 1) 
          isSignal_ = (dt.at(1) == "SIGNAL") ? true : false;       
      }
    }
    else if (key == "era")
      era_ = stoi(value.c_str());    
    else if (key == "evtWtSum")
      evtWtSum_ = value;
    else if (key == "dataset")
      dataset_ = value;
    else if (key == "readGenInfo")
      readGenInfo_ = stoi(value.c_str()) > 0 ? true : false;
    else if (key == "logFile")
      logFile_ = value;
    else if (key == "eventFile")
      evFile_  = value;
    else if (key == "selEventFile")
      selEvFile_  = value;
    else if (key == "logOption")
      logOption_ = strtol(value.c_str(), NULL, 2);
    else if (key == "maxEvent")
      maxEvt_ = stoi(value.c_str());
    else if (key == "nFiles")
      nFiles_ = stoi(value.c_str());
    else if (key == "startEvent")
      firstEvt_ = stoi(value.c_str());
    else if (key == "endEvent")
      lastEvt_ = stoi(value.c_str());
    else if (key == "histFile")
      histFile_ = value;
    //else if (key == "fakehistFile")
    //  fakehistFile_ = value;
    else if (key == "inputFile")
      AnaUtil::buildList(tokens, fileList_);
    else if (key == "muonIdSFRootFile")
      SFHandler_.muonIdSFRootFile_ = value;
    else if (key == "looseMuonIdSFhistName")
      SFHandler_.looseMuonIdSFhistName_ = value;
    else if (key == "medMuonIdSFhistName")
      SFHandler_.medMuonIdSFhistName_ = value;
    else if (key == "tightMuonIdSFhistName")
      SFHandler_.tightMuonIdSFhistName_ = value;
    else if (key == "electronLooseIdSFRootFile")
      SFHandler_.electronLooseIdSFRootFile_ = value;
    else if (key == "looseEleIdSFhistName")
      SFHandler_.looseEleIdSFhistName_ = value;
    else if (key == "electronTightIdSFRootFile")
      SFHandler_.electronTightIdSFRootFile_ = value;
    else if (key == "tightEleIdSFhistName")
      SFHandler_.tightEleIdSFhistName_ = value;
    else if (key == "muonTightIsoSFRootFile")
      SFHandler_.muonTightIsoSFRootFile_ = value;
    else if (key == "tightMuIsoSFhistName")
      SFHandler_.tightMuIsoSFhistName_ = value;
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
    else if (key == "SingleMuonForFake")
      AnaUtil::buildList(tokens, singleMuonHltForFakePathList_);
    else if (key == "SingleElectronForFake")
      AnaUtil::buildList(tokens, singleElectronHltForFakePathList_);
    else if (key == "eventId" && tokens.size() == 4)
      AnaUtil::buildMap(tokens, eventIdMap_);
    else {
      AnaUtil::storeCuts(tokens, hmap_);
    }
  }
  // Close the file
  fin.close();

  if (!isSignal_) readGenInfo_ = false;

  // Build the chain of root files
  for (const auto& fname: fileList_) {
    cout << "==> INFO. Adding input file " << fname << " to TChain " << endl;
    ++nFiles;
    int nevt = setInputFile(fname);
    std::cout<<">>>Right now the chain contains :: "<<nevt<<" events\n";
    if (nFiles_ > 0 && nFiles > nFiles_) break;
    // if (maxEvt_ > 0 && nevt >= maxEvt_) break; // -- Useful for debug -- //
  }
  if (!nFiles) {
    cerr << "==> WARN. Input Root file list is empty! exiting ..." << endl;
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
     << "      logOption: " << logOption_ << endl
     << "       maxEvent: " << maxEvt_ 
     << endl;

  // InputFiles
  if (chain_ != nullptr) {
    TObjArray *fileElements = chain_->GetListOfFiles();
    os << "==> INFO. nFiles: " << fileElements->GetEntries() 
       << ", Files to analyse:" 
       << endl;
    TIter next(fileElements);
    TChainElement* chEl = nullptr;
    while (( chEl = dynamic_cast<TChainElement*>(next()) ))
      os << chEl->GetTitle() 
         << endl;
  }
  else
    AnaUtil::showList(fileList_, "==> INFO. inputFiles:", os);

  // EventID 
  AnaUtil::showMap<string, int>(eventIdMap_, "Event List:", os);

  // Cuts
  AnaUtil::showCuts(hmap_, os);
}
