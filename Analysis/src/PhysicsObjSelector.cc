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
#include <sstream>
#include <utility> 
#include <typeinfo>

#include "TFile.h"
#include "TH1K.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TVector2.h"
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include "AnaUtil.h"
#include "PhysicsObjSelector.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::ostringstream;
using std::vector;
using std::map;
using std::pair;
using std::setprecision;
using std::setw;

PhysicsObjSelector::PhysicsObjSelector()
  : AnaBase()
{}

bool PhysicsObjSelector::beginJob() {
  if (!AnaBase::beginJob()) return false;

  // Set Branch Addresses :: VERY VERY IMPORTANT
  if (!AnaBase::init()) {
    std::cout<<"ERROR: Branches not Found!!!"<<std::endl;
  }

  histf()->cd();
  histf()->mkdir("ObjectSelection");
    
  return true;
}

void PhysicsObjSelector::endJob() {
  objectEfficiency();
  AnaBase::endJob();
}
bool PhysicsObjSelector::readJob(const string& jobFile, int& nFiles)
{
  if (!AnaBase::readJob(jobFile, nFiles)) return false;
  return true;
}
void PhysicsObjSelector::objectEfficiency() {
  // Object selection Efficiency
  histf()->cd();
  histf()->cd("ObjectSelection");

  // $  : objects are saved at these levels
  // -- : to show the cutflow
  // Muon
  vector<string> muLabels {
    "nRawMuons",
    " --$afterPreSel cuts",
    "   --pT>10 & jet cleaning [NA]",
    "     --$isFakeable",
    "       --$isIsoFakeable",
    "     --$isTight",
    "       --$isTightIso"  
  };
  AnaUtil::showEfficiency("muCutFlow", muLabels, "Muon Selection", "Muons");  

  // Electron
  vector<string> eleLabels {
    "nRawElectrons",
    " --preSelElectrons",
    "   --$mucleaning",
    " --pT>10,ConvVeto,NoLostHits",
    "   --$fakeable",
    "   --tight",
    "     --$mCleaning"
  };
  AnaUtil::showEfficiency("eleCutFlow", eleLabels, "Electron Selection", "Electrons");  

  // Jet
  vector<string> jetLabels {
    "nRawJets",
    " --looseSelection",
    "   --$puId cut",
    "     --$lepton cleaning"
  };
  AnaUtil::showEfficiency("jetCutFlow", jetLabels, "Jet Selection", "Jets");  

  // Tau
  vector<string> tauLabels {
    "nRawTaus",
    " --$tauSelCuts",
    "   --$lepClean"
  };
  AnaUtil::showEfficiency("tauCutFlow", tauLabels, "Tau Selection", "Taus");  

  // FatJet
  vector<string> fatJetLabels {
    "nRawFatJets",
    " --after jetId, pT & eta cuts",
    "   --hasValidSubJets",
    "     --30 < mSoftDrop < 210",
    "       --tau21 < 0.75",
    "         --leptonCleaning",
    "           --isBTagged"  
  };
  AnaUtil::showEfficiency("fatJetCutFlow", fatJetLabels, "FatJet Selection", "FatJets");  


}
void PhysicsObjSelector::bookHistograms() {
  histf()->cd();
  histf()->cd("ObjectSelection");

  new TH1D("muCutFlow", "Muon CutFlow", 7, -0.5, 6.5);
  new TH1D("eleCutFlow", "Electron CutFlow", 7, -0.5, 6.5);
  new TH1D("jetCutFlow", "Jet CutFlow", 4, -0.5, 3.5);
  new TH1D("tauCutFlow", "Tau CutFlow", 3, -0.5, 2.5);
  new TH1D("fatJetCutFlow", "FatJet CutFlow", 7, -0.5, 6.5);
}

// clear lists
void PhysicsObjSelector::clear() {
  eventList_.clear();
  metList_.clear();

  preSelMuList_.clear();
  fakeableMuList_.clear();
  fakeableIsoMuList_.clear();
  tightMuList_.clear();
  tightIsoMuList_.clear();

  preSelEleList_.clear();
  fakeableIsoEleList_.clear();
  tightIsoEleList_.clear();

  preSelJetList_.clear();
  leptonCleanJetList_.clear();
  looseBJetList_.clear();
  bJetList_.clear();

  tauList_.clear();
  leptonCleanTauList_.clear();

  fatJetList_.clear();
  cleanFatJetList_.clear();
  bTaggedFatJetList_.clear();

  subJetList_.clear();

  genParticleList_.clear();
  lheParticleList_.clear();

  searchedMu_  = false; 
  searchedEle_ = false; 
  searchedJet_ = false;
  searchedFatJet_ = false;
}

double PhysicsObjSelector::getGenSumW(){
  return *genEventSumw_ -> Get();
}

bool PhysicsObjSelector::findEventInfo() {
  vhtm::Event ev;
  ev.run    = *run_->Get();
  ev.event  = *event_->Get();
  ev.lumis  = *lumis_->Get();
  if (isMC()){
    ev.genEvWt = *genEvWt->Get();
    ev.PUWeight = *PU_Weight->Get();
    if (readGenInfo()) {
      ev.nGenParticles = *nGenPart->Get();
      ev.nLHEParticles = *nLHEPart->Get();
    }
    //ev.nLHEJets = *LHEnJets->Get();
    ev.btagWeight_CSVV2 = *btagWeight_CSVV2->Get();
  }
  ev.nPV     = *nPV->Get();
  ev.nGoodPV = *nGoodPV->Get();
  ev.PVscore = *PV_score->Get();
  ev.PVchi2  = *PV_chi2->Get();
  ev.PVndf   = *PV_ndf->Get();
  ev.PVx     = *PV_xPos->Get();
  ev.PVy     = *PV_yPos->Get();
  ev.PVz     = *PV_zPos->Get();

  eventList_.push_back(ev);

  return true;
}

std::vector < bool > PhysicsObjSelector::getDoubleMuonHLTscores() {
  std::vector < TTreeReaderValue<bool>* > ptrs   = AnaBase::getDoubleMuonHLTptrs();
  std::vector < bool > hltScores;
  if (ptrs.size() > 0)
    for (TTreeReaderValue<bool>* hlt : ptrs) {
      hltScores.push_back(*hlt->Get());
    }
  return hltScores;
}
std::vector < bool > PhysicsObjSelector::getSingleMuonHLTscores() {
  std::vector < TTreeReaderValue<bool>* > ptrs   = AnaBase::getSingleMuonHLTptrs();
  std::vector < bool > hltScores;
  if (ptrs.size() > 0)
    for (TTreeReaderValue<bool>* hlt : ptrs) {
      hltScores.push_back(*hlt->Get());
    }
  return hltScores;
}
std::vector < bool > PhysicsObjSelector::getDoubleEgHLTscores() {
  std::vector < TTreeReaderValue<bool>* > ptrs   = AnaBase::getDoubleEgHLTptrs();
  std::vector < bool > hltScores;
  if (ptrs.size() > 0)
    for (TTreeReaderValue<bool>* hlt : ptrs) {
      hltScores.push_back(*hlt->Get());
    }
  return hltScores;
}
std::vector < bool > PhysicsObjSelector::getSingleElectronHLTscores() {
  std::vector < TTreeReaderValue<bool>* > ptrs   = AnaBase::getSingleElectronHLTptrs();
  std::vector < bool > hltScores;
  if (ptrs.size() > 0)
    for (TTreeReaderValue<bool>* hlt : ptrs) {
      hltScores.push_back(*hlt->Get());
    }
  return hltScores;
}
std::vector < bool > PhysicsObjSelector::getMuonEgHLTscores() {
  std::vector < TTreeReaderValue<bool>* > ptrs   = AnaBase::getMuonEgHLTptrs();
  std::vector < bool > hltScores;
  if (ptrs.size() > 0)
    for (TTreeReaderValue<bool>* hlt : ptrs) {
      hltScores.push_back(*hlt->Get());
    }
  return hltScores;
}

bool PhysicsObjSelector::findGenPartInfo() {
  histf()->cd();
  histf()->cd("ObjectSelection");

  if (isMC() && readGenInfo()) {
    for (size_t j = 0; j < (*nGenPart->Get()); ++j){      
      vhtm::GenParticle gp;
      
      gp.index = j;
      gp.pt   = GenPart_pt->At(j);
      gp.eta  = GenPart_eta->At(j);
      gp.phi  = GenPart_phi->At(j);
      gp.mass = GenPart_mass->At(j);
      gp.pdgId= GenPart_pdgId->At(j);
      gp.motherIdx = GenPart_motherIdx->At(j);
      gp.status = GenPart_status->At(j);
      gp.statusFlags = GenPart_statusFlags->At(j);
    
      genParticleList_.push_back(gp);
    }
  }
  return true;
}

bool PhysicsObjSelector::findLHEPartInfo() {
  histf()->cd();
  histf()->cd("ObjectSelection");

  if (isMC() && readGenInfo()){
    for (size_t j = 0; j < (*nLHEPart->Get()); ++j){
      vhtm::LHEParticle gp;
      
      gp.index = j;
      gp.pt   = LHEPart_pt->At(j);
      gp.eta  = LHEPart_eta->At(j);
      gp.phi  = LHEPart_phi->At(j);
      gp.mass = LHEPart_mass->At(j);
      gp.pdgId= LHEPart_pdgId->At(j);
      
      lheParticleList_.push_back(gp);
    }
  }
  return true;
}


void PhysicsObjSelector::findObjects() {
  muonSelector();  
  electronSelector();
  jetSelector();
  fatJetSelector();
  subJetSelector();
  metSelector();
  tauSelector();
}

// Muon selection
void PhysicsObjSelector::muonSelector() {
  histf()->cd();
  histf()->cd("ObjectSelection");
  for (size_t i = 0; i < (*nMuon->Get()); ++i){
    AnaUtil::fillHist1D ("muCutFlow", 0, 1.0);
    // Loose muon selection
    if (Muon_corrpt->At(i) < 5.0) continue;
    if (std::fabs(Muon_eta->At(i)) > 2.4) continue;
    if (std::fabs(Muon_dxy->At(i)) > 0.05) continue;
    if (std::fabs(Muon_dz->At(i))  > 0.1) continue;
    if (Muon_miniPFRelIso_all->At(i) > 0.4) continue;
    if (std::fabs(Muon_sip3d->At(i)) > 8) continue;
    if (!Muon_LooseId->At(i)) continue;

    AnaUtil::fillHist1D ("muCutFlow", 1, 1.0);

    vhtm::Muon mu;
    mu.index   = i;
    mu.pt      = Muon_corrpt->At(i);
    mu.eta     = Muon_eta->At(i);
    mu.phi     = Muon_phi->At(i);
    mu.mass    = Muon_mass->At(i);
    mu.charge  = Muon_charge->At(i);
    mu.jetIdx  = Muon_jetIdx->At(i);
    if (isMC()){
      mu.genIdx = Muon_genPartIdx->At(i);
      mu.genFlv = Muon_genPartFlv->At(i); // genFlv = 1  : promt muon, genFlv = 15 : from tau decay
    }
    mu.mediumId  = Muon_MediumId->At(i);
    mu.tightId   = Muon_TightId->At(i);
    mu.highPtId  = Muon_highPtId->At(i);
    mu.pfRelIso03_all = Muon_pfRelIso03_all->At(i);
    mu.pfRelIso04_all = Muon_pfRelIso04_all->At(i);

    preSelMuList_.push_back(mu);

    // Fakeable Muon Selection
    if (Muon_corrpt->At(i) < 10.0) continue;
    //if (Muon_jetIdx->At(i) != -1) continue; // cleaning against jets
    AnaUtil::fillHist1D ("muCutFlow", 2, 1.0);

    if (!Muon_TightId->At(i)) {
      AnaUtil::fillHist1D ("muCutFlow", 3, 1.0);
      fakeableMuList_.push_back(mu);
      if (Muon_pfRelIso04_all->At(i) <= 0.15) {
	AnaUtil::fillHist1D ("muCutFlow", 4, 1.0);
	fakeableIsoMuList_.push_back(mu);
      }
    }
    
    // Tight Muon selection
    else if (Muon_TightId->At(i)) {
      AnaUtil::fillHist1D ("muCutFlow", 5, 1.0);
      tightMuList_.push_back(mu);
      if (Muon_pfRelIso04_all->At(i) <= 0.15) {
	AnaUtil::fillHist1D ("muCutFlow", 6, 1.0);
	tightIsoMuList_.push_back(mu);
      }
    }
  }
  searchedMu_ = true;
}


// Electron selecton
void PhysicsObjSelector::electronSelector() {
  if (!searchedMu_) std::cout<<">>>Muons are not selected yet!!!\n";
  histf()->cd();
  histf()->cd("ObjectSelection");
  for (size_t i = 0; i < (*nElectron->Get()); ++i){
    AnaUtil::fillHist1D ("eleCutFlow", 0, 1.0);

    if (Electron_pt->At(i) < 7) continue;
    if (std::fabs(Electron_eta->At(i)) > 2.5) continue;
    if (std::fabs(Electron_dxy->At(i)) > 0.05) continue;
    if (std::fabs(Electron_dz->At(i)) > 0.1) continue;
    if (Electron_sip3d->At(i) > 8) continue;
    if (!Electron_mvaFall17V2noIso_WPL->At(i)) continue;
    if (Electron_lostHits->At(i) > 1) continue;

    AnaUtil::fillHist1D ("eleCutFlow", 1, 1.0);

    vhtm::Electron el;
    el.index   = i;
    el.pt      = Electron_pt->At(i);
    el.eta     = Electron_eta->At(i);
    el.phi     = Electron_phi->At(i);
    el.mass    = Electron_mass->At(i);
    el.charge  = Electron_charge->At(i);
    el.jetIdx  = Electron_jetIdx->At(i);
    el.phoIdx  = Electron_phoIdx->At(i);
    if (isMC()){
      el.genIdx = Electron_genPartIdx->At(i);
      el.genFlv = Electron_genPartFlv->At(i);
    }

    if (!thisElectronIsMuon(el, true, false)) {
      preSelEleList_.push_back(el);
      AnaUtil::fillHist1D ("eleCutFlow", 2, 1.0); 
    }
    if (Electron_pt->At(i) < 10) continue;
    if (!Electron_convVeto->At(i)) continue;;
    if (Electron_lostHits->At(i) != 0) continue;
    AnaUtil::fillHist1D ("eleCutFlow", 3, 1.0); 

    if (!Electron_mvaFall17V2Iso_WP90->At(i)) {
      AnaUtil::fillHist1D ("eleCutFlow", 4, 1.0);
      fakeableIsoEleList_.push_back(el);
    }
    if (Electron_mvaFall17V2Iso_WP90->At(i)) {
      AnaUtil::fillHist1D ("eleCutFlow", 5, 1.0);
      if (!thisElectronIsMuon(el, false, true)) {
	tightIsoEleList_.push_back(el);
	AnaUtil::fillHist1D ("eleCutFlow", 6, 1.0);
      }
    }
  }
  searchedEle_ = true;
}

// Jet selection
void PhysicsObjSelector::jetSelector() {
  if (!(searchedMu_ && searchedEle_)) std::cout<<">>>Muon and Electron are not selected yet!!!\n";
  histf()->cd();
  histf()->cd("ObjectSelection");
  for (size_t i = 0; i < (*nJet->Get()); ++i){
    AnaUtil::fillHist1D ("jetCutFlow", 0, 1.0);
    if (!(Jet_jetId->At(i) & 2 
	  && Jet_nomPt->At(i) >= 25
	  && std::fabs(Jet_eta->At(i)) <= 2.5)
	) continue;
    AnaUtil::fillHist1D ("jetCutFlow", 1, 1.0);

    //Apply PuID on Jets with pT < 50 GeV (Need to be discussed!!!)
    if (Jet_nomPt->At(i) <= 50 && !((Jet_puId->At(i) >> 2) & 1))  continue;
    AnaUtil::fillHist1D ("jetCutFlow", 2, 1.0);

    vhtm::Jet jet;
    jet.index   = i;
    jet.pt      = Jet_nomPt->At(i);
    jet.eta     = Jet_eta->At(i);
    jet.phi     = Jet_phi->At(i);
    jet.mass    = Jet_mass->At(i);
    jet.muIdx1  = Jet_muIdx1->At(i);
    jet.muIdx2  = Jet_muIdx2->At(i);
    jet.elIdx1  = Jet_elIdx1->At(i);
    jet.elIdx2  = Jet_elIdx2->At(i);
    jet.btagDeepFlavB   = Jet_btagDeepFlavB->At(i);
    jet.btagCSVV2  = Jet_btagCSVV2->At(i);

    preSelJetList_.push_back(jet);

    if (!jetLeptonCleaning(jet)) continue;
    AnaUtil::fillHist1D ("jetCutFlow", 3, 1.0);
    leptonCleanJetList_.push_back(jet);

    //making b_jet collection
    if (Jet_btagDeepFlavB->At(i) >= 0.0521) looseBJetList_.push_back(jet);
    if (Jet_btagDeepFlavB->At(i) >= 0.3033) bJetList_.push_back(jet); 
  }
  searchedJet_ = true;
}


// MET selection
void PhysicsObjSelector::metSelector() {
  histf()->cd();
  histf()->cd("ObjectSelection");

  vhtm::MET met;
  met.pt    = *Met_nomPt->Get();
  met.phi   = *Met_nomPhi->Get();
  met.signf = *Met_significance->Get();
  met.sumEt = *Met_sumEt->Get();

  metList_.push_back(met);
}

// Tau selection
void PhysicsObjSelector::tauSelector() {
  if (!searchedJet_) jetSelector();
  histf()->cd();
  histf()->cd("ObjectSelection");

  for (size_t it = 0; it < (*nTau->Get()); ++it){
    AnaUtil::fillHist1D ("tauCutFlow", 0, 1.0);

    if (Tau_pt->At(it) < 20) continue;
    if (std::fabs(Tau_eta->At(it)) < 2.3) continue;
    if (std::fabs(Tau_dxy->At(it)) > 0.1) continue;
    if (std::fabs(Tau_dz->At(it)) > 0.2) continue;
    if (!Tau_idDecayModeNewDMs->At(it)) continue;
    if (!(Tau_decayMode->At(it) == 0 || Tau_decayMode->At(it) == 1 || Tau_decayMode->At(it) == 2 || Tau_decayMode->At(it) == 10 || Tau_decayMode->At(it) == 11)) continue;
    if (!(Tau_idDeepTau2017v2VSjet->At(it) >> 4 & 0x1 && Tau_idDeepTau2017v2VSe->At(it) >> 0 & 0x1 && Tau_idDeepTau2017v2VSmu->At(it) >> 0 & 0x1)) continue;

    AnaUtil::fillHist1D ("tauCutFlow", 1, 1.0);
    
    vhtm::Tau ta;
    ta.index   = it;
    ta.pt      = Tau_pt->At(it);
    ta.eta     = Tau_eta->At(it);
    ta.phi     = Tau_phi->At(it);
    ta.charge  = Tau_charge->At(it);
    ta.mass    = Tau_mass->At(it);
    ta.jetIdx  = Tau_jetIdx->At(it);

    tauList_.push_back(ta);
    
    if (!tauLeptonCleaning(ta)) continue;
    AnaUtil::fillHist1D ("tauCutFlow", 2, 1.0);
    leptonCleanTauList_.push_back(ta);
  }
}


// FatJet selection
void PhysicsObjSelector::fatJetSelector() {
  histf()->cd();
  histf()->cd("ObjectSelection");
  for (size_t i = 0; i < (*nFatJet->Get()); ++i){
    bool hasValidSubJets {false};
    AnaUtil::fillHist1D ("fatJetCutFlow", 0, 1.0);
    // preSelection
    if (!(FatJet_jetId->At(i) & 2)) continue;
    if (FatJet_pt->At(i) < 200) continue;
    if (std::fabs(FatJet_eta->At(i)) > 2.4) continue;
    AnaUtil::fillHist1D ("fatJetCutFlow", 1, 1.0);
    if ((FatJet_subJetIdx1->At(i) <= *nSubJet->Get() && (SubJet_pt->At(FatJet_subJetIdx1->At(i)) >= 20) && std::fabs(SubJet_eta->At(FatJet_subJetIdx1->At(i))) <= 2.4)
	&& (FatJet_subJetIdx2->At(i) <= *nSubJet->Get() && (SubJet_pt->At(FatJet_subJetIdx2->At(i)) >= 20) && std::fabs(SubJet_eta->At(FatJet_subJetIdx2->At(i))) <= 2.4)) 
      hasValidSubJets = true;
    if (!hasValidSubJets) continue;
    AnaUtil::fillHist1D ("fatJetCutFlow", 2, 1.0);
    if (FatJet_msoftdrop->At(i) < 30 || FatJet_msoftdrop->At(i) > 210) continue;
    AnaUtil::fillHist1D ("fatJetCutFlow", 3, 1.0);
    if (FatJet_tau2->At(i)/FatJet_tau1->At(i) > 0.75) continue;
    AnaUtil::fillHist1D ("fatJetCutFlow", 4, 1.0);

    vhtm::FatJet fj;
    fj.index                       = i;
    fj.pt                          = FatJet_pt->At(i);
    fj.eta                         = FatJet_eta->At(i);
    fj.phi                         = FatJet_phi->At(i);
    fj.mass                        = FatJet_mass->At(i);
    fj.softDropMass                = FatJet_msoftdrop->At(i);
    fj.n2b1                        = FatJet_n2b1->At(i);
    fj.n3b1                        = FatJet_n3b1->At(i);
    //fj.hadronFlavour               = FatJet_hadronFlavour->At(i);
    //fj.nBHadrons                   = FatJet_nBHadrons->At(i);
    //fj.nCHadrons                   = FatJet_nCHadrons->At(i);
    fj.rawFactor                   = FatJet_rawFactor->At(i);
    fj.subJetIdx1                  = FatJet_subJetIdx1->At(i);
    fj.subJetIdx2                  = FatJet_subJetIdx2->At(i);
    fj.btagDeepB                   = FatJet_btagDeepB->At(i);    
    fj.btagCSVV2                   = FatJet_btagCSVV2->At(i);
    fj.tau1                        = FatJet_tau1->At(i);
    fj.tau2                        = FatJet_tau2->At(i);
    fj.tau3                        = FatJet_tau3->At(i);
    fj.tau4                        = FatJet_tau4->At(i);
    fj.deepTag_WvsQCD              = FatJet_deepTag_WvsQCD->At(i);
    fj.deepTag_ZvsQCD              = FatJet_deepTag_ZvsQCD->At(i);
    fj.deepTag_TvsQCD              = FatJet_deepTag_TvsQCD->At(i);
    fj.deepTagMD_WvsQCD            = FatJet_deepTagMD_WvsQCD->At(i);
    fj.deepTagMD_ZvsQCD            = FatJet_deepTagMD_ZvsQCD->At(i);
    fj.deepTagMD_TvsQCD            = FatJet_deepTagMD_TvsQCD->At(i);
    //fj.electronIdx3SJ              = FatJet_electronIdx3SJ->At(i);
    //fj.muonIdx3SJ                  = FatJet_muonIdx3SJ->At(i);
    
    fatJetList_.push_back(fj);

    // fatJet cleaning wrt leptons
    if (!fatJetLeptonCleaning(fj)) continue;
    AnaUtil::fillHist1D ("fatJetCutFlow", 5, 1.0);
    cleanFatJetList_.push_back(fj);

    // these fatJets have b-tagged jets
    bool is_bTagged = ((SubJet_pt->At(FatJet_subJetIdx1->At(i)) >= 30 && SubJet_btagDeepB->At(FatJet_subJetIdx1->At(i)) > 0.4941) 
		       || (SubJet_pt->At(FatJet_subJetIdx2->At(i)) >= 30 && SubJet_btagDeepB->At(FatJet_subJetIdx2->At(i)) > 0.4941)) ? true : false;
    if (!is_bTagged) continue; 
    AnaUtil::fillHist1D ("fatJetCutFlow", 6, 1.0);
    bTaggedFatJetList_.push_back(fj);
  }
  searchedFatJet_ = true;
}

// FatJet selection
void PhysicsObjSelector::subJetSelector() {
  histf()->cd();
  histf()->cd("ObjectSelection");

  for (size_t i = 0; i < (*nSubJet->Get()); ++i){
    vhtm::SubJet sj;
    sj.index     = i;
    sj.pt        = SubJet_pt->At(i);
    sj.eta       = SubJet_eta->At(i);
    sj.phi       = SubJet_phi->At(i);
    sj.mass      = SubJet_mass->At(i);
    sj.btagDeepB = SubJet_btagDeepB->At(i);
    sj.rawFactor = SubJet_rawFactor->At(i);
    
    subJetList_.push_back(sj);
  }
}

bool PhysicsObjSelector::jetLeptonCleaning(const vhtm::Jet& jet) const {
  bool isMuon {false};
  bool isElectron {false};
  TLorentzVector jp4(AnaUtil::getP4(jet));
  // Medium Isolated muons
  for (const auto& mu: tightIsoMuList_){
    if (jp4.DeltaR(AnaUtil::getP4(mu)) <= 0.4) {
      isMuon = true;
      break;
    }
    if (isMuon) return false;
  }
  // LooseMVA Isolated electrons
  for (const auto& el: tightIsoEleList_){
    if (jp4.DeltaR(AnaUtil::getP4(el)) <= 0.4) {
      isElectron = true;
      break;
    }
    if (isElectron) return false;
  }
  return true;
}
// fat-jet lepton cleaning
bool PhysicsObjSelector::fatJetLeptonCleaning(const vhtm::FatJet& jet) const {
  bool isMuon {false};
  bool isElectron {false};
  TLorentzVector jp4(AnaUtil::getP4(jet));
  // Medium Isolated muons
  for (const auto& mu: tightIsoMuList_){
    if (jp4.DeltaR(AnaUtil::getP4(mu)) <= 0.8) {
      isMuon = true;
      break;
    }
    if (isMuon) return false;
  }
  // LooseMVA Isolated electrons
  for (const auto& el: tightIsoEleList_){
    if (jp4.DeltaR(AnaUtil::getP4(el)) <= 0.8) {
      isElectron = true;
      break;
    }
    if (isElectron) return false;
  }
  return true;
}

bool PhysicsObjSelector::tauLeptonCleaning(const vhtm::Tau& tau) const {
  bool isMuon {false};
  bool isElectron {false};
  TLorentzVector tp4(AnaUtil::getP4(tau));
  // Medium Isolated muons
  for (const auto& mu: tightIsoMuList_){
    if (tp4.DeltaR(AnaUtil::getP4(mu)) <= 0.3) {
      isMuon = true;
      break;
    }
    if (isMuon) return false;
  }
  // LooseMVA Isolated electrons
  for (const auto& el: tightIsoEleList_){
    if (tp4.DeltaR(AnaUtil::getP4(el)) <= 0.3) {
      isElectron = true;
      break;
    }
    if (isElectron) return false;
  }
  return true;
}

bool PhysicsObjSelector::thisElectronIsMuon(const vhtm::Electron& ele, bool VsLooseMuons,  bool VsTightMuons) const {
  bool isMuon {false};
  TLorentzVector elep4(AnaUtil::getP4(ele));
  if (VsLooseMuons) {
    for (const auto& mu: preSelMuList_){
      if (elep4.DeltaR(AnaUtil::getP4(mu)) <= 0.3) {
	isMuon  = true;
	break;
      }
      if (isMuon) return true;
    }
  }
  else if (VsTightMuons) {
    for (const auto& mu: tightMuList_){
      if (elep4.DeltaR(AnaUtil::getP4(mu)) <= 0.3) {
	isMuon  = true;
	break;
      }
      if (isMuon) return true;
    }
  }
  return false;
}
/*
bool PhysicsObjSelector::jetLeptonCleaning(const vhtm::Jet& jet) const {
  bool isMuon {false};
  bool isElectron {false};
  TLorentzVector jp4(AnaUtil::getP4(jet));
  // Medium Isolated muons
  for (const auto& mu: tightIsoMuList_){
    if (jet.muIdx1 == static_cast<int>(mu.index)) {
      AnaUtil::fillHist1D ("jmu1DR", jp4.DeltaR(AnaUtil::getP4(mu)), 1.0);
      isMuon = true;
      break;
    }
    else if (jet.muIdx2 == static_cast<int>(mu.index)) {
      AnaUtil::fillHist1D ("jmu2DR", jp4.DeltaR(AnaUtil::getP4(mu)), 1.0);
      isMuon = true;
      break;
    }
    if (isMuon) return false;
  }
  // LooseMVA Isolated electrons
  for (const auto& el: tightIsoEleList_){
    if (jet.elIdx1 == static_cast<int>(el.index)) {
      AnaUtil::fillHist1D ("jel1DR", jp4.DeltaR(AnaUtil::getP4(el)), 1.0);
      isElectron = true;
      break;
    }
    else if (jet.elIdx2 == static_cast<int>(el.index)) {
      AnaUtil::fillHist1D ("jel2DR", jp4.DeltaR(AnaUtil::getP4(el)), 1.0);
      isElectron = true;
      break;
    }
    if (isElectron) return false;
  }
  return true;
}
*/

/*
void PhysicsObjSelector::dumpEverything(int evNo, ostream& os) const {
  os << std::setprecision(3);

  // Event                                     
  os <<">>>Event Number: "<<evNo<<"\n";

  // Muons    
  if (*nMu->Get() > 0) {
    os << " -- # Muons: " << *nMu->Get() << endl;
    os << "  indx      pT     eta     phi  charge      dxy       dz  global tracker      PF         SIP   tightCharge   Loose    Medium    Tight   HighPtId   relIso_chr   relIso_all"
       << endl;
    for (size_t i = 0; i < (*nMu->Get()); ++i) {
      os << setw(6)  << i
         << setw(8)  << Mu_corrpt->At(i)
         << setw(8)  << Mu_eta->At(i)
         << setw(8)  << Mu_phi->At(i)
         << setw(8)  << Mu_charge->At(i)
         << setw(9)  << Mu_dxy->At(i)
         << setw(9)  << Mu_dz->At(i)
         << setw(8)  << (Mu_isGlobal->At(i) ? "T" : "F")
         << setw(8)  << (Mu_isPFcand->At(i) ? "T" : "F")
         << setw(8)  << (Mu_isTracker->At(i) ? "T" : "F")
         << setw(14) << Mu_sip3d->At(i)
	 << setw(8) << Mu_tightCharge->At(i)
         << setw(10) << (Mu_LooseId->At(i) ? "T" : "F")
         << setw(10) << (Mu_MediumId->At(i) ? "T" : "F")
         << setw(10) << (Mu_TightId->At(i) ? "T" : "F")
         << setw(10) << (Mu_highPtId->At(i) ? "T" : "F")
	 << setw(12) << Mu_pfRelIso03_chg->At(i)
         << setw(10) << Mu_pfRelIso03_all->At(i)
         << endl;
    }
    os <<"   Muon |||https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2?rev=26|||"
       <<"\n"
       <<"   |||---------------------------------------------------------------------------|||"
       <<endl;
  }

  // Electrons                                                                                 
  if (*nEle->Get() > 0) {
    os << " -- # Electrons: " << *nEle->Get() << endl;
    os << "  indx      pT     eta     phi  charge     dxy      dz     misHit       SIP3D   Loose(WP90) Tight(WP80)  relIso_chr  relIso_all"
       << endl;
    for (size_t ie = 0; ie < (*nEle->Get()); ++ie) {
      os << setw(6)  << ie
         << setw(8)  << Ele_pt->At(ie)
         << setw(8)  << Ele_eta->At(ie)
         << setw(8)  << Ele_phi->At(ie)
         << setw(8)  << Ele_charge->At(ie)
         << setw(8)  << Ele_dxy->At(ie)
         << setw(8)  << Ele_dz->At(ie)
         << setw(8)  << static_cast<int>(Ele_missHits->At(ie))
         << setw(16) << Ele_sip3d->At(ie)
         << setw(8)  << (Ele_mvaSpring16GP_WP90->At(ie) ? "T" : "F")
         << setw(13) << (Ele_mvaSpring16GP_WP80->At(ie) ? "T" : "F")
         << setw(13)  << Ele_pfRelIso03_chg->At(ie)
         << setw(12)  << Ele_pfRelIso03_all->At(ie)
         << endl;
    }
    os <<"   Electron |||https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRunIIRecommendations#Spring16_mva_Summer16_cut_based|||"
       <<"\n"
       <<"   Electron |||https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2#MVA_recipes_for_2016_data_and_Sp|||"
       <<"\n"
       <<"   Electron |||https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2?rev=%2039|||"
       <<"\n"
       <<"   |||---------------------------------------------------------------------------|||"
       <<endl;

  }

  // Taus                                                                            
  if (*nTau->Get()> 0) {
    os << " -- # Taus: " << *nTau->Get() << endl;
    os << "  indx       pT      eta      phi charge decayMode idDecayMode   isolation LmuVeto LeleVeto TmuVeto TeleVeto"
       << endl;
    for (size_t it = 0; it < (*nTau->Get()); ++it) {
      os << setw(6) << it
         << setw(9) << Tau_pt->At(it)
         << setw(9) << Tau_eta->At(it)
         << setw(9) << Tau_phi->At(it)
         << setw(7) << Tau_charge->At(it)
	 << setw(10)<< Tau_decayMode->At(it)
         << setw(10) << ((Tau_idDecayMode->At(it)) ? "T" : "F")
         << setw(10) << static_cast<int>(Tau_idMVAoldDM->At(it))
         << setw(9) << ((Tau_idAntiMu->At(it) == 1) ? "T" : "F")
         << setw(8) << ((Tau_idAntiEle->At(it) == 2) ? "T" : "F")
         << setw(9) << ((Tau_idAntiMu->At(it) == 2) ? "T" : "F")
         << setw(8) << ((Tau_idAntiEle->At(it) == 8) ? "T" : "F")
         << endl;
    }
    os <<"   Tau |||IsolationMVArun2v1DBoldDMwLT ID wp (2015): 1 = VLoose, 2 = Loose, 4 = Medium, 8 = Tight, 16 = VTight, 32 = VVTight |||"
       <<"\n"
       <<"   Tau |||Anti-electron MVA discriminator V6: bitmask 1 = VLoose, 2 = Loose, 4 = Medium, 8 = Tight, 16 = VTight|||"
       <<"\n"
       <<"   Tau |||Anti-muon discriminator V3: : bitmask 1 = Loose, 2 = Tight|||"
       <<"\n"
       <<"   |||---------------------------------------------------------------------------|||"
       <<endl;

  }

  //Jets
  if (*nJet->Get() > 0) {
    os << " -- # Jets: " << *nJet->Get() << endl;
    os << "  indx       pT      eta      phi NConst   nMu      nEle     CHF   CEMF      NHF     NEMF     puID   bDisc  looseId  tightId tightLepVeto"
       << endl;
    for (size_t j = 0; j < (*nJet->Get()); ++j) {
      os << setw(6)  << j
         << setw(9)  << Jet_nomPt->At(j)
         << setw(9)  << Jet_eta->At(j)
         << setw(9)  << Jet_phi->At(j)
         << setw(7)  << Jet_nConstituents->At(j)
         << setw(7)  << Jet_nMuons->At(j)
         << setw(7)  << Jet_nElectrons->At(j)
         << setw(9)  << Jet_chHEF->At(j)
         << setw(9)  << Jet_chEmEF->At(j)
         << setw(9)  << Jet_neHEF->At(j)
         << setw(9)  << Jet_neEmEF->At(j)
         << setw(5)  << Jet_puId->At(j)
         << setw(8)  << Jet_btagCSVV2->At(j) 
         << setw(9)  << ((Jet_jetId->At(j) >= 1) ? "T" : "F")
	 << setw(9)  << ((Jet_jetId->At(j) >= 3) ? "T" : "F")
         << setw(5)  << ((Jet_jetId->At(j) == 7) ? "T" : "F")
         << endl;
    }
    os <<"   Jet |||https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_2016|||"
       <<"\n"
       <<"   Jet |||https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD#Jets|||"
       <<"\n"
       <<"   Jet |||jetId==1: pass loose ID, fail tight, fail tightLepVeto|||"
       <<"\n"
       <<"   Jet |||jetId==3 means: pass loose and tight ID, fail tightLepVeto|||"
       <<"\n"
       <<"   Jet |||jetId==7 means: pass loose, tight, tightLepVeto ID|||"
       <<"\n"
       <<"   Jet ||| puId is valid when jetPt < 50 GeV {0: fail all PU ID | 4: pass loose ID, fail medium, fail tight | 6: pass loose and medium ID, fail tight | 7: pass all"
       <<"\n"
       <<"   |||---------------------------------------------------------------------------|||"
       <<"\n"
       <<endl;

  }
}
*/
