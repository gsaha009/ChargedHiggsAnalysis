#include "configana.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include <fstream>

#include "PhysicsObjSelector.h"
#include "AnaUtil.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::ostream;
using std::ostringstream;
using std::vector;
using std::map;
using std::pair;
using std::setprecision;
using std::setw;
using std::fabs;

PhysicsObjSelector::PhysicsObjSelector()
  : AnaBase()
{}
bool PhysicsObjSelector::beginJob() {
  if (!AnaBase::beginJob()) return false;

  histf()->cd();
  histf()->mkdir("ObjectSelection");
    
  return true;
}
void PhysicsObjSelector::bookHistograms() {
  histf()->cd();
  histf()->cd("ObjectSelection");

  new TH1D("muCutFlow", "Muon CutFlow", 10, -0.5, 9.5);
  new TH1D("eleCutFlow", "Electron CutFlow", 13, -0.5, 12.5);
  new TH1D("jetCutFlow", "Jet CutFlow", 6, -0.5, 5.5);
  new TH1D("tauCutFlow", "Tau CutFlow", 9, -0.5, 8.5);
  new TH1D("fatJetCutFlow", "FatJet CutFlow", 9, -0.5, 8.5);
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

  // show object selection cutflow
  // Muon
  {
    ostringstream ptCutTag;
    ptCutTag << "pT >= " << AnaUtil::cutValue(muonCutMap(), "pt") << " GeV";

    ostringstream etaCutTag;
    etaCutTag << "|eta| < " << AnaUtil::cutValue(muonCutMap(), "eta");

    ostringstream dxyCutTag;
    dxyCutTag << "|dxy| < " << AnaUtil::cutValue(muonCutMap(), "dxy") << " cm";

    ostringstream dzCutTag;
    dzCutTag << "|dz| < " << AnaUtil::cutValue(muonCutMap(), "dz") << " cm";

    ostringstream sipCutTag;
    sipCutTag << "|sip3d| <= " << AnaUtil::cutValue(muonCutMap(), "SIP3D");

    ostringstream ptFakeCutTag;
    ptFakeCutTag << "pTFake >= " << AnaUtil::cutValue(muonCutMap(), "ptFake") << " GeV";

    ostringstream isoCutTag;
    isoCutTag << "pfRelIso <  " << AnaUtil::cutValue(muonCutMap(), "reliso");

    vector<string> muLabels {
      "nMuonCand",
      ptCutTag.str(),
      etaCutTag.str(),
      dxyCutTag.str(),
      dzCutTag.str(),
      sipCutTag.str(),
      "LooseId",
      ptFakeCutTag.str(),
      isoCutTag.str(),
      "MediumId"
    };
    AnaUtil::showEfficiency("muCutFlow", muLabels, "Muon Selection", "Muons");  
  }
  // Electron
  {
    ostringstream ptCutTag;
    ptCutTag << "pT >= " << AnaUtil::cutValue(electronCutMap(), "pt") << " GeV";
    ostringstream etaCutTag;
    etaCutTag << "|eta| < " << AnaUtil::cutValue(electronCutMap(), "eta");

    ostringstream dxyCutTag;
    dxyCutTag << "|dxy| < " << AnaUtil::cutValue(electronCutMap(), "dxy") << " cm";

    ostringstream dzCutTag;
    dzCutTag << "|dz| < " << AnaUtil::cutValue(electronCutMap(), "dz") << " cm";

    ostringstream sipCutTag;
    sipCutTag << "|sip3d| <= " << AnaUtil::cutValue(electronCutMap(), "SIP3D");

    ostringstream ptFakeCutTag;
    ptFakeCutTag << "pTFake >= " << AnaUtil::cutValue(electronCutMap(), "ptFake") << " GeV";

    ostringstream isoCutTag;
    isoCutTag << "pfRelIso <  " << AnaUtil::cutValue(electronCutMap(), "reliso");

    vector<string> eleLabels {
      "nElectronCand",
      ptCutTag.str(),
      etaCutTag.str(),
      dxyCutTag.str(),
      dzCutTag.str(),
      sipCutTag.str(),
      "mvaFall17V2Iso_WP80",
      "cleaned against tight muons",
      ptFakeCutTag.str(),
      "ConvVeto",
      "no LostHits",
      isoCutTag.str(),
      "mvaFall17V2Iso_WP90"
    };
    AnaUtil::showEfficiency("eleCutFlow", eleLabels, "Electron Selection", "Electrons");  
  }
  // Jet
  {
    ostringstream ptCutTag;
    ptCutTag << "pT >= " << AnaUtil::cutValue(jetCutMap(), "pt") << " GeV";

    ostringstream etaCutTag;
    etaCutTag << "|eta| < " << AnaUtil::cutValue(jetCutMap(), "eta");

    vector<string> jetLabels {
      "nJetCand",
      ptCutTag.str(),
      etaCutTag.str(),
      "puId cut",
      "lepton cleaned",
      "b-Tagged"
    };
    AnaUtil::showEfficiency("jetCutFlow", jetLabels, "Jet Selection", "Jets");  
  }
  // Tau
  {
    ostringstream ptCutTag;
    ptCutTag << "pT >= " << AnaUtil::cutValue(tauCutMap(), "pt") << " GeV";

    ostringstream etaCutTag;
    etaCutTag << "|eta| < " << AnaUtil::cutValue(tauCutMap(), "eta");

    ostringstream dxyCutTag;
    dxyCutTag << "|dxy| < " << AnaUtil::cutValue(tauCutMap(), "dxy") << " cm";

    ostringstream dzCutTag;
    dzCutTag << "|dz| < " << AnaUtil::cutValue(tauCutMap(), "dz") << " cm";

    vector<string> tauLabels {
      "nTauCand",
      ptCutTag.str(),
      etaCutTag.str(),
      dxyCutTag.str(),
      dzCutTag.str(),
      "Hadronic Tau decay mode",
      "Allowed decay modes",
      "lepton and jet veto",
      "lepCleaned"
    };
    AnaUtil::showEfficiency("tauCutFlow", tauLabels, "Tau Selection", "Taus");  
  }
  // FatJet
  {
    ostringstream ptCutTag;
    ptCutTag << "pT >= " << AnaUtil::cutValue(fatJetCutMap(), "pt") << " GeV";

    ostringstream etaCutTag;
    etaCutTag << "|eta| < " << AnaUtil::cutValue(fatJetCutMap(), "eta");

    ostringstream msdCutTag;
    msdCutTag << AnaUtil::cutValue(fatJetCutMap(), "msoftdropMin") << " <= msoftdtop <= " 
              << AnaUtil::cutValue(fatJetCutMap(), "msoftdropMax") << " GeV";

    ostringstream trCutTag;
    trCutTag << " fjetTauRatio <= " << AnaUtil::cutValue(fatJetCutMap(), "fjetRatio");

    vector<string> fatJetLabels {
      "nFatJetCand",
      "jetId",
      ptCutTag.str(),
      etaCutTag.str(),
      "hasValidSubJets",
      msdCutTag.str(),
      trCutTag.str(),
      "leptonCleaning",
      "b-Tagged"  
    };
    AnaUtil::showEfficiency("fatJetCutFlow", fatJetLabels, "FatJet Selection", "FatJets");  
  }
}
// clear lists
void PhysicsObjSelector::clear() {
  eventList_.clear();
  metList_.clear();

  // Muons
  preSelMuList_.clear();
  fakeableMuList_.clear();
  tightMuList_.clear();

  // Electrons
  preSelEleList_.clear();
  fakeableEleList_.clear();
  tightEleList_.clear();

  // Jets
  preSelJetList_.clear();
  leptonCleanJetList_.clear();
  leptonCleanJetListOutsideAk8_.clear();
  looseBJetList_.clear();
  bJetList_.clear();

  // Taus
  tauList_.clear();
  leptonCleanTauList_.clear();

  // fat Jets
  fatJetList_.clear();
  cleanFatJetList_.clear();
  bTaggedFatJetList_.clear();

  // subject
  subJetList_.clear();

  // Truth
  genParticleList_.clear();
  lheParticleList_.clear();

  // other flags
  searchedMu_  = false; 
  searchedEle_ = false; 
  searchedJet_ = false;
  searchedFatJet_ = false;
}
double PhysicsObjSelector::getGenSumW() {
  return *genEventSumw_->Get();
}
bool PhysicsObjSelector::findEventInfo() {
  vhtm::Event ev;
  ev.run    = *run_->Get();
  ev.event  = *event_->Get();
  ev.lumis  = *lumis_->Get();
  if (isMC()) {
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
bool PhysicsObjSelector::findGenPartInfo() {
  for (size_t j = 0; j < *nGenPart->Get(); ++j) {      
    vhtm::GenParticle gp;
    
    gp.index = j;
    gp.pt   = GenPart_pt->At(j);
    gp.eta  = GenPart_eta->At(j);
    gp.phi  = GenPart_phi->At(j);
    gp.mass = GenPart_mass->At(j);
    gp.pdgId       = GenPart_pdgId->At(j);
    gp.motherIdx   = GenPart_motherIdx->At(j);
    gp.status      = GenPart_status->At(j);
    gp.statusFlags = GenPart_statusFlags->At(j);
    
    genParticleList_.push_back(gp);
  }
  return true;
}

bool PhysicsObjSelector::findLHEPartInfo() {
  for (size_t j = 0; j < *nLHEPart->Get(); ++j) {
    vhtm::LHEParticle gp;
      
    gp.index = j;
    gp.pt    = LHEPart_pt->At(j);
    gp.eta   = LHEPart_eta->At(j);
    gp.phi   = LHEPart_phi->At(j);
    gp.mass  = LHEPart_mass->At(j);
    gp.pdgId = LHEPart_pdgId->At(j);
    
    lheParticleList_.push_back(gp);
  }
  return true;
}
void PhysicsObjSelector::findObjects() {
  histf()->cd();
  histf()->cd("ObjectSelection");

  // the order is usually important
  muonSelector();
  electronSelector();
  fatJetSelector(); // order important
  jetSelector();
  subJetSelector();
  metSelector();
  tauSelector();
}
// Muon selection (already sorted in pT)
void PhysicsObjSelector::muonSelector() {
  static double ptMin     = AnaUtil::cutValue(muonCutMap(), "pt"),
                etaMax    = AnaUtil::cutValue(muonCutMap(), "eta"),
                dxyMax    = AnaUtil::cutValue(muonCutMap(), "dxy"),
                dzMax     = AnaUtil::cutValue(muonCutMap(), "dz"),
                sip3dMax  = AnaUtil::cutValue(muonCutMap(), "SIP3D"),
                ptFakeMin = AnaUtil::cutValue(muonCutMap(), "ptFake"),
                relisoMax = AnaUtil::cutValue(muonCutMap(), "reliso");

  for (size_t j = 0; j < *nMuon->Get(); ++j) {
    AnaUtil::fillHist1D("muCutFlow", 0);
  
    if (getEra() == 2018 && Muon_corrpt->At(j) <= ptMin) continue;
    AnaUtil::fillHist1D("muCutFlow", 1);

    if (fabs(Muon_eta->At(j)) >= etaMax) continue;
    AnaUtil::fillHist1D("muCutFlow", 2);

    if (getEra() == 2018 && fabs(Muon_dxy->At(j)) >= dxyMax) continue;
    if (getEra() == 2016 && fabs(Muon_dxy->At(j)) >= 0.05) continue;
    AnaUtil::fillHist1D("muCutFlow", 3);

    if (getEra() == 2018 && fabs(Muon_dz->At(j)) >= dzMax) continue;
    if (getEra() == 2016 && fabs(Muon_dz->At(j)) >= 0.1) continue;
    AnaUtil::fillHist1D("muCutFlow", 4);

    if (fabs(Muon_sip3d->At(j)) > sip3dMax) continue;
    AnaUtil::fillHist1D("muCutFlow", 5);

    if (!Muon_LooseId->At(j)) continue;
    AnaUtil::fillHist1D("muCutFlow", 6);

    vhtm::Muon mu;
    mu.index   = j;
    mu.pt      = Muon_corrpt->At(j);
    mu.eta     = Muon_eta->At(j);
    mu.phi     = Muon_phi->At(j);
    mu.mass    = Muon_mass->At(j);
    mu.charge  = Muon_charge->At(j);
    mu.jetIdx  = Muon_jetIdx->At(j);
    if (isMC()) {
      mu.genIdx = Muon_genPartIdx->At(j);
      // genFlv = 1  : promt muon, genFlv = 15 : from tau decay
      mu.genFlv = Muon_genPartFlv->At(j); 
    }
    mu.mediumId  = Muon_MediumId->At(j);
    mu.tightId   = Muon_TightId->At(j);
    mu.highPtId  = Muon_highPtId->At(j);
    mu.pfRelIso03_all = Muon_pfRelIso03_all->At(j);
    mu.pfRelIso04_all = Muon_pfRelIso04_all->At(j);
    mu.pdgId     = Muon_pdgId->At(j);

    preSelMuList_.push_back(mu);

    // Fakeable Muon Selection
    if (Muon_corrpt->At(j) < ptFakeMin) continue;
    AnaUtil::fillHist1D("muCutFlow", 7);

    // PF Isolation
    if (getEra() == 2018 && Muon_pfRelIso03_all->At(j) > relisoMax) continue;
    if (getEra() == 2016 && Muon_pfRelIso03_all->At(j) > 0.15) continue;
    AnaUtil::fillHist1D("muCutFlow", 8);
    fakeableMuList_.push_back(mu);

    // Tight Muon selection
    //if (!Muon_TightId->At(j)) continue;
    if (!Muon_MediumId->At(j)) continue;
    AnaUtil::fillHist1D("muCutFlow", 9);

    tightMuList_.push_back(mu);
  }
  searchedMu_ = true;
}
// Electron selecton (already sorted in pT)
void PhysicsObjSelector::electronSelector() {
  if (!searchedMu_) {
    cout << ">>> Muons are not selected yet!" << endl;
    muonSelector();
  }
  static double ptMin     = AnaUtil::cutValue(electronCutMap(), "pt"),
                etaMax    = AnaUtil::cutValue(electronCutMap(), "eta"),
                dxyMax    = AnaUtil::cutValue(electronCutMap(), "dxy"),
                dzMax     = AnaUtil::cutValue(electronCutMap(), "dz"),
                sip3dMax  = AnaUtil::cutValue(electronCutMap(), "SIP3D"),
                ptFakeMin = AnaUtil::cutValue(electronCutMap(), "ptFake"),
                relisoMax = AnaUtil::cutValue(electronCutMap(), "reliso");

  for (size_t j = 0; j < *nElectron->Get(); ++j) {
    AnaUtil::fillHist1D("eleCutFlow", 0);

    if (Electron_pt->At(j) < ptMin) continue;
    AnaUtil::fillHist1D("eleCutFlow", 1);

    if (fabs(Electron_eta->At(j)) > etaMax) continue;
    AnaUtil::fillHist1D("eleCutFlow", 2);

    if (getEra() == 2018 && fabs(Electron_dxy->At(j)) > dxyMax) continue;
    if (getEra() == 2016)
      if ((std::fabs(Electron_eta->At(j)) <= 1.479
	   && std::fabs(Electron_dxy->At(j)) >= 0.05)||
	  ((std::fabs(Electron_eta->At(j)) > 1.479)
	   && std::fabs(Electron_dxy->At(j)) >= 0.1)) continue;
    AnaUtil::fillHist1D("eleCutFlow", 3);

    if (getEra() == 2018 && fabs(Electron_dz->At(j)) > dzMax) continue;
    if (getEra() == 2016)
      if ((std::fabs(Electron_eta->At(j)) <= 1.479
	   && std::fabs(Electron_dz->At(j)) >= 0.1)  ||
	  ((std::fabs(Electron_eta->At(j)) > 1.479)
	   && std::fabs(Electron_dz->At(j)) >= 0.2)) continue;
    AnaUtil::fillHist1D("eleCutFlow", 4);
    
    if (Electron_sip3d->At(j) > sip3dMax) continue;
    AnaUtil::fillHist1D("eleCutFlow", 5);

    if (getEra() == 2018 && !Electron_mvaFall17V2noIso_WP80->At(j)) continue;
    if (getEra() == 2016 && !Electron_mvaSpring16GP_WP80->At(j)) continue;
    AnaUtil::fillHist1D("eleCutFlow", 6);

    vhtm::Electron el;
    el.index   = j;
    el.pt      = Electron_pt->At(j);
    el.eta     = Electron_eta->At(j);
    el.SCdEta  = Electron_deltaEtaSC->At(j);
    el.SCeta   = Electron_eta->At(j) + Electron_deltaEtaSC->At(j); /// ???
    el.phi     = Electron_phi->At(j);
    el.mass    = Electron_mass->At(j);
    el.charge  = Electron_charge->At(j);
    el.jetIdx  = Electron_jetIdx->At(j);
    el.phoIdx  = Electron_phoIdx->At(j);
    if (isMC()) {
      el.genIdx = Electron_genPartIdx->At(j);
      el.genFlv = Electron_genPartFlv->At(j);
    }
    el.pdgId   = Electron_pdgId->At(j);

    //    if (thisElectronIsMuon(el, false, true)) continue;
    if (getEra() == 2018 && thisElectronIsMuon(el, tightMuList_)) continue;
    AnaUtil::fillHist1D("eleCutFlow", 7, 1.0); 
    preSelEleList_.push_back(el);

    if (Electron_pt->At(j) < ptFakeMin) continue;
    AnaUtil::fillHist1D("eleCutFlow", 8); 

    //if (!Electron_convVeto->At(j)) continue;;
    AnaUtil::fillHist1D("eleCutFlow", 9); 

    //if (Electron_lostHits->At(j) != 0) continue;
    AnaUtil::fillHist1D("eleCutFlow", 10); 

    if (getEra() == 2018 && Electron_pfRelIso03_all->At(j) > relisoMax) continue;
    if (getEra() == 2016 && Electron_pfRelIso03_all->At(j) > 0.2) continue;
    AnaUtil::fillHist1D("eleCutFlow", 11);     
    fakeableEleList_.push_back(el);

    if (getEra() == 2018 && !Electron_mvaFall17V2Iso_WP90->At(j)) continue;
    if (getEra() == 2016 && !Electron_mvaSpring16GP_WP90->At(j)) continue;
    AnaUtil::fillHist1D("eleCutFlow", 12);
    tightEleList_.push_back(el);
  }
  searchedEle_ = true;
}
// Jet selection
void PhysicsObjSelector::jetSelector() {
  if (!searchedMu_) {
    cerr << ">>> Muons have not been selected yet!" << endl;
    muonSelector();
  }  
  if (!searchedEle_) {
    cerr << ">>> Electron have not been selected yet!" << endl;
    electronSelector();
  }
  if (!searchedFatJet_) {
    cerr << ">>> Ak4 Jets have not been selected yet!" << endl;
    fatJetSelector();
  }
  static double ptMin      = AnaUtil::cutValue(jetCutMap(), "pt"),
                etaMax     = AnaUtil::cutValue(jetCutMap(), "eta"),
                lbjetScore = AnaUtil::cutValue(jetCutMap(), "lbjetScore"),
                bjetScore  = AnaUtil::cutValue(jetCutMap(), "bjetScore"),
                deltaRMax  = AnaUtil::cutValue(jetCutMap(), "deltaR");

  for (size_t j = 0; j < *nJet->Get(); ++j) {
    AnaUtil::fillHist1D("jetCutFlow", 0);

    if (!(Jet_jetId->At(j) & 0x2)) continue;
    AnaUtil::fillHist1D("jetCutFlow", 1);

    if (Jet_nomPt->At(j) < ptMin) continue;
    AnaUtil::fillHist1D("jetCutFlow", 2);

    if (fabs(Jet_eta->At(j)) > etaMax) continue;
    AnaUtil::fillHist1D("jetCutFlow", 3);

    // Apply PuID on Jets with pT < 50 GeV (Need to be discussed!!!)
    if (Jet_nomPt->At(j) <= 50 && !(Jet_puId->At(j) >> 2 & 0x1))  continue;
    AnaUtil::fillHist1D("jetCutFlow", 4);

    vhtm::Jet jet;
    jet.index   = j;
    jet.pt      = Jet_nomPt->At(j);
    jet.eta     = Jet_eta->At(j);
    jet.phi     = Jet_phi->At(j);
    jet.mass    = Jet_mass->At(j);
    jet.muIdx1  = Jet_muIdx1->At(j);
    jet.muIdx2  = Jet_muIdx2->At(j);
    jet.elIdx1  = Jet_elIdx1->At(j);
    jet.elIdx2  = Jet_elIdx2->At(j);
    jet.btagDeepFlavB = Jet_btagDeepFlavB->At(j);
    jet.btagCSVV2 = Jet_btagCSVV2->At(j);

    preSelJetList_.push_back(jet);

    if (!jetLeptonCleaning(jet)) continue;
    AnaUtil::fillHist1D("jetCutFlow", 5);
    leptonCleanJetList_.push_back(jet);

    // making b_jet collection
    if (Jet_btagDeepFlavB->At(j) >= lbjetScore)
      looseBJetList_.push_back(jet);

    if (Jet_btagDeepFlavB->At(j) >= bjetScore) {
      bJetList_.push_back(jet);
      AnaUtil::fillHist1D("jetCutFlow", 6);
    }
    // ak8 cleaned jet list
    TLorentzVector jetp4 = AnaUtil::getP4(jet);
    if (cleanFatJetList_.size() == 1) {
      if (jetp4.DeltaR(AnaUtil::getP4(cleanFatJetList_[0])) > deltaRMax)
	leptonCleanJetListOutsideAk8_.push_back(jet);
    }
    else if (cleanFatJetList_.size() >= 2) {
      if (jetp4.DeltaR(AnaUtil::getP4(cleanFatJetList_[0])) > deltaRMax &&
	  jetp4.DeltaR(AnaUtil::getP4(cleanFatJetList_[1])) > deltaRMax)
	leptonCleanJetListOutsideAk8_.push_back(jet);
    }
  }
  searchedJet_ = true;
}
// MET selection
void PhysicsObjSelector::metSelector() {
  vhtm::MET met;
  met.pt    = isSignal() ? *Met_nomPt->Get()  : *Met_pt->Get();
  met.phi   = isSignal() ? *Met_nomPhi->Get() : *Met_phi->Get();
  met.signf = *Met_significance->Get();
  met.sumEt = *Met_sumEt->Get();

  metList_.push_back(met);
}
// Tau selection
void PhysicsObjSelector::tauSelector() {
  if (!searchedJet_) jetSelector();
  static double ptMin  = AnaUtil::cutValue(tauCutMap(), "pt"),
                etaMax = AnaUtil::cutValue(tauCutMap(), "eta"),
                dxyMax = AnaUtil::cutValue(tauCutMap(), "dxy"),
                dzMax  = AnaUtil::cutValue(tauCutMap(), "dz");

  for (size_t j = 0; j < *nTau->Get(); ++j) {
    AnaUtil::fillHist1D("tauCutFlow", 0);

    if (Tau_pt->At(j) < ptMin) continue;
    AnaUtil::fillHist1D("tauCutFlow", 1);

    if (fabs(Tau_eta->At(j)) > etaMax) continue;
    AnaUtil::fillHist1D("tauCutFlow", 2);

    if (fabs(Tau_dxy->At(j)) > dxyMax) continue;
    AnaUtil::fillHist1D("tauCutFlow", 3);

    if (fabs(Tau_dz->At(j)) > dzMax) continue;
    AnaUtil::fillHist1D("tauCutFlow", 4);

    if (!Tau_idDecayModeNewDMs->At(j)) continue;
    AnaUtil::fillHist1D("tauCutFlow", 5);

    vector<int> dmodes {0, 1, 2, 10, 11};
    if (std::find(dmodes.begin(), dmodes.end(), Tau_decayMode->At(j)) == dmodes.end()) continue;
    AnaUtil::fillHist1D("tauCutFlow", 6);

    if (isSignal()) {
      if (!(Tau_idDeepTau2017v2VSjet->At(j) >> 4 & 0x1 &&
	    Tau_idDeepTau2017v2VSe->At(j)   & 0x1 &&
	    Tau_idDeepTau2017v2VSmu->At(j)  & 0x1)) continue;
    }
    else {
      if (!(Tau_idDeepTau2017v2p1VSjet->At(j) >> 4 & 0x1 &&
	    Tau_idDeepTau2017v2p1VSe->At(j)   & 0x1 &&
	    Tau_idDeepTau2017v2p1VSmu->At(j)  & 0x1)) continue;
    }
    AnaUtil::fillHist1D("tauCutFlow", 7);

    vhtm::Tau obj;
    obj.index   = j;
    obj.pt      = Tau_pt->At(j);
    obj.eta     = Tau_eta->At(j);
    obj.phi     = Tau_phi->At(j);
    obj.charge  = Tau_charge->At(j);
    obj.mass    = Tau_mass->At(j);
    obj.jetIdx  = Tau_jetIdx->At(j);

    tauList_.push_back(obj);
    
    if (!tauLeptonCleaning(obj)) continue;
    AnaUtil::fillHist1D("tauCutFlow", 8);
    leptonCleanTauList_.push_back(obj);
  }
}
// FatJet selection
void PhysicsObjSelector::fatJetSelector() {
  static double ptMin    = AnaUtil::cutValue(fatJetCutMap(), "pt"),
         etaMax          = AnaUtil::cutValue(fatJetCutMap(), "eta"),
         subJetPtMin     = AnaUtil::cutValue(fatJetCutMap(), "subJetPt"),
         subJetEtaMax    = AnaUtil::cutValue(fatJetCutMap(), "subJetEta"),
         msdMin          = AnaUtil::cutValue(fatJetCutMap(), "msoftdropMin"),
         msdMax          = AnaUtil::cutValue(fatJetCutMap(), "msoftdropMax"),
         fjetRatioMax    = AnaUtil::cutValue(fatJetCutMap(), "fjetRatio"),
         bSubJetPtMin    = AnaUtil::cutValue(fatJetCutMap(), "bSubJetPt"),
         bSubJetScoreMin = AnaUtil::cutValue(fatJetCutMap(), "bSubJetScore");

  for (size_t j = 0; j < *nFatJet->Get(); ++j) {
    AnaUtil::fillHist1D("fatJetCutFlow", 0);

    if (!(FatJet_jetId->At(j) & 0x2)) continue;
    AnaUtil::fillHist1D("fatJetCutFlow", 1);

    if (FatJet_pt->At(j) < ptMin) continue;
    AnaUtil::fillHist1D("fatJetCutFlow", 2);

    if (fabs(FatJet_eta->At(j)) > etaMax) continue;
    AnaUtil::fillHist1D("fatJetCutFlow", 3);

    int nsjet = static_cast<int>(*nSubJet->Get());
    int idx1 = FatJet_subJetIdx1->At(j);
    int idx2 = FatJet_subJetIdx2->At(j);
    bool hasValidSubJets = (
          (idx1 <= nsjet && SubJet_pt->At(idx1) >= subJetPtMin 
	                 && fabs(SubJet_eta->At(idx1)) <= subJetEtaMax)
       && (idx2 <= nsjet && SubJet_pt->At(idx2) >= subJetPtMin 
	                 && fabs(SubJet_eta->At(idx2)) <= subJetEtaMax)
    );
    if (!hasValidSubJets) continue;
    AnaUtil::fillHist1D("fatJetCutFlow", 4);

    float msd = FatJet_msoftdrop->At(j);
    if (msd < msdMin || msd > msdMax) continue;
    AnaUtil::fillHist1D("fatJetCutFlow", 5);

    // will the ratio be always valid?
    if (FatJet_tau2->At(j)/FatJet_tau1->At(j) > fjetRatioMax) continue;
    AnaUtil::fillHist1D("fatJetCutFlow", 6);

    vhtm::FatJet jet;
    jet.index   = j;
    jet.pt      = FatJet_pt->At(j);
    jet.eta     = FatJet_eta->At(j);
    jet.phi     = FatJet_phi->At(j);
    jet.mass    = FatJet_mass->At(j);
    jet.softDropMass  = msd;
    jet.n2b1    = FatJet_n2b1->At(j);
    jet.n3b1    = FatJet_n3b1->At(j);
    // jet.hadronFlavour = FatJet_hadronFlavour->At(j);
    // jet.nBHadrons = FatJet_nBHadrons->At(j);
    // jet.nCHadrons = FatJet_nCHadrons->At(j);
    jet.rawFactor = FatJet_rawFactor->At(j);
    jet.subJetIdx1 = FatJet_subJetIdx1->At(j);
    jet.subJetIdx2 = FatJet_subJetIdx2->At(j);
    jet.btagDeepB = FatJet_btagDeepB->At(j);
    jet.btagCSVV2 = FatJet_btagCSVV2->At(j);
    jet.tau1 = FatJet_tau1->At(j);
    jet.tau2 = FatJet_tau2->At(j);
    jet.tau3 = FatJet_tau3->At(j);
    jet.tau4 = FatJet_tau4->At(j);
    jet.deepTag_WvsQCD = FatJet_deepTag_WvsQCD->At(j);
    jet.deepTag_ZvsQCD = FatJet_deepTag_ZvsQCD->At(j);
    jet.deepTag_TvsQCD = FatJet_deepTag_TvsQCD->At(j);
    jet.deepTagMD_WvsQCD = FatJet_deepTagMD_WvsQCD->At(j);
    jet.deepTagMD_ZvsQCD = FatJet_deepTagMD_ZvsQCD->At(j);
    jet.deepTagMD_TvsQCD = FatJet_deepTagMD_TvsQCD->At(j);
    //jet.electronIdx3SJ = FatJet_electronIdx3SJ->At(j);
    //jet.muonIdx3SJ     = FatJet_muonIdx3SJ->At(j);
    
    // fatJet cleaning wrt leptons
    if (!fatJetLeptonCleaning(jet)) continue;
    AnaUtil::fillHist1D("fatJetCutFlow", 7);
    cleanFatJetList_.push_back(jet);

    // these fatJets have b-tagged jets
    bool is_bTagged = (
       (SubJet_pt->At(idx1) >= bSubJetPtMin && SubJet_btagDeepB->At(idx1) > bSubJetScoreMin) ||
       (SubJet_pt->At(idx2) >= bSubJetPtMin && SubJet_btagDeepB->At(idx2) > bSubJetScoreMin) 
    ) ? true : false;
    if (!is_bTagged) continue; 
    AnaUtil::fillHist1D("fatJetCutFlow", 8);
    bTaggedFatJetList_.push_back(jet);
  }
  searchedFatJet_ = true;
}
// SubJet selection
void PhysicsObjSelector::subJetSelector() {
  for (size_t j = 0; j < *nSubJet->Get(); ++j) {
    vhtm::SubJet jet;
    jet.index     = j;
    jet.pt        = SubJet_pt->At(j);
    jet.eta       = SubJet_eta->At(j);
    jet.phi       = SubJet_phi->At(j);
    jet.mass      = SubJet_mass->At(j);
    jet.btagDeepB = SubJet_btagDeepB->At(j);
    jet.rawFactor = SubJet_rawFactor->At(j);
    
    subJetList_.push_back(jet);
  }
}
bool PhysicsObjSelector::jetLeptonCleaning(const vhtm::Jet& jet, double minDR) const {
  TLorentzVector jp4(AnaUtil::getP4(jet));
  for (const auto& mu: fakeableMuList_)  if (jp4.DeltaR(AnaUtil::getP4(mu)) <= minDR) return false;
  for (const auto& el: fakeableEleList_) if (jp4.DeltaR(AnaUtil::getP4(el)) <= minDR) return false;

  return true;
}
// fat-jet lepton cleaning
bool PhysicsObjSelector::fatJetLeptonCleaning(const vhtm::FatJet& jet, double minDR) const {
  TLorentzVector jp4(AnaUtil::getP4(jet));
  for (const auto& mu: tightMuList_)  if (jp4.DeltaR(AnaUtil::getP4(mu)) <= minDR) return false;
  for (const auto& el: tightEleList_) if (jp4.DeltaR(AnaUtil::getP4(el)) <= minDR) return false;

  return true;
}
bool PhysicsObjSelector::tauLeptonCleaning(const vhtm::Tau& tau, double minDR) const {
  TLorentzVector tp4(AnaUtil::getP4(tau));
  for (const auto& mu: tightMuList_)  if (tp4.DeltaR(AnaUtil::getP4(mu)) <= minDR) return false;
  for (const auto& el: tightEleList_) if (tp4.DeltaR(AnaUtil::getP4(el)) <= minDR) return false;

  return true;
}
/*
bool PhysicsObjSelector::thisElectronIsMuon(const vhtm::Electron& ele, bool VsLooseMuons, bool VsTightMuons, double minDR) const {
  TLorentzVector elep4(AnaUtil::getP4(ele));
  if (VsLooseMuons) {
    for (const auto& mu: preSelMuList_) if (elep4.DeltaR(AnaUtil::getP4(mu)) <= minDR) return true;
  }
  else if (VsTightMuons) {
    for (const auto& mu: tightMuList_)  if (elep4.DeltaR(AnaUtil::getP4(mu)) <= minDR) return true;
  }
  return false;
}
*/
bool PhysicsObjSelector::thisElectronIsMuon(const vhtm::Electron& ele, std::vector<vhtm::Muon> muonList, double minDR) const {
  TLorentzVector elep4(AnaUtil::getP4(ele));
  for (const auto& mu: muonList)  
    if (elep4.DeltaR(AnaUtil::getP4(mu)) <= minDR) 
      return true;
  return false;
}
void PhysicsObjSelector::dumpEvent(int evNo, ostream& os) const {
  os << std::setprecision(3);

  // Event                                     
  os << ">>> Event Number: " << evNo << endl;

  // Primary Vertex
  os << " -- Primary Vertex information: " << endl;
  os << "   nPV  nGoodPV  PVscore   PVchi2    PVndf      PVx      PVy      PVz" << endl;
  os << setw(6) << *nPV->Get()
     << setw(9) << *nGoodPV->Get()
     << setw(9) << static_cast<int>(*PV_score->Get())
     << setw(9) << *PV_chi2->Get()
     << setw(9) << *PV_ndf->Get()
     << setw(9) << *PV_xPos->Get()
     << setw(9) << *PV_yPos->Get()
     << setw(9) << *PV_zPos->Get()
     << endl;

  // Muons 
  size_t nMuonCand = *nMuon->Get();
  if (nMuonCand > 0) {
    os << " -- #Muons: " << nMuonCand << endl;
    os << "  indx      pT     eta     phi  charge      dxy       dz  global tracker      PF      SIP3D tightChg Loose"
       << " Medium Tight HighPt   IsoChg   IsoAll jetIdx"
       << endl;
    for (size_t j = 0; j < nMuonCand; ++j) {
      os << setw(6)  << j
         << setw(8)  << Muon_corrpt->At(j)
         << setw(8)  << Muon_eta->At(j)
         << setw(8)  << Muon_phi->At(j)
         << setw(8)  << Muon_charge->At(j)
         << setw(9)  << Muon_dxy->At(j)
         << setw(9)  << Muon_dz->At(j)
         << setw(8)  << (Muon_isGlobal->At(j) ? "T" : "F")
         << setw(8)  << (Muon_isPFcand->At(j) ? "T" : "F")
         << setw(8)  << (Muon_isTracker->At(j) ? "T" : "F")
         << setw(11) << Muon_sip3d->At(j)
	 << setw(9)  << Muon_tightCharge->At(j)
         << setw(6)  << (Muon_LooseId->At(j) ? "T" : "F")
         << setw(7)  << (Muon_MediumId->At(j) ? "T" : "F")
         << setw(6)  << (Muon_TightId->At(j) ? "T" : "F")
         << setw(7)  << (Muon_highPtId->At(j) ? "T" : "F")
	 << setw(9)  << Muon_pfRelIso03_chg->At(j)
         << setw(9)  << Muon_pfRelIso03_all->At(j)
         << setw(7)  << Muon_jetIdx->At(j)
         << endl;
    }
  }
  // Electrons                                                                                 
  size_t nElectronCand = *nElectron->Get();
  if (nElectronCand > 0) {
    os << " -- #Electrons: " << nElectronCand << endl;
    os << "  indx      pT     eta     phi  charge     dxy      dz  lostHit      SIP3D  LWP90  TWP80   IsoChg   IsoAll"
       << endl;
    for (size_t j = 0; j < nElectronCand; ++j) {
      os << setw(6)  << j
         << setw(8)  << Electron_pt->At(j)
         << setw(8)  << Electron_eta->At(j)
         << setw(8)  << Electron_phi->At(j)
         << setw(8)  << Electron_charge->At(j)
         << setw(8)  << Electron_dxy->At(j)
         << setw(8)  << Electron_dz->At(j)
         << setw(9)  << static_cast<int>(Electron_lostHits->At(j))
         << setw(11) << Electron_sip3d->At(j)
         << setw(7)  << (Electron_mvaFall17V2Iso_WP90->At(j) ? "T" : "F")
         << setw(7)  << (Electron_mvaFall17V2Iso_WP80->At(j) ? "T" : "F")
         << setw(9)  << Electron_pfRelIso03_chg->At(j)
         << setw(9)  << Electron_pfRelIso03_all->At(j)
         << endl;
    }
  }
  // Taus                                                 
  size_t nTauCand = *nTau->Get();
  if (nTauCand > 0) {
    os << " -- #Taus: " << nTauCand << endl;
    os << "  indx       pT      eta      phi charge      dxy      dz    DM   idDM  jetVeto   muVeto  eleVeto"
       << endl;
    for (size_t j = 0; j < nTauCand; ++j) {
      os << setw(6) << j
         << setw(9) << Tau_pt->At(j)
         << setw(9) << Tau_eta->At(j)
         << setw(9) << Tau_phi->At(j)
         << setw(7) << Tau_charge->At(j)
         << setw(9) << Tau_dxy->At(j)
         << setw(9) << Tau_dz->At(j)
         << setw(6) << Tau_idDecayMode->At(j)
	 << setw(6) << Tau_idDecayModeNewDMs->At(j)
         << setw(9) << (isSignal() ? (Tau_idDeepTau2017v2VSjet->At(j) >> 4 & 0x1) : (Tau_idDeepTau2017v2p1VSjet->At(j) >> 4 & 0x1))
         << setw(9) << (isSignal() ? (Tau_idDeepTau2017v2VSe->At(j)        & 0x1) : (Tau_idDeepTau2017v2p1VSe->At(j)        & 0x1))
         << setw(9) << (isSignal() ? (Tau_idDeepTau2017v2VSmu->At(j)       & 0x1) : (Tau_idDeepTau2017v2p1VSmu->At(j)       & 0x1))
         << endl;
    }
  }
  // Jets
  size_t nJetCand = *nJet->Get();
  if (nJetCand > 0) {
    os << " -- #Jets: " << nJetCand << endl;
    os << "  indx       pT      eta      phi  puId     mass  muIdx1  muIdx2  elIdx1  elIdx2  bTagCSVV2  bTagDFB"
       << endl;
    for (size_t j = 0; j < nJetCand; ++j) {
      os << setw(6)  << j
         << setw(9)  << Jet_nomPt->At(j)
         << setw(9)  << Jet_eta->At(j)
         << setw(9)  << Jet_phi->At(j)
         << setw(6)  << Jet_puId->At(j)
         << setw(9)  << Jet_mass->At(j)
         << setw(8)  << Jet_muIdx1->At(j)
         << setw(8)  << Jet_muIdx2->At(j)
         << setw(8)  << Jet_elIdx1->At(j)
         << setw(8)  << Jet_elIdx2->At(j)
         << setw(11) << Jet_btagCSVV2->At(j)
         << setw(9)  << Jet_btagDeepFlavB->At(j) 
         << endl;
    }
  }
  // FatJets
  size_t nFatJetCand = *nFatJet->Get();
  if (nFatJetCand > 0) {
    os << " -- #FatJets: " << nFatJetCand << endl;
    os << "  indx  jetId       pT      eta      phi  nsubjet  sjIdx1  sjIdx2   sdmass     tau1     tau2  bTagCSVV2  bTagDeepB"
       << endl;
    for (size_t j = 0; j < nFatJetCand; ++j) {
      os << setw(6)  << j
         << setw(7)  << FatJet_jetId->At(j)
         << setw(9)  << FatJet_pt->At(j)
         << setw(9)  << FatJet_eta->At(j)
         << setw(9)  << FatJet_phi->At(j)
         << setw(9)  << static_cast<int>(*nSubJet->Get()) 
         << setw(8)  << FatJet_subJetIdx1->At(j)
         << setw(8)  << FatJet_subJetIdx2->At(j)
         << setw(9)  << FatJet_msoftdrop->At(j)
         << setw(9)  << FatJet_tau1->At(j)
         << setw(9)  << FatJet_tau2->At(j)
         << setw(11) << FatJet_btagCSVV2->At(j)
         << setw(11) << FatJet_btagDeepB->At(j) 
         << endl;
    }
  }
  // FatJets
  size_t nSubJetCand = *nSubJet->Get();
  if (nFatJetCand > 0) {
    os << " -- #SubJets: " << nSubJetCand << endl;
    os << "  indx       pT      eta      phi     mass  bTagDeepB  rawFactor"
       << endl;
    for (size_t j = 0; j < nSubJetCand; ++j) {
      os << setw(6)  << j
         << setw(9)  << SubJet_pt->At(j)
         << setw(9)  << SubJet_eta->At(j)
         << setw(9)  << SubJet_phi->At(j)
         << setw(9)  << SubJet_mass->At(j)
         << setw(11) << SubJet_btagDeepB->At(j) 
         << setw(11) << SubJet_rawFactor->At(j) 
         << endl;
    }
  }
  // Add MET
  os << " -- Met information: " << endl;
  os << "       pT      phi      sig    sumEt" << endl;
  os << setw(9) << (isSignal() ? *Met_nomPt->Get()  : *Met_pt->Get())
     << setw(9) << (isSignal() ? *Met_nomPhi->Get() : *Met_phi->Get())
     << setw(9) << *Met_significance->Get()
     << setw(9) << *Met_sumEt->Get()
     << endl;
}
