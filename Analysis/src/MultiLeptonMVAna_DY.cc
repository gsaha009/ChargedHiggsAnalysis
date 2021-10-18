#include "configana.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <numeric>
#include <sstream>
#include <fstream>

#include "MultiLeptonMVAna_DY.h"
#include "LeptonCand.h"

#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TRandom.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TH1K.h"

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
using AnaUtil::fillHist1D;
// -----------
// Constructor
// -----------
MultiLeptonMVAna_DY::MultiLeptonMVAna_DY()
  : PhysicsObjSelector()
   {}

// ----------
// Destructor
// ----------
MultiLeptonMVAna_DY::~MultiLeptonMVAna_DY() 
{
}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool MultiLeptonMVAna_DY::beginJob() 
{ 
  if (!PhysicsObjSelector::beginJob()) return false;

  bookHistograms();

  if (_createMVATree)
    skimObj_ = std::make_unique<MVASkim>(_mvaInputFile);
  else if (_readMVA)
    mvaObj_ = std::make_unique<MVAnalysis>(_MVAnetwork, _MVAxmlFile);

  return true;
}
// ---------------
// Book histograms
// ---------------
void MultiLeptonMVAna_DY::bookHistograms()
{
  PhysicsObjSelector::bookHistograms();
  HistBooker hbooker;
  histf()->cd();
  hbooker.bookHistograms_DY(isMC());
  histf()->ls();
}

// -------------------------------
// Clear vectors before event loop
// -------------------------------
void MultiLeptonMVAna_DY::clearLists() {
  PhysicsObjSelector::clear();
}
// -------------------
// The main event loop
// -------------------
void MultiLeptonMVAna_DY::eventLoop()
{
  Options op;
  op.verbose = true;
  op.usesbit = true;  // Crucial
  op.printselected = false;
  
  //if (!beginJob()) std::cout<<"See your beginJob() function!!!\n";
  int nPrint = std::max(100000, nEvents()/10000);
  
  // --------------------
  // Start the event loop
  // --------------------
  string lastFile;
  int fevt = (firstEvent() > -1) ? firstEvent() : 0;
  int levt = (lastEvent() > -1) ? lastEvent() : nEvents();
  cout << ">>> Event range: [" << fevt << ", " << levt -1 << "]" << endl;
  
  //--------------------------------------------------//
  //---Sum of GenEvWtSum over all the trees in chain--//
  //--------------------------------------------------//
  if (isMC()) {
    while (treeReaderRun()->Next()) {
      evtWeightSum_ += getGenSumW();
    }
    lumiFac = lumiWt(evtWeightSum_, true);
    cout <<setprecision(5)<< ">>> evtWeightSum: "<<evtWeightSum_
	 <<"\t"<<">>> lumiWt: "<<lumiFac
	 <<" <<<<"<<endl;
    histf()->cd();
    AnaUtil::fillHist1D("EventWtSum", 0, evtWeightSum_); // very important
  }
  //--------------------------------------------------//
  //--------------------Event Loop--------------------//
  //--------------------------------------------------//
  int ev = 0;
  while (treeReader()->Next()) {
    ev++;
    // reset analysis related lists for each event
    clearLists(); 
    if (!findEventInfo()) {
      cerr << ">>> Event Info Not Found! ev = " 
	   << ev << ", go to the next event!" 
	   << endl;
      continue;
    }
    
    if (isMC() && readGenInfo())
      if (!findGenPartInfo()) {
	cerr << ">>> Gen Info Not Found, ev = " 
	     << ev << ", go to the next event!" 
	     << endl;
	continue;
      }
    
    const vhtm::Event& evt = getEventList().at(0);
    
    string currentFile(gSystem->BaseName(chain()->GetCurrentFile()->GetName()));
    if (ev%nPrint == 0 || firstEvent() > -1)
      cout << "Tree# " << setw(4) << chain()->GetTreeNumber()
	   << " ==> " << currentFile
	   << " <<< Run# " << setw(8) << evt.run
	   << " Lumis# " << setw(6) << evt.lumis
	   << " Event# " << setw(12) << evt.event << " >>> "
	   << " Events proc. " << setw(8) << ((firstEvent() > -1) ? ev - firstEvent() : ev)
	   << endl;
    
    // Select a set of events by [run, event]
    if (useEventList_ && eventIdMap().size()) {
      std::ostringstream mkey;
      mkey << evt.run << "-" << evt.lumis << "-" << evt.event;
      if (eventIdMap().find(mkey.str()) == eventIdMap().end()) continue;
    }

    histf()->cd(); //required
    
    // For data or for MC without pileup
    float puWt   = 1.0;
    float evWt   = 1.0;
    if (isMC()) {
      evWt = evt.genEvWt;
      puWt = evt.PUWeight;
    }
    float MCweight = evWt * puWt;
    if (ev == 1) 
      cout << ">>> eventNo: "  << setw(8) << ev << endl
           << setprecision(3)
	   << "     eventWt: " << setw(8) << evWt << endl
	   << "    pileUpWt: " << setw(8) << puWt << endl
	   << "      lumiWt: " << setw(8) << lumiFac << endl
	   << "      bTagWt: " << setw(8) << evt.btagWeight_CSVV2 << endl
	   << "     totalWt: " << setw(8) << MCweight
	   << endl;
    
    AnaUtil::fillHist1D("evtCutFlow", 0);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 0, MCweight);
    
    bool passDoubleMuonHLT = AnaBase::isTriggered(getDoubleMuonHLTpaths(),  AnaBase::getHLTscores(getDoubleMuonHLTptrs()));
    bool passSingleEleHLT  = AnaBase::isTriggered(getSingleElectronHLTpaths(), AnaBase::getHLTscores(getSingleElectronHLTptrs()));
    bool passDoubleEgHLT   = AnaBase::isTriggered(getDoubleEgHLTpaths(), AnaBase::getHLTscores(getDoubleEgHLTptrs()));
    bool passSingleMuonHLT = AnaBase::isTriggered(getSingleMuonHLTpaths(), AnaBase::getHLTscores(getSingleMuonHLTptrs()));
    bool passMuonEgHLT     = AnaBase::isTriggered(getMuonEgHLTpaths(), AnaBase::getHLTscores(getMuonEgHLTptrs()));

    // Duplicate Event removal for data
    if (!isMC() && AnaBase::isDuplicate(passDoubleMuonHLT, 
					passDoubleEgHLT, 
					passMuonEgHLT, 
    					passSingleMuonHLT, 
					passSingleEleHLT, 
					getDatasetName(), 
					getEra())) continue;
    AnaUtil::fillHist1D("evtCutFlow", 1);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 1, MCweight);

    // Should have at-least one good PV
    if (evt.nGoodPV < 1) continue;
    AnaUtil::fillHist1D("evtCutFlow", 2);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 2, MCweight);

    // Prepare object collections
    findObjects();

    histf()->cd(); //required

    // Access Selected Objects
    const auto& fakeableMuColl     = getFakeableMuList();
    const auto& tightMuColl        = getTightMuList();
    const auto& fakeableElColl     = getFakeableEleList();
    const auto& tightElColl        = getTightEleList();
    const auto& jetColl            = getCleanJetList();
    const vhtm::MET& met           = getMETList().at(0);

    //P A C K I N G  L E P T O N S
    std::vector<LeptonCand>fakeableLepColl;
    if (fakeableMuColl.size() > 0) packLeptons <vhtm::Muon> (fakeableMuColl, fakeableLepColl, isMC());
    if (fakeableElColl.size() > 0) packLeptons <vhtm::Electron> (fakeableElColl, fakeableLepColl, isMC());
    std::sort(std::begin(fakeableLepColl), std::end(fakeableLepColl), PtComparator<LeptonCand>()); //sorting fakeable lepton candidates

    if (fakeableLepColl.size() < 2) continue;
    AnaUtil::fillHist1D("evtCutFlow", 3);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 3, MCweight);

    auto lep1 = fakeableLepColl[0];
    auto lep2 = fakeableLepColl[1];

    int  hltCounter   = 0;
    if (getEra() == 2018) {
      if (lep1.flavour == 1 && lep2.flavour == 1) {
	if (passDoubleMuonHLT || passSingleMuonHLT) hltCounter++;
      }
      else if (lep1.flavour == 2 && lep2.flavour == 2) {
	if (passDoubleEgHLT || passSingleEleHLT) hltCounter++;
      }
      else if ((lep1.flavour == 1 && lep2.flavour == 2) || (lep1.flavour == 2 && lep2.flavour == 1)) {
	if (passSingleMuonHLT || passSingleEleHLT || passMuonEgHLT) hltCounter++;
      }
      else std::cout<<"wrong flavour\n";
      if (hltCounter == 0) continue;
    }
    AnaUtil::fillHist1D("evtCutFlow", 4);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 4, MCweight);

    // looking for Z
    std::vector<ZSpace::ZCandidate>zcandList;
    if (tightMuColl.size() >= 2)  ZSelector(tightMuColl, zcandList);
    if (tightElColl.size() >= 2)  ZSelector(tightElColl, zcandList);
    std::sort(std::begin(zcandList), std::end(zcandList), MassDiffComparator<ZSpace::ZCandidate>());

    if(zcandList.size() == 0) continue;
    AnaUtil::fillHist1D("evtCutFlow", 5);                                                             
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 5, MCweight);

    auto& zcand = zcandList[0];
    if (zcand.flavour == 0)      AnaUtil::fillHist1D("MuMu_SR_DY_lep1lep2InvM",   zcand.mass, MCweight);
    else if (zcand.flavour == 1) AnaUtil::fillHist1D("EleEle_SR_DY_lep1lep2InvM", zcand.mass, MCweight);

    if(zcand.massDiff > 15) continue;
    AnaUtil::fillHist1D("evtCutFlow", 6);                                                             
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 6, MCweight);

    if (met.pt > 50) continue;
    AnaUtil::fillHist1D("evtCutFlow", 7);                                                             
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 7, MCweight);

    //std::cout<<"Z_flavour : "<<zcand.flavour<<"\n";
    bool isMuMu   = (zcand.flavour == 0) ? true : false;
    bool isEleEle = (zcand.flavour == 1) ? true : false;

    TLorentzVector lep1p4 = zcand.l1P4;
    TLorentzVector lep2p4 = zcand.l2P4;

    auto zlep1 = fakeableLepColl[zcand.l1Index];
    auto zlep2 = fakeableLepColl[zcand.l2Index];

    // ---------------------------------------------------------------------------------------------------------------- //
    // --------------------------------------- !!! Weight Factory and Flags !!! --------------------------------------- //
    // ---------------------------------------------------------------------------------------------------------------- //    
    if (isMC()) {
      if (isEleEle)
	MCweight *= SFHandler().getIdSF("Loose", zlep1.pt, zlep1.SCeta, "Electron") 
	  * SFHandler().getIdSF("Loose", zlep2.pt, zlep2.SCeta, "Electron");
      else if (isMuMu)  
	MCweight *= SFHandler().getIdSF("Loose", zlep1.pt, zlep1.eta, "Muon") 
	  * SFHandler().getIdSF("Loose", zlep2.pt, zlep2.eta, "Muon")
	  * SFHandler().getIsoSF("Tight", zlep1.pt, zlep1.eta, "Muon")
	  * SFHandler().getIsoSF("Tight", zlep2.pt, zlep2.eta, "Muon");
    }
    bool isSR = true;
    
    // Region Flags
    std::map <std::string, bool> regionFlags {
      {"_SR_", isSR}
    };
    // Channel Flags ... very important
    std::map <std::string, bool> channelFlags { 
      {"EleEle", isEleEle},
      {"MuMu",  isMuMu}
    };
    
    //dumpEvent(evt.event);
    float dr_l1l2   = lep1p4.DeltaR(lep2p4);
    float pt_l1l2   = (lep1p4+lep2p4).Pt();
    float dphi_l1l2 = TVector2::Phi_mpi_pi(lep1p4.Phi() - lep2p4.Phi());
    float deta_l1l2 = lep1p4.Eta() - lep2p4.Eta();
    float invM_l1l2 = (lep1p4+lep2p4).M();

    fillHist1D("Lep1pt",lep1p4.Pt(),channelFlags, regionFlags, MCweight);
    fillHist1D("Lep2pt",lep2p4.Pt(),channelFlags, regionFlags, MCweight);
    fillHist1D("DiLepPt", pt_l1l2,  channelFlags, regionFlags, MCweight);
    fillHist1D("Lep1eta",lep1p4.Eta(),channelFlags, regionFlags, MCweight);
    fillHist1D("Lep2eta",lep2p4.Eta(),channelFlags, regionFlags, MCweight);
    fillHist1D("Lep1phi",lep1p4.Phi(),channelFlags, regionFlags, MCweight);
    fillHist1D("Lep2phi",lep2p4.Phi(),channelFlags, regionFlags, MCweight);
    fillHist1D("DR_lep1lep2", dr_l1l2, channelFlags, regionFlags, MCweight);
    fillHist1D("DPhi_lep1lep2", dphi_l1l2, channelFlags, regionFlags, MCweight);
    fillHist1D("DEta_lep1lep2", deta_l1l2, channelFlags, regionFlags, MCweight);
    fillHist1D("InvM_l1l2", invM_l1l2, channelFlags, regionFlags, MCweight);

    fillHist1D("MetPt", met.pt, channelFlags, regionFlags, MCweight);
    fillHist1D("MetPhi", met.phi, channelFlags, regionFlags, MCweight);
    fillHist1D("NoAk4Jets", jetColl.size(), channelFlags, regionFlags, MCweight);

    for (size_t i=0; i < jetColl.size(); ++i) {
      auto& jet1 = jetColl[i];
      std::string jetPt_hname  = "Ak4Jet"+std::to_string(i+1)+"Pt";
      fillHist1D (jetPt_hname.c_str(), jet1.pt, channelFlags, regionFlags, MCweight);
      std::string jetEta_hname  = "Ak4Jet"+std::to_string(i+1)+"Eta";
      fillHist1D (jetEta_hname.c_str(), jet1.eta, channelFlags, regionFlags, MCweight);
      std::string jetPhi_hname  = "Ak4Jet"+std::to_string(i+1)+"Phi";
      fillHist1D (jetPhi_hname.c_str(), jet1.phi, channelFlags, regionFlags, MCweight);
    }

    if (!isMC()) selEvLog() << evt.run << " " << evt.lumis << " " << evt.event << std::endl;
  } 
  //Event loop ends
}


bool MultiLeptonMVAna_DY::hasZcandidate(const std::vector<LeptonCand>& lepColl) {
  bool hasZToLL {false};
  for (size_t i = 0; i < lepColl.size(); ++i){
    auto& lep1 = lepColl.at(i);
    TLorentzVector l1p4 = AnaUtil::getP4(lep1);
    for (size_t j = i+1; j < lepColl.size(); ++j){
      auto& lep2 = lepColl.at(j);
      TLorentzVector l2p4 = AnaUtil::getP4(lep2);
      double lepInvM = (l1p4+l2p4).M();
      if(lep1.flavour == lep2.flavour) {
	if (lep1.charge * lep2.charge < 0.0 && std::fabs(lepInvM - 91.1876) < 10.0) {
	  hasZToLL = true;
	  if (hasZToLL) break;
	}
      }
    }
    if (hasZToLL) return true;
  }
  return false;
}

bool MultiLeptonMVAna_DY::hasLowMassResonance(const std::vector<LeptonCand>& lepColl) {
  bool hasLowMassRes {false};
  for (size_t i = 0; i < lepColl.size(); ++i){
    auto& lep1 = lepColl.at(i);
    TLorentzVector l1p4 = AnaUtil::getP4(lep1);
    for (size_t j = i+1; j < lepColl.size(); ++j){
      auto& lep2 = lepColl.at(j);
      TLorentzVector l2p4 = AnaUtil::getP4(lep2);
      double lepInvM = (l1p4+l2p4).M();
      if (lepInvM < 12.0) {
	hasLowMassRes = true;
	if (hasLowMassRes) break;
      }
    }
    if (hasLowMassRes) return true;
  }
  return false;
}

bool MultiLeptonMVAna_DY::isPrompt(const LeptonCand& lep) {
  //std::cout<<lep.genFlv<<std::endl;
  if (lep.genFlv == 1  || // prompt : final state 
      lep.genFlv == 15) // from tau
    //lep.genFlv == 22)   // from photon conversion (for ele only)
    return true;
  return false;
}

void MultiLeptonMVAna_DY::endJob() {
  PhysicsObjSelector::endJob();
  
  histf()->cd();
  vector<string> evLabels {
    "Events processed",
      "dataMasking",
      "has GoodPV",
      "nLoose leptons >= 2",
      "pass HLT",
      //"nTight leptons >= 2",
      "nZ_cand > 0",
      "found Zmass +-15 GeV",
      "met < 50 GeV"
      };
  AnaUtil::SetEvtCutFlowBinLabels("evtCutFlow", evLabels);
  
  AnaUtil::showEfficiency("evtCutFlow", evLabels, "Event Selection [Unweighted]");  
  if (isMC()) {
    cout << endl
         << "evtWeightSum: " << setw(10) << setprecision(0) << evtWeightSum_ << endl
         << "      lumiWt: " << setw(10) << setprecision(5) << lumiFac
         << endl;
    AnaUtil::showEfficiency("evtCutFlowWt", evLabels, "Event Selection (Weighted)", "Events");
  }
}

void MultiLeptonMVAna_DY::closeFiles() {
  AnaBase::closeFiles();
  // Take care of local stuff
  if (skimObj_ != nullptr) skimObj_->close(); 
}

// -------------------------------------------------------------------------------
// Poor man's way of a datacard. Each line between the 'START' and 'END' tags
// is read in turn, split into words, where the first element is the 'key' and
// the rest the value(s). If more than one values are present they are 
// stored in a vector. No safety mechanism is in place. Any line with an unknown 
// key is skipped. Comments lines should start with either '#' or '//', preferably
// in the first column. Empty lines are skipped. The file containing the datacards 
// is passed as the only argument of the program, there is no default
// -------------------------------------------------------------------------------
bool MultiLeptonMVAna_DY::readJob(const string& jobFile, int& nFiles)
{
  if (!PhysicsObjSelector::readJob(jobFile, nFiles)) return false;
  
  // Open the same file containing the datacards again to read analysis specific cards
  std::ifstream fin(jobFile.c_str(), std::ios::in);    
  if (!fin) {
    cerr << "==> Input File: " << jobFile << " could not be opened!" << endl;
    return false;
  }
  
  eventFilelist_.clear();  

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
    assert(tokens.size() > 1);
    const string& key   = tokens[0];
    const string& value = tokens[1];
    if (key == "useEventList")
      useEventList_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "dumpGenInfo")
      dumpGenInfo_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "inputEventFile")
      eventFilelist_.push_back(value.c_str());
    else if (key == "dumpEventMax")
      dumpEventCount_ = std::stoi(value.c_str());
    else if (key == "readMVA")
      _readMVA = (atoi(value.c_str()) > 0) ? true : false;
    else if (key == "MVAnetwork")
      _MVAnetwork = value;
    else if (key == "MVAxmlFile")
      _MVAxmlFile = value;
    else if (key == "createMVATree")
      _createMVATree = (atoi(value.c_str()) > 0) ? true : false;
    else if (key == "mvaInputFile")
      _mvaInputFile = value;
  }
  // Close the file
  fin.close();
  
  printJob();

  return true;
}
void MultiLeptonMVAna_DY::printJob(ostream& os) const
{
  AnaBase::printJob(os);
  os << endl;
  
  os << "   useEventList: " << std::boolalpha << useEventList_ << endl
     << " dumpEventCount: " << dumpEventCount_ << endl
     << "   dumpEventMax: " << dumpEventCount_ << endl;
}
