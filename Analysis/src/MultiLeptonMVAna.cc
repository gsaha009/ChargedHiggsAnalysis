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
#include <memory>

#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TRandom.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TH1K.h"
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include "configana.h"
#include "LeptonCand.h"
#include "MultiLeptonMVAna.h"

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

// -----------
// Constructor
// -----------
MultiLeptonMVAna::MultiLeptonMVAna()
  : PhysicsObjSelector()
   {}

// ----------
// Destructor
// ----------
MultiLeptonMVAna::~MultiLeptonMVAna() 
{
}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool MultiLeptonMVAna::beginJob() 
{ 
  if (!PhysicsObjSelector::beginJob()) return false;

  histf()->cd();
  histf()->mkdir("TMVAnalysis");
  fakehistf()->cd();
  fakehistf()->mkdir("TMVAnalysis");

  bookHistograms();

  if (_createMVATree) {
#ifdef TRUE_CPP14
    skimObj_ = std::make_unique <MVASkim> (_mvaInputFile);
#else
    skimObj_ = std::unique_ptr <MVASkim>(new MVASkim (_mvaInputFile));
#endif
    if (!skimObj_) return false;
  }

  else if (_readMVA) {
#ifdef TRUE_CPP14
    mvaObj_ = std::make_unique<MVAnalysis>(_MVAnetwork, _MVAxmlFile);
#else
    mvaObj_ = std::unique_ptr <MVAnalysis>(new MVAnalysis (_MVAnetwork, _MVAxmlFile));
#endif
    if (!mvaObj_) return false;
  }

#ifdef  SKIP_DUPLICATE_ALL
  eventIdStore_.clear();
#endif

  return true;
}

// ---------------
// Book histograms
// ---------------
void MultiLeptonMVAna::bookHistograms()
{
  PhysicsObjSelector::bookHistograms();
  histf()->cd();
  histf()->cd("TMVAnalysis");

  // book histograms to be filled at different stages
  new TH1D("nPV", "Number of vertices", 60, 0, 60);
  new TH1D("nGoodPV", "Number of Good vertices", 60, 0, 60);
  if (isMC()) {
    new TH1D("puWeight", "PU wight", 100, -20., 20.);
    new TH1D("evWeight", "Evt Wt", 100, -20., 20.);
    new TH1D("TotWt", "Total Wt", 100, -20., 20.);
    new TH1D("pu_evtWt", "PU*Event weight factor (MC events)", 20, -10., 10.);
    new TH1D("nLHEJets", "nJets in ME level", 10, -0.5, 9.5);
  }
  new TH1D("bTagWt", "CSSV2", 100, 0., 10.);
  new TH1D("evType", "0,1,2=Same Charged Leptons :: -1=Opposite charged Leptons", 7, -1.5, 5.5);
  
  new TH1D("evtCutFlow", "Event CutFlow", 4, -0.5, 3.5);
  if (isMC()) new TH1D("evtCutFlowWt", "Event CutFlow (Weighted)", 4, -0.5, 3.5);
  new TH1D("nMuons", "nTightIsoMuons", 10, -0.5, 9.5);
  new TH1D("nElectrons", "nTightIsoElectrons", 10, -0.5, 9.5);


  histf()->cd();
  histf()->ls();

  fakehistf()->cd();
  fakehistf()->cd("TMVAnalysis");
  new TH1D("nMuons", "nTightIsoMuons", 10, -0.5, 9.5);
  new TH1D("nElectrons", "nTightIsoElectrons", 10, -0.5, 9.5);

  fakehistf()->cd();
  fakehistf()->ls();
}

// -------------------------------
// Clear vectors before event loop
// -------------------------------

void MultiLeptonMVAna::clearLists() {
  PhysicsObjSelector::clear();
}
// -------------------
// The main event loop
// -------------------

void MultiLeptonMVAna::eventLoop()
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
  int ev = 0;

  //--------------------------------------------------//
  //---Sum of GenEvWtSum over all the trees in chain--//
  //--------------------------------------------------//

  if (isMC()) {
    while (treeReaderRun()->Next()) {
      evtWeightSum_ += getGenSumW();
    }
    lumiFac = lumiWt(evtWeightSum_, true);
    cout <<">>> lumiWt: "<<lumiFac<<" <<<<"<<endl;
  }
  //--------------------------------------------------//
  //--------------------Event Loop--------------------//
  //--------------------------------------------------//
  while (treeReader()->Next()) {
    ev++;
    if (ev > levt) {
      std::cout<<"Event Loop Finished!!!"<<std::endl;
      break;
    }

    clearEvent(); // reset tree variables 
    clearLists(); // reset analysis related lists for each event

    string currentFile(gSystem->BaseName(chain()->GetCurrentFile()->GetName()));
    if (!findEventInfo()) std::cout<<">>>Event Info Not Found :( :( :(\n";


    if (isMC() && readGenInfo())
      if (!findGenPartInfo()) 
	std::cout<<">>>Gen Info Not Found :( :( :(\n";

    const auto& eventlist = getEventList();
    const vhtm::Event& evt = eventlist.at(0);

    // Show the status
    if (ev%nPrint == 0 || firstEvent() > -1)
    //    if (ev%nPrint == 0)
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
      if (eventIdMap().find(mkey.str()) != eventIdMap().end()) continue;
    }

    histf()->cd(); //required
    histf()->cd("TMVAnalysis");
    
    //if (isMC()) AnaUtil::fillHist1D ("nLHEJets", evt.nLHEJets, 1.0);
    // For data or for MC without pileup
    float puWt   = 1.0; //for Data
    float evWt   = 1.0; //for Data
    float allWt  = 1.0; //for Data

    AnaUtil::fillHist1D("bTagWt", evt.btagWeight_CSVV2);

    if (isMC()) {
      evWt = evt.genEvWt;
      AnaUtil::fillHist1D ("evWeight", evWt, 1.0);
      if (usePUWt()) {
	puWt = evt.PUWeight;
	AnaUtil::fillHist1D("puWeight", puWt);      
      }
      allWt = evWt * puWt * evt.btagWeight_CSVV2 * lumiFac;
      AnaUtil::fillHist1D ("TotWt", allWt, 1.0);
    }
    if (op.verbose && ev == 1) 
      cout << " <<<eventNo# " << setw(8) << ev <<"\n"
	   << " eventWt# " << setw(8) << evWt <<"\n"
	   << " pileUpWt# " << setw(8) << puWt <<"\n"
	   << " lumiWt# " << setw(8) << lumiFac <<"\n"
	   << " bTagWt# " << setw(8) << evt.btagWeight_CSVV2 << "\n"
	   << " totalWt# " << setw(8) << allWt << " >>> "
	   << endl;

    //    else evtWeightSum_ += evWt;

    AnaUtil::fillHist1D("evType", evt.evType, allWt);

    AnaUtil::fillHist1D("evtCutFlow", 0);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 0, allWt);

    AnaUtil::fillHist1D("nPV", evt.nPV, allWt);
    AnaUtil::fillHist1D("nGoodPV", evt.nGoodPV, allWt);
    if (evt.nGoodPV < 1) continue;
    AnaUtil::fillHist1D("evtCutFlow", 1);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 1, allWt);

    bool trigDoubleMuonHLT {false};
    bool trigSingleMuonHLT {false};
    bool trigDoubleEgHLT {false};
    bool trigSingleEleHLT {false};
    bool trigMuonEgHLT {false};

    // HLT Selection
    auto doubleMuHLTpaths  = AnaBase::getDoubleMuonHLTpaths();
    auto doubleMuHLTscores = PhysicsObjSelector::getDoubleMuonHLTscores();
    for (size_t i = 0; i < doubleMuHLTpaths.size(); ++i){
      if (doubleMuHLTscores[i]) {
	//std::cout<<doubleMuHLTpaths[i]<<"\t"<<doubleMuHLTscores[i]<<"\n";
	trigDoubleMuonHLT = true;
      }
    }
    auto doubleEgHLTpaths  = AnaBase::getDoubleEgHLTpaths();
    auto doubleEgHLTscores = PhysicsObjSelector::getDoubleEgHLTscores();
    for (size_t i = 0; i < doubleEgHLTpaths.size(); ++i){
      if (doubleEgHLTscores[i]) {
	//std::cout<<doubleEgHLTpaths[i]<<"\t"<<doubleEgHLTscores[i]<<"\n";
	trigDoubleEgHLT = true;
      }
    }
    auto muonEgHLTpaths  = AnaBase::getMuonEgHLTpaths();
    auto muonEgHLTscores = PhysicsObjSelector::getMuonEgHLTscores();
    for (size_t i = 0; i < muonEgHLTpaths.size(); ++i){
      if (muonEgHLTscores[i]) {
	//std::cout<<muonEgHLTpaths[i]<<"\t"<<muonEgHLTscores[i]<<"\n";
	trigMuonEgHLT = true;
      }
    }
    auto singleMuHLTpaths  = AnaBase::getSingleMuonHLTpaths();
    auto singleMuHLTscores = PhysicsObjSelector::getSingleMuonHLTscores();
    for (size_t i = 0; i < singleMuHLTpaths.size(); ++i){
      if (singleMuHLTscores[i]) {
	//std::cout<<singleMuHLTpaths[i]<<"\t"<<singleMuHLTscores[i]<<"\n";
	trigSingleMuonHLT = true;
      }
    }
    auto singleEleHLTpaths  = AnaBase::getSingleElectronHLTpaths();
    auto singleEleHLTscores = PhysicsObjSelector::getSingleElectronHLTscores();
    for (size_t i = 0; i < singleEleHLTpaths.size(); ++i){
      if (singleEleHLTscores[i]) {
	//std::cout<<singleEleHLTpaths[i]<<"\t"<<singleEleHLTscores[i]<<"\n";
	trigSingleEleHLT = true;
      }
    }
    //std::cout<<trigDoubleMuonHLT<<"  "<<trigSingleMuonHLT<<"  "<<trigDoubleEgHLT<<"  "<<trigSingleEleHLT<<"  "<<trigMuonEgHLT<<"  "<<"\n";
    if (!(trigDoubleMuonHLT || trigSingleMuonHLT || trigDoubleEgHLT || trigSingleEleHLT || trigMuonEgHLT)) continue;
    AnaUtil::fillHist1D("evtCutFlow", 2);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 2, allWt);
    
    //Making the object collections ready!!!
    findObjects();

    //if (ev < 100) dumpEverything (ev, fLog());

    histf()->cd(); //required
    histf()->cd("TMVAnalysis");

    //Access Selected Objects
    const auto& elColl      = getTightIsoEleList();
    const auto& muColl      = getTightIsoMuList();
    const auto& jetColl     = getCleanJetList();
    const auto& bJetColl    = getBJetList();
    const auto& fatJetColl  = getCleanFatJetList();
    const auto& fatbJetColl = getBTaggedFatJetList();
    //const auto& tauColl     = getIsoTauList();
    const vhtm::MET& met    = getMETList().at(0);

    //P A C K I N G  L E P T O N S in LepCandList_
    std::vector<LeptonCand>lepColl;
    if (muColl.size() > 0) packLeptons<vhtm::Muon>(muColl, lepColl);
    if (elColl.size() > 0)  packLeptons<vhtm::Electron>(elColl, lepColl);
    AnaUtil::fillHist1D("nLeptons", lepColl.size(), allWt);

    std::sort(std::begin(lepColl), std::end(lepColl), PtComparator<LeptonCand>()); //sorting lepton candidates

    if (lepColl.size() !=2) continue;
    AnaUtil::fillHist1D("evtCutFlow", 3);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 3, allWt);


    if (!isMC()) selEvLog() << evt.run << " " << evt.lumis << " " << evt.event << std::endl;
   // Print only the first n events; n configurable
   //  if (isMC() && dumpEventCount_ > -1 && ++nEventSel >= dumpEventCount_) continue;
  }
  // Analysis over
}

bool MultiLeptonMVAna::hasZcandidate(const std::vector<LeptonCand>& lepColl, float wt) {
  bool hasZToLLtp {false};
  bool hasZToLLfp {false};
  for (size_t i = 0; i < lepColl.size(); ++i){
    auto& lep1 = lepColl.at(i);
    TLorentzVector l1p4 = AnaUtil::getP4(lep1);
    for (size_t j = i+1; j < lepColl.size(); ++j){
      auto& lep2 = lepColl.at(j);
      TLorentzVector l2p4 = AnaUtil::getP4(lep2);
      double lepInvM = (l1p4+l2p4).M();
      AnaUtil::fillHist1D("l1l2InvM_soCsoF_pc", lepInvM, wt);
      if(lep1.flavour == lep2.flavour) {
	if (lep1.charge * lep2.charge > 0.0) {
	  AnaUtil::fillHist1D("l1l2InvM_sCsF", lepInvM, wt);
	  if ((lepInvM > 84. && lepInvM < 98.) || lepInvM < 12.)  {
	    AnaUtil::fillHist1D ("ZcandProb", 0, wt);
	    hasZToLLtp = true;
	  }
	}
	else if (lep1.charge * lep2.charge < 0.0) {
	  AnaUtil::fillHist1D("l1l2InvM_oCsF", lepInvM, wt);
	  if ((lepInvM > 76. && lepInvM < 106.) || lepInvM < 12.) {
	    AnaUtil::fillHist1D ("ZcandProb", 1, wt);
	    hasZToLLfp = true;
	  }
	}
      }
      if (hasZToLLtp || hasZToLLfp) return true;
    }
  }
  return false;
}


void MultiLeptonMVAna::endJob() {
  PhysicsObjSelector::endJob();
  
  histf()->cd();
  histf()->cd("TMVAnalysis");
  vector<string> evLabels {
    "0) Events processed   : ",
      "1) has GoodPV       : ",
      "2) HLT pass         : ",
      "3) nLeptons >= 2    : "
      };
  AnaUtil::showEfficiency("evtCutFlow", evLabels, "Event Selection");  
  if (isMC()) {
    cout << endl
         << "evtWeightSum: " << setw(10) << setprecision(0) << evtWeightSum_ << endl
         << "      lumiWt: " << setw(10) << setprecision(5) << lumiFac
         << endl;
    AnaUtil::showEfficiency("evtCutFlowWt", evLabels, "Event Selection (Weighted)", "Events");
  }
  else AnaUtil::showEfficiency("evtCutFlow", evLabels, "Event Selection", "Events");
}

void MultiLeptonMVAna::closeFiles() {
  AnaBase::closeFiles();
  // Take care of local stuff first 
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
bool MultiLeptonMVAna::readJob(const string& jobFile, int& nFiles)
{
  if (!PhysicsObjSelector::readJob(jobFile, nFiles)) return false;
  
  // Open the same file containing the datacards again to read analysis specific cards
  ifstream fin(jobFile.c_str(), std::ios::in);    
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
    std::cout << line << std::endl;
    assert(tokens.size() > 1);
    const string& key   = tokens[0];
    const string& value = tokens[1];
    if (key == "useEventList")
      useEventList_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "skipDuplicate")
      skipDuplicate_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "dumpGenInfo")
      dumpGenInfo_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "inputEventFile")
      eventFilelist_.push_back(value.c_str());
    else if (key == "syncDumpFile")
      dumpFilename_ = value.c_str();
    else if (key == "dumpEventMax")
      dumpEventCount_ = std::stoi(value.c_str());
    else if (key == "selectPartons")
      selectPM_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "nMEPartons")
      nMEPartons_ = std::stoi(value.c_str());
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
  
  if (!dumpFilename_.empty()) {
    syncDumpf_.open(dumpFilename_, std::ios::out);
    if (!syncDumpf_) {
      cerr << "Output File: " << dumpFilename_ << " could not be opened!" << endl;
      return false;
    }
  }  

  printJob();
  return true;
}
void MultiLeptonMVAna::printJob(ostream& os) const
{
  AnaBase::printJob(os);
  os << endl;
  /*
  os << "   useEventList: " << std::boolalpha << useEventList_ << endl
     << "  skipDuplicate: " << std::boolalpha << skipDuplicate_ << endl
     << " dumpEventCount: " << dumpEventCount_ << endl
     << "   syncDumpFile: " << dumpFilename_ << endl
     << "   dumpEventMax: " << dumpEventCount_ << endl
     << "  selectPartons: " << std::boolalpha << selectPM_ << endl
     << "     nMEPartons: " << nMEPartons_ << endl;*/
}
