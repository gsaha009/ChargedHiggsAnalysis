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
#include "TNtuple.h"
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include "configana.h"
#include "LeptonCand.h"
#include "MultiLeptonFakeEstimation.h"

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
MultiLeptonFakeEstimation::MultiLeptonFakeEstimation()
  : PhysicsObjSelector()
   {}

// ----------
// Destructor
// ----------
MultiLeptonFakeEstimation::~MultiLeptonFakeEstimation() 
{
}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool MultiLeptonFakeEstimation::beginJob() 
{ 
  if (!PhysicsObjSelector::beginJob()) return false;

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
void MultiLeptonFakeEstimation::bookHistograms()
{
  PhysicsObjSelector::bookHistograms();
  histf()->cd();

  // book basic histograms to be filled at different stages
  new TH1D("evtCutFlow", "Event CutFlow", 8, -0.5, 7.5);
  if (isMC()) new TH1D("evtCutFlowWt", "Event CutFlow (Weighted)", 8, -0.5, 7.5);

  histf()->ls();
}

// -------------------------------
// Clear vectors before event loop
// -------------------------------

void MultiLeptonFakeEstimation::clearLists() {
  PhysicsObjSelector::clear();
}
// -------------------
// The main event loop
// -------------------

void MultiLeptonFakeEstimation::eventLoop()
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
  int dumpIdx = 0;

  //--------------------------------------------------//
  //---Sum of GenEvWtSum over all the trees in chain--//
  //--------------------------------------------------//

  if (isMC()) {
    while (treeReaderRun()->Next()) {
      evtWeightSum_ += getGenSumW();
    }
    lumiFac = lumiWt(evtWeightSum_, true);
    cout <<setprecision(5)<<">>> evtWeightSum: "<<evtWeightSum_<<"\t"<<">>> lumiWt: "<<lumiFac<<" <<<<"<<endl;
    histf()->cd();
    AnaUtil::fillHist1DBasic("EventWtSum", 0, evtWeightSum_);
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
    
    // For data or for MC without pileup
    float puWt   = 1.0; //for Data
    float evWt   = 1.0; //for Data
    float MCweight  = 1.0; //for Data

    if (isMC()) {
      evWt = evt.genEvWt;
      puWt = evt.PUWeight;
      MCweight = evWt * puWt;
    }
    if (op.verbose && ev == 1) 
      cout << " <<<eventNo# " << setw(8) << ev <<"\n"
	   << " eventWt# " << setw(8) << evWt <<"\n"
	   << " pileUpWt# " << setw(8) << puWt <<"\n"
	   << " lumiWt# " << setw(8) << lumiFac <<"\n"
	   << " bTagWt# " << setw(8) << evt.btagWeight_CSVV2 << "\n"
	   << " totalWt# " << setw(8) << MCweight << " >>> "
	   << endl;

    AnaUtil::fillHist1DBasic("evtCutFlow", 0);
    AnaUtil::fillHist1DBasic("evtCutFlowWt", 0, MCweight*lumiFac, isMC());

    if (evt.nGoodPV < 1) continue;
    AnaUtil::fillHist1DBasic("evtCutFlow", 1);
    AnaUtil::fillHist1DBasic("evtCutFlowWt", 1, MCweight*lumiFac, isMC());

    bool trigSingleMuonHLT {false};
    bool trigSingleEleHLT {false};

    auto singleMuHLTpaths  = AnaBase::getSingleMuonHLTForFakepaths();
    auto singleMuHLTscores = PhysicsObjSelector::getSingleMuonHLTForFakescores();
    for (size_t i = 0; i < singleMuHLTpaths.size(); ++i){
      if (singleMuHLTscores[i]) {
	trigSingleMuonHLT = true;
	break;
      }
    }
    auto singleEleHLTpaths  = AnaBase::getSingleElectronHLTForFakepaths();
    auto singleEleHLTscores = PhysicsObjSelector::getSingleElectronHLTForFakescores();
    for (size_t i = 0; i < singleEleHLTpaths.size(); ++i){
      if (singleEleHLTscores[i]) {
	trigSingleEleHLT = true;
	break;
      }
    }

    // Duplicate Event removal [FOR DATA]
    /*
    if (!isMC()) {
      bool HLT_Scores[2] = {trigSingleMuonHLT, trigSingleEleHLT};
      std::string dataset = getDatasetName();
      if (dataset == "DoubleMuon" && !HLT_Scores[0])
	continue;
      else if (dataset == "DoubleEG" && (!HLT_Scores[1] || HLT_Scores[0]))
	continue;
      else if (dataset == "MuonEG" && (!HLT_Scores[2] || (HLT_Scores[0] || HLT_Scores[1])))
	continue;
      else if (dataset == "SingleMuon" && (!HLT_Scores[3] || (HLT_Scores[0] || HLT_Scores[1] || HLT_Scores[2])))
	continue;
      else if (dataset == "SingleElectron" && (!HLT_Scores[4] || (HLT_Scores[0] || HLT_Scores[1] || HLT_Scores[2] || HLT_Scores[3])))
	continue;
    } 
    */
    //Making the object collections ready!!!
    findObjects();

    //if (ev < 100) dumpEverything (ev, fLog()); // need to fix this function !
    histf()->cd(); //required

    //Access Selected Objects
    const auto& preselElColl       = getPreSelEleList();
    const auto& preselMuColl       = getPreSelMuList();
    const auto& fakeableElColl     = getFakeableEleList();
    const auto& fakeableMuColl     = getFakeableMuList();
    const auto& tightElColl        = getTightEleList();
    const auto& tightMuColl        = getTightMuList();
    const auto& jetColl            = getCleanJetList();
    const auto& bJetColl           = getBJetList();
    const auto& fatJetColl         = getCleanFatJetList();
    const auto& fatbJetColl        = getBTaggedFatJetList();
    const auto& tauColl            = getLepCleanTauList();
    const vhtm::MET& met           = getMETList().at(0);


    if (fakeableElColl.size() != 1) continue;
    AnaUtil::fillHist1DBasic("evtCutFlow", 2);
    AnaUtil::fillHist1DBasic("evtCutFlowWt", 2, MCweight*lumiFac, isMC());

    if (fakeableMuColl.size() > 0) continue;
    AnaUtil::fillHist1DBasic("evtCutFlow", 3);
    AnaUtil::fillHist1DBasic("evtCutFlowWt", 3, MCweight*lumiFac, isMC());

    auto lep = fakeableElColl[0];

    if (!trigSingleEleHLT) continue;
    AnaUtil::fillHist1DBasic("evtCutFlow", 4);
    AnaUtil::fillHist1DBasic("evtCutFlowWt", 4, MCweight*lumiFac, isMC());

    if (!(lep.genFlv == 1  || lep.genFlv == 15 || lep.genFlv == 22))  continue;
    AnaUtil::fillHist1DBasic("evtCutFlow", 5);
    AnaUtil::fillHist1DBasic("evtCutFlowWt", 5, MCweight*lumiFac, isMC());

    if (met.pt > 20) continue;
    AnaUtil::fillHist1DBasic("evtCutFlow", 6);
    AnaUtil::fillHist1DBasic("evtCutFlowWt", 6, MCweight*lumiFac, isMC());

    TLorentzVector metp4;
    metp4.SetPtEtaPhiE(met.pt,0.0,met.phi,met.pt);

    double mT = Calculate_MT(AnaUtil::getP4(lep), metp4);
    if (mT > 20) continue;
    AnaUtil::fillHist1DBasic("evtCutFlow", 7);
    AnaUtil::fillHist1DBasic("evtCutFlowWt", 7, MCweight*lumiFac, isMC());

    std::map <std::string, bool> channelFlags {};
    channelFlags.insert(std::make_pair("Loose", true));
    channelFlags.insert(std::make_pair("Tight", tightElColl.size() == 1));

    AnaUtil::fillHist1D("pt", lep.pt, 100, 0, 300, "Control", channelFlags, 1.0);
    AnaUtil::fillHist1D("eta", lep.eta, 20, -3.2, 3.2, "Control", channelFlags, 1.0);

    if (!isMC()) selEvLog() << evt.run << " " << evt.lumis << " " << evt.event << std::endl;
  } // Event loop ends
}


bool MultiLeptonFakeEstimation::hasZcandidate(const std::vector<LeptonCand>& lepColl) {
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

bool MultiLeptonFakeEstimation::hasLowMassResonance(const std::vector<LeptonCand>& lepColl) {
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

bool MultiLeptonFakeEstimation::isPrompt(LeptonCand lep) {
  if (!isMC()) return true;
  else {
    if (lep.genFlv == 1  || // prompt : final state 
	lep.genFlv == 15 || // from tau
	lep.genFlv == 22)   // from photon conversion (for ele only)
      return true;
  }
  return false;
}

void MultiLeptonFakeEstimation::endJob() {
  PhysicsObjSelector::endJob();
  
  histf()->cd();
  //histf()->cd("Analysis");
  vector<string> evLabels {
    "Events processed                    : ",
      "has GoodPV                          : ",
      "only 1 fakeable electron            : ",
      "no muons                            : ",
      "pass HLT                            : ",
      "is prompt                           : ",
      "met < 20 GeV                        : ",
      "mT < 20 GeV                         : "
      };
  
  AnaUtil::showEfficiency("evtCutFlow", evLabels, "Event Selection [Unweighted]");  
  if (isMC()) {
    cout << endl
         << "evtWeightSum: " << setw(10) << setprecision(0) << evtWeightSum_ << endl
         << "      lumiWt: " << setw(10) << setprecision(5) << lumiFac
         << endl;
    AnaUtil::showEfficiency("evtCutFlowWt", evLabels, "Event Selection (Weighted)", "Events");
  }
}

void MultiLeptonFakeEstimation::closeFiles() {
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
bool MultiLeptonFakeEstimation::readJob(const string& jobFile, int& nFiles)
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
void MultiLeptonFakeEstimation::printJob(ostream& os) const
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
