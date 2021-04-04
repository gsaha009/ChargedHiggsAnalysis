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

  //histf()->cd();
  //histf()->mkdir("Analysis");
  //fakehistf()->cd();
  //fakehistf()->mkdir("Analysis");

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

  // book basic histograms to be filled at different stages
  new TH1D("evtCutFlow", "Event CutFlow", 12, -0.5, 11.5);
  if (isMC()) new TH1D("evtCutFlowWt", "Event CutFlow (Weighted)", 12, -0.5, 11.5);
  new TH1D("SR_yield", "Yield in Signal Region", 3, -0.5, 2.5);
  if (isMC()) new TH1D("SR_yieldWt", "Yield in Signal Region Weighted", 3, -0.5, 2.5);
  new TH1D("FR_yield", "Yield in Fake Extrapolated Region", 3, -0.5, 2.5);
  if (isMC()) new TH1D("FR_yieldWt", "Yield in Fake Extrapolated Region Weighted", 3, -0.5, 2.5);

  histf()->ls();
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
  int dumpIdx = 0;

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
    //histf()->cd("Analysis");
    
    // For data or for MC without pileup
    float puWt   = 1.0; //for Data
    float evWt   = 1.0; //for Data
    float MCweight  = 1.0; //for Data

    if (isMC()) {
      evWt = evt.genEvWt;
      //AnaUtil::fillHist1D ("evWeight", evWt, 1.0);
      puWt = evt.PUWeight;
      //AnaUtil::fillHist1D("puWeight", puWt, 1.0);      
      MCweight = evWt * puWt;
      //AnaUtil::fillHist1D ("TotWt", MCweight, 1.0);
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
    if (isMC()) AnaUtil::fillHist1DBasic("evtCutFlowWt", 0, MCweight*lumiFac);

    //AnaUtil::fillHist1D("nPV", evt.nPV, MCweight);
    //AnaUtil::fillHist1D("nGoodPV", evt.nGoodPV, MCweight);
    if (evt.nGoodPV < 1) continue;
    AnaUtil::fillHist1DBasic("evtCutFlow", 1);
    if (isMC()) AnaUtil::fillHist1DBasic("evtCutFlowWt", 1, MCweight*lumiFac);

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
	trigDoubleMuonHLT = true;
        break;
      }
    }
    auto doubleEgHLTpaths  = AnaBase::getDoubleEgHLTpaths();
    auto doubleEgHLTscores = PhysicsObjSelector::getDoubleEgHLTscores();
    for (size_t i = 0; i < doubleEgHLTpaths.size(); ++i){
      if (doubleEgHLTscores[i]) {
	trigDoubleEgHLT = true;
	break;
      }
    }
    auto muonEgHLTpaths  = AnaBase::getMuonEgHLTpaths();
    auto muonEgHLTscores = PhysicsObjSelector::getMuonEgHLTscores();
    for (size_t i = 0; i < muonEgHLTpaths.size(); ++i){
      if (muonEgHLTscores[i]) {
	trigMuonEgHLT = true;
	break;
      }
    }
    auto singleMuHLTpaths  = AnaBase::getSingleMuonHLTpaths();
    auto singleMuHLTscores = PhysicsObjSelector::getSingleMuonHLTscores();
    for (size_t i = 0; i < singleMuHLTpaths.size(); ++i){
      if (singleMuHLTscores[i]) {
	trigSingleMuonHLT = true;
	break;
      }
    }
    auto singleEleHLTpaths  = AnaBase::getSingleElectronHLTpaths();
    auto singleEleHLTscores = PhysicsObjSelector::getSingleElectronHLTscores();
    for (size_t i = 0; i < singleEleHLTpaths.size(); ++i){
      if (singleEleHLTscores[i]) {
	trigSingleEleHLT = true;
	break;
      }
    }
    // Duplicate Event removal [FOR DATA]
    // array if 5 HLT Scores
    if (!isMC()) {
      bool HLT_Scores[5] = {trigDoubleMuonHLT, trigDoubleEgHLT, trigMuonEgHLT, trigSingleMuonHLT, trigSingleEleHLT};
      std::string dataset = getDatasetName();
      if (dataset == "DoubleMuon" && !HLT_Scores[0]) continue;
      else if (dataset == "DoubleEG" && (!HLT_Scores[1] || HLT_Scores[0])) continue;
      else if (dataset == "MuonEG" && (!HLT_Scores[2] || (HLT_Scores[0] || HLT_Scores[1]))) continue;
      else if (dataset == "SingleMuon" && (!HLT_Scores[3] || (HLT_Scores[0] || HLT_Scores[1] || HLT_Scores[2]))) continue;
      else if (dataset == "SingleElectron" && (!HLT_Scores[4] || (HLT_Scores[0] || HLT_Scores[1] || HLT_Scores[2] || HLT_Scores[3]))) continue;
    }  
    //Making the object collections ready!!!
    findObjects();

    //if (ev < 100) dumpEverything (ev, fLog());
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
    const auto& jetColl_ak8Cleaned = getAk8CleanJetList();

    //P A C K I N G  L E P T O N S
    std::vector<LeptonCand>preselLepColl;
    if (preselMuColl.size() > 0) packLeptons <vhtm::Muon> (preselMuColl, preselLepColl);
    if (preselElColl.size() > 0) packLeptons <vhtm::Electron> (preselElColl, preselLepColl);
    std::sort(std::begin(preselLepColl), std::end(preselLepColl), PtComparator<LeptonCand>());

    std::vector<LeptonCand>fakeableLepColl;
    if (fakeableMuColl.size() > 0) packLeptons <vhtm::Muon> (fakeableMuColl, fakeableLepColl);
    if (fakeableElColl.size() > 0) packLeptons <vhtm::Electron> (fakeableElColl, fakeableLepColl);
    std::sort(std::begin(fakeableLepColl), std::end(fakeableLepColl), PtComparator<LeptonCand>()); //sorting fakeable lepton candidates

    std::vector<LeptonCand>tightLepColl;
    if (tightMuColl.size() > 0) packLeptons <vhtm::Muon> (tightMuColl, tightLepColl);
    if (tightElColl.size() > 0) packLeptons <vhtm::Electron> (tightElColl, tightLepColl);
    std::sort(std::begin(tightLepColl), std::end(tightLepColl), PtComparator<LeptonCand>()); //sorting fakeable lepton candidates

    if (fakeableLepColl.size() < 2) continue;
    AnaUtil::fillHist1DBasic("evtCutFlow", 2);
    if (isMC()) AnaUtil::fillHist1DBasic("evtCutFlowWt", 2, MCweight*lumiFac);

    if (!(fakeableLepColl[0].pt > 25 && fakeableLepColl[1].pt > 20)) continue;
    AnaUtil::fillHist1DBasic("evtCutFlow", 3);
    if (isMC()) AnaUtil::fillHist1DBasic("evtCutFlowWt", 3, MCweight*lumiFac);
    
    // HLT Condition
    bool ismumu {false}; 
    bool iselel {false};
    bool iselmu {false};
    
    if (fakeableLepColl[0].flavour == 1 && fakeableLepColl[1].flavour == 1) 
      {if (trigDoubleMuonHLT || trigSingleMuonHLT) ismumu = true;}
    
    else if (fakeableLepColl[0].flavour == 2 && fakeableLepColl[1].flavour == 2)  
      {if (trigDoubleEgHLT || trigSingleEleHLT) iselel = true;}

    else if ((fakeableLepColl[0].flavour == 1 && fakeableLepColl[1].flavour == 2) 
	     || (fakeableLepColl[0].flavour == 2 && fakeableLepColl[1].flavour == 1)) 
      {if (trigSingleMuonHLT || trigSingleEleHLT || trigMuonEgHLT) iselmu = true;}

    if (!(ismumu || iselel || iselmu)) continue;
    AnaUtil::fillHist1DBasic("evtCutFlow", 4);
    if (isMC()) AnaUtil::fillHist1DBasic("evtCutFlowWt", 4, MCweight*lumiFac);

    // Leptons must be of same sign
    if (fakeableLepColl[0].charge * fakeableLepColl[1].charge < 0) continue;
    AnaUtil::fillHist1DBasic("evtCutFlow", 5);
    if (isMC()) AnaUtil::fillHist1DBasic("evtCutFlowWt", 5, MCweight*lumiFac);
    
    // no low mass resonance
    if (hasLowMassResonance(preselLepColl)) continue;
    AnaUtil::fillHist1DBasic("evtCutFlow", 6);
    if (isMC()) AnaUtil::fillHist1DBasic("evtCutFlowWt", 6, MCweight*lumiFac);

    // no Z like candidate (Z->ll)
    if (hasZcandidate(preselLepColl)) continue;
    AnaUtil::fillHist1DBasic("evtCutFlow", 7);
    if (isMC()) AnaUtil::fillHist1DBasic("evtCutFlowWt", 7, MCweight*lumiFac);

    // Max 2 tight leptons
    if (tightLepColl.size() > 2) continue;
    AnaUtil::fillHist1DBasic("evtCutFlow", 8);
    if (isMC()) AnaUtil::fillHist1DBasic("evtCutFlowWt", 8, MCweight*lumiFac);

    // tau veto
    if (tauColl.size() > 0) continue;
    AnaUtil::fillHist1DBasic("evtCutFlow", 9);
    if (isMC()) AnaUtil::fillHist1DBasic("evtCutFlowWt", 9, MCweight*lumiFac);

    // -------------------------------------------------------------------- //
    // To select the SR and Fake Region. Two leading leptons must be tight.
    // If both of the fakeable leptons are electron and both pass the tight
    // Id, tightEleIDSF*tightEleIdSF will be applied. If one of the loose 
    // leptons is Muon & the other is Electron passing the tight Id, the SF
    // would be tightEleIDSF*tightMuIdSF.
    // For fake region, both of the leptons can fail the tight condition or
    // one of them fails. Thats why several bools have been taken to define
    // correct weight for events of every possible category.
    // -------------------------------------------------------------------- //

    bool isSR {false}; // to define main SR, where both of the fakeable leptons are tight
    bool isFR {false}; // to define FR, opposite of main SR
    size_t tmu = 0, fmu = 0;
    size_t tel = 0, fel = 0;
    bool leadIsTightMuon {false};
    bool leadIsTightEle {false};

    for (size_t i = 0; i < 2; ++i) {
      auto indx = fakeableLepColl[i].index;
      auto flav = fakeableLepColl[i].flavour;
      if (flav == 1) {
	fmu++;
	for (auto& mu : tightMuColl) 
	  if (indx == mu.index) {
	    tmu++;
	    if (i == 0) leadIsTightMuon = true;
	  }
      }
      else if (flav == 2) {
	fel++;
	for (auto& el : tightElColl)
	  if (indx == el.index) {
	    tel++;
	    if (i == 0) leadIsTightEle = true;
	  }
      }
    }
    
    isSR = ((tmu + tel) == 2) ? true : false;
    isFR = (isSR) ? false : true;

    auto lep1 = fakeableLepColl[0];
    auto lep2 = fakeableLepColl[1];

    // -------------------------------------------------------------------------- //
    // ---------------------- !!! Weight Factory !!! ---------------------------- //
    // -------------------------------------------------------------------------- //
    bool isSR_EleEle = (isSR && tel == 2) ? true : false;
    bool isSR_MuMu   = (isSR && tmu == 2) ? true : false;
    bool isSR_EleMu  = (isSR && tel == 1 && tmu == 1) ? true : false;

    bool isFR_EleEle = (isFR && fel == 2) ? true : false;
    bool isFR_MuMu   = (isFR && fmu == 2) ? true : false;
    bool isFR_EleMu  = (isFR && fel == 1 && fmu ==1) ? true : false;

    // For Signal region
    if (isMC() && isSR) {
      if (isSR_EleEle)     MCweight = MCweight 
			     * AnaBase::getIdSF("Tight", lep1.pt, lep1.SCeta, "Electron") 
			     * AnaBase::getIdSF("Tight", lep2.pt, lep2.SCeta, "Electron");
      else if (isSR_MuMu)  MCweight = MCweight 
			     * AnaBase::getIdSF("Medium", lep1.pt, lep1.eta, "Muon") 
			     * AnaBase::getIdSF("Medium", lep2.pt, lep2.eta, "Muon")
			     * AnaBase::getIsoSF("Tight", lep1.pt, lep1.eta, "Muon")
			     * AnaBase::getIsoSF("Tight", lep2.pt, lep2.eta, "Muon");
      else if (isSR_EleMu) {
	if (leadIsTightMuon) MCweight = MCweight 
			       * AnaBase::getIdSF("Medium", lep1.pt, lep1.eta, "Muon")
			       * AnaBase::getIsoSF("Tight", lep1.pt, lep1.eta, "Muon")
			       * AnaBase::getIdSF("Tight", lep2.pt, lep2.SCeta, "Electron");
	else MCweight = MCweight 
	       * AnaBase::getIdSF("Tight", lep1.pt, lep1.SCeta, "Electron") 
	       * AnaBase::getIdSF("Medium", lep2.pt, lep2.eta, "Muon");
      }
    }
    // For Fake extrapolated region
    else if (isMC() && isFR) {
      if (isFR_EleEle) {
	if (tel == 0) MCweight = MCweight 
			* 0.01 
			* 0.01 
			* AnaBase::getIdSF("Loose", lep1.pt, lep1.SCeta, "Electron") 
			* AnaBase::getIdSF("Loose", lep2.pt, lep2.SCeta, "Electron");
	else MCweight = (leadIsTightEle) ? 
	       MCweight 
	       * 0.01 
	       * AnaBase::getIdSF("Tight", lep1.pt, lep1.SCeta, "Electron") 
	       * AnaBase::getIdSF("Loose", lep2.pt, lep2.SCeta, "Electron") 
	       : MCweight 
	       * 0.01 
	       * AnaBase::getIdSF("Loose", lep1.pt, lep1.SCeta, "Electron") 
	       * AnaBase::getIdSF("Tight", lep2.pt, lep2.SCeta, "Electron"); 
      }
      else if (isFR_MuMu) { 
	if (tmu == 0) MCweight = MCweight 
			* 0.01 
			* 0.01 
			* AnaBase::getIdSF("Loose", lep1.pt, lep1.eta, "Muon") 
			* AnaBase::getIdSF("Loose", lep2.pt, lep2.eta, "Muon");
	else MCweight = (leadIsTightMuon) ? 
	       MCweight 
	       * 0.01 
	       * AnaBase::getIdSF("Medium", lep1.pt, lep1.eta, "Muon")
	       * AnaBase::getIsoSF("Tight", lep1.pt, lep1.eta, "Muon")
	       * AnaBase::getIdSF("Loose", lep2.pt, lep2.eta, "Muon")
	       : MCweight 
	       * 0.01 
	       * AnaBase::getIdSF("Loose", lep1.pt, lep1.eta, "Muon") 
	       * AnaBase::getIdSF("Medium", lep2.pt, lep2.eta, "Muon")
	       * AnaBase::getIdSF("Tight", lep2.pt, lep2.eta, "Muon");
      }
      else if (isFR_EleMu) {
	if (tel == 0 && tmu == 0) MCweight = (lep1.flavour == 1) ? 
				    MCweight 
				    * 0.01 
				    * 0.01 
				    * AnaBase::getIdSF("Loose", lep1.pt, lep1.eta, "Muon") 
				    * AnaBase::getIdSF("Loose", lep2.pt, lep2.SCeta, "Electron")
				    : MCweight
				    * 0.01 
				    * 0.01 
				    * AnaBase::getIdSF("Loose", lep1.pt, lep1.SCeta, "Electron") 
				    * AnaBase::getIdSF("Loose", lep2.pt, lep2.eta, "Muon");
	
	else if (tel == 1 && tmu == 0) 
	  MCweight = (leadIsTightEle) ? MCweight 
	    * 0.01 
	    * AnaBase::getIdSF("Tight", lep1.pt, lep1.SCeta, "Electron") 
	    * AnaBase::getIdSF("Loose", lep2.pt, lep2.eta, "Muon") 
	    : MCweight 
	    * 0.01 
	    * AnaBase::getIdSF("Loose", lep1.pt, lep1.eta, "Muon") 
	    * AnaBase::getIdSF("Tight", lep2.pt, lep2.SCeta, "Electron");
	else if (tel == 0 && tmu == 1) 
	  MCweight = (leadIsTightMuon) ? MCweight 
	    * 0.01 
	    * AnaBase::getIdSF("Medium", lep1.pt, lep1.eta, "Muon") 
	    * AnaBase::getIsoSF("Tight", lep1.pt, lep1.eta, "Muon") 
	    * AnaBase::getIdSF("Loose", lep2.pt, lep2.SCeta, "Electron")
	    : MCweight 
	    * 0.01 
	    * AnaBase::getIdSF("Loose", lep1.pt, lep1.SCeta, "Electron") 
	    * AnaBase::getIdSF("Medium", lep2.pt, lep2.eta, "Muon") 
	    * AnaBase::getIsoSF("Tight", lep2.pt, lep2.eta, "Muon"); 
      }
    }

    if (dumpIdx < 5) {
      dumpIdx++;
      std::cout<<"Event : "<<ev<<"\n";
      std::cout<<"fel :"<<fel<<"\t"
	       <<"fmu :"<<fmu<<"\t"
	       <<"tel :"<<tel<<"\t"
	       <<"tmu :"<<tmu<<"\t"
	       <<"lep1Flv :"<<lep1.flavour<<"\t"
	       <<"lep2Flv :"<<lep2.flavour<<"\t"
	       <<"leadIsTightMuon :"<<leadIsTightMuon<<"\t"
	       <<"leadIsTightEle :"<<leadIsTightEle<<"\n";
      std::cout<<"isSR :"<<isSR<<"\t"
	       <<"isFR :"<<isFR<<"\t"
	       <<"isSR_EleEle :"<<isSR_EleEle<<"\t"
	       <<"isSR_EleMu :"<<isSR_EleMu<<"\t"
	       <<"isSR_MuMu :"<<isSR_MuMu<<"\t"
	       <<"isFR_EleEle :"<<isFR_EleEle<<"\t"
	       <<"isFR_EleMu :"<<isFR_EleMu<<"\t"
	       <<"isFR_MuMu :"<<isFR_MuMu<<"\n";
      std::cout<<setprecision(10)<<"MC_weight : "<<MCweight<<"\t"<<"Lumi_weight : "<<MCweight*lumiFac<<"\n";
    }

    std::string channels[3] = {"EleEle", "EleMu", "MuMu"};
    bool SR_flags[3]        = {isSR_EleEle, isSR_EleMu, isSR_MuMu};
    bool FR_flags[3]        = {isFR_EleEle, isFR_EleMu, isFR_MuMu};
    // -------------------------------------------------------------------------- //
    // MET Distribution
    AnaUtil::fillHist1D("metPt", met.pt, 100, 0, 500, "SR", channels, MCweight, SR_flags);
    if (met.pt < 60) continue;
    AnaUtil::fillHist1DBasic("evtCutFlow", 10);
    if (isMC()) AnaUtil::fillHist1DBasic("evtCutFlowWt", 10, MCweight*lumiFac);
        
    if (isSR) {
      AnaUtil::fillHist1DBasic("evtCutFlow", 11);
      if (isMC()) AnaUtil::fillHist1DBasic("evtCutFlowWt", 11, MCweight*lumiFac);
    }

    // To get the yields for different regions
    if (isSR_EleEle){
      AnaUtil::fillHist1DBasic("SR_yield", 0);
      if (isMC()) AnaUtil::fillHist1DBasic("SR_yieldWt", 0, MCweight*lumiFac);
    }
    else if (isSR_EleMu){
      AnaUtil::fillHist1DBasic("SR_yield", 1);
      if (isMC()) AnaUtil::fillHist1DBasic("SR_yieldWt", 1, MCweight*lumiFac);
    }
    else if (isSR_MuMu){
      AnaUtil::fillHist1DBasic("SR_yield", 2);
      if (isMC()) AnaUtil::fillHist1DBasic("SR_yieldWt", 2, MCweight*lumiFac);
    }
    else if (isFR_EleEle){
      AnaUtil::fillHist1DBasic("FR_yield", 0);
      if (isMC()) AnaUtil::fillHist1DBasic("FR_yieldWt", 0, MCweight*lumiFac);
    }
    else if (isFR_EleMu){
      AnaUtil::fillHist1DBasic("FR_yield", 1);
      if (isMC()) AnaUtil::fillHist1DBasic("FR_yieldWt", 1, MCweight*lumiFac);
    }
    else if (isFR_MuMu){
      AnaUtil::fillHist1DBasic("FR_yield", 2);
      if (isMC()) AnaUtil::fillHist1DBasic("FR_yieldWt", 2, MCweight*lumiFac);
    }

    
    bool isResolved_WZ = (jetColl.size() >= 3 && fatJetColl.size() == 0 && bJetColl.size()==0) ? true : false;
    bool isBoosted_WZ  = (fatJetColl.size() >= 1 && bJetColl.size()==0 && fatbJetColl.size()==0) ? true : false;
  
    if (isResolved_WZ) {
      AnaUtil::fillHist1D("nAk4Jets_Resolved_WZ", jetColl.size(), 10, -0.5, 9.5, "SR", channels, MCweight, SR_flags);
      if (!isSignal()) AnaUtil::fillHist1D("nAk4Jets_Resolved_WZ", jetColl.size(), 10, -0.5, 9.5, "FR", channels, MCweight, FR_flags);
      AnaUtil::fillHist1D("ak4Jet1Pt_Resolved_WZ", jetColl[0].pt, 20, 0, 300, "SR", channels, MCweight, SR_flags);
      if (!isSignal()) AnaUtil::fillHist1D("ak4Jet1Pt_Resolved_WZ", jetColl[0].pt, 20, 0, 300, "FR", channels, MCweight, FR_flags);
      AnaUtil::fillHist1D("ak4Jet2Pt_Resolved_WZ", jetColl[1].pt, 20, 0, 300, "SR", channels, MCweight, SR_flags);
      if (!isSignal()) AnaUtil::fillHist1D("ak4Jet2Pt_Resolved_WZ", jetColl[1].pt, 20, 0, 300, "FR", channels, MCweight, FR_flags);
      AnaUtil::fillHist1D("ak4Jet3Pt_Resolved_WZ", jetColl[2].pt, 20, 0, 300, "SR", channels, MCweight, SR_flags);
      if (!isSignal()) AnaUtil::fillHist1D("ak4Jet3Pt_Resolved_WZ", jetColl[2].pt, 20, 0, 300, "FR", channels, MCweight, FR_flags);
      if (jetColl.size() > 3) {
	AnaUtil::fillHist1D("ak4Jet4Pt_Resolved_WZ", jetColl[3].pt, 20, 0, 300, "SR", channels, MCweight, SR_flags);
	if (!isSignal()) AnaUtil::fillHist1D("ak4Jet4Pt_Resolved_WZ", jetColl[3].pt, 20, 0, 300, "FR", channels, MCweight, FR_flags);
      }
      // jet inv mass
      float jetsInvM = (jetColl.size() > 3) ? (AnaUtil::getP4(jetColl[0]) + AnaUtil::getP4(jetColl[1]) + 
					       AnaUtil::getP4(jetColl[2]) + AnaUtil::getP4(jetColl[3])).M()
	: (AnaUtil::getP4(jetColl[0]) + AnaUtil::getP4(jetColl[1]) + AnaUtil::getP4(jetColl[2])).M();
      AnaUtil::fillHist1D("jetsInvMass_Resolved_WZ", jetsInvM, 100, 0, 500, "SR", channels, MCweight, SR_flags);
      if (!isSignal()) AnaUtil::fillHist1D("jetsInvMass_Resolved_WZ", jetsInvM, 100, 0, 500, "FR", channels, MCweight, FR_flags);
    }

    if (isBoosted_WZ) {
      AnaUtil::fillHist1D("nAk4Jets_Boosted_WZ", jetColl_ak8Cleaned.size(), 10, -0.5, 9.5, "SR", channels, MCweight, SR_flags);
      if (!isSignal()) AnaUtil::fillHist1D("nAk4Jets_Boosted_WZ", jetColl_ak8Cleaned.size(), 10, -0.5, 9.5, "FR", channels, MCweight, FR_flags);
      if (jetColl_ak8Cleaned.size() >= 1){ 
	AnaUtil::fillHist1D("ak4Jet1Pt_Boosted_WZ", jetColl_ak8Cleaned[0].pt, 50, 0, 1000, "SR", channels, MCweight, SR_flags);
	if (!isSignal()) AnaUtil::fillHist1D("ak4Jet1Pt_Boosted_WZ", jetColl_ak8Cleaned[0].pt, 50, 0, 1000, "FR", channels, MCweight, FR_flags);
	if (jetColl_ak8Cleaned.size() >= 2){ 
	  AnaUtil::fillHist1D("ak4Jet2Pt_Boosted_WZ", jetColl_ak8Cleaned[1].pt, 50, 0, 1000, "SR", channels, MCweight, SR_flags);
	  if (!isSignal()) AnaUtil::fillHist1D("ak4Jet2Pt_Boosted_WZ", jetColl_ak8Cleaned[1].pt, 50, 0, 1000, "FR", channels, MCweight, FR_flags);
	}
      }
      if (fatJetColl.size() == 1) {
	AnaUtil::fillHist1D("nAk4Jets_has1FatJet_Boosted_WZ", jetColl_ak8Cleaned.size(), 10, -0.5, 9.5, "SR", channels, MCweight, SR_flags);
	if (!isSignal()) AnaUtil::fillHist1D("nAk4Jets_has1FatJet_Boosted_WZ", jetColl_ak8Cleaned.size(), 10, -0.5, 9.5, "FR", channels, MCweight, FR_flags);
      }
      if (fatJetColl.size() >= 2) {
	AnaUtil::fillHist1D("nAk4Jets_has2OrMoreFatJet_Boosted_WZ", jetColl_ak8Cleaned.size(), 10, -0.5, 9.5, "SR", channels, MCweight, SR_flags);
	if (!isSignal()) AnaUtil::fillHist1D("nAk4Jets_has2OrMoreFatJet_Boosted_WZ", jetColl_ak8Cleaned.size(), 10, -0.5, 9.5, "FR", channels, MCweight, FR_flags);
      }
    } 
    if (!isMC()) selEvLog() << evt.run << " " << evt.lumis << " " << evt.event << std::endl;
   // Print only the first n events; n configurable
   //  if (isMC() && dumpEventCount_ > -1 && ++nEventSel >= dumpEventCount_) continue;
  }
  // Analysis over
}


bool MultiLeptonMVAna::hasZcandidate(const std::vector<LeptonCand>& lepColl) {
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

bool MultiLeptonMVAna::hasLowMassResonance(const std::vector<LeptonCand>& lepColl) {
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


void MultiLeptonMVAna::endJob() {
  PhysicsObjSelector::endJob();
  
  histf()->cd();
  //histf()->cd("Analysis");
  vector<string> evLabels {
      "Events processed                    : ",
      "has GoodPV                          : ",
      "at least 2 fakeable leptons         : ",
      "lep1pt > 25 and lep2pt > 20         : ",
      "pass HLT                            : ",
      "OS fakeable leptons                 : ",
      "low mass resonance veto             : ",
      "Z mass resonance veto               : ",
      "has max 2 tight leptons             : ",
      "tau veto                            : ",
      "met > 60                            : ",
      "is SR                               : "
      };
  vector<string> yieldLabels {
      "EleEle",
      "EleMu",
      "MuMu"
	};

  AnaUtil::showEfficiency("evtCutFlow", evLabels, "Event Selection [Unweighted]");  
  AnaUtil::showYield("SR_yield",   yieldLabels, "Prompt Contribution in Signal Region [unweighted]", "Yield");
  if (!isSignal()) AnaUtil::showYield("FR_yield",   yieldLabels, "Fake Extrapolation in Signal Region [unweighted]", "Yield");
  if (isMC()) {
    cout << endl
         << "evtWeightSum: " << setw(10) << setprecision(0) << evtWeightSum_ << endl
         << "      lumiWt: " << setw(10) << setprecision(5) << lumiFac
         << endl;
    //AnaUtil::showEfficiency("evtCutFlowWt", evLabels, "Event Selection (Weighted)", "Events");
    AnaUtil::showYield("SR_yieldWt",   yieldLabels, "Prompt Contribution in Signal Region [Lumi weighted]", "Yield");
    if (!isSignal()) AnaUtil::showYield("FR_yieldWt",   yieldLabels, "Fake Extrapolation in Signal Region [Lumi weighted]", "Yield");
  }
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
