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
  new TH1D("evtCutFlow", "Event CutFlow", 15, -0.5, 14.5);
  if (isMC()) new TH1D("evtCutFlowWt", "Event CutFlow (Weighted)", 15, -0.5, 14.5);
  new TH1D("SR_yield", "Yield in Signal Region", 9, -0.5, 8.5);
  if (isMC()) new TH1D("SR_yieldWt", "Yield in Signal Region Weighted", 9, -0.5, 8.5);
  new TH1D("FR_yield", "Yield in Fake Extrapolated Region", 9, -0.5, 8.5);
  if (isMC()) new TH1D("FR_yieldWt", "Yield in Fake Extrapolated Region Weighted", 9, -0.5, 8.5);
  if (isMC()) new TH1D("EventWtSum", "Event Weight Sum", 1, -0.5, 0.5);

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

    // Here one will only get to know if any of the <Type> triggers give true or false
    // If you want to know the HLT_path also, you can modify the isTriggered function in such a way
    // so that it returns the HLT_paths as well as the values which are giving true values
    bool passDoubleMuonHLT = AnaUtil::isTriggered(AnaBase::getDoubleMuonHLTpaths(),PhysicsObjSelector::getDoubleMuonHLTscores());
    bool passSingleEleHLT  = AnaUtil::isTriggered(AnaBase::getSingleElectronHLTpaths(),PhysicsObjSelector::getSingleElectronHLTscores());
    bool passDoubleEgHLT   = AnaUtil::isTriggered(AnaBase::getDoubleEgHLTpaths(),PhysicsObjSelector::getDoubleEgHLTscores());
    bool passSingleMuonHLT = AnaUtil::isTriggered(AnaBase::getSingleMuonHLTpaths(),PhysicsObjSelector::getSingleMuonHLTscores());
    bool passMuonEgHLT     = AnaUtil::isTriggered(AnaBase::getMuonEgHLTpaths(),PhysicsObjSelector::getMuonEgHLTscores());

    // Duplicate Event removal [FOR DATA ONLY]
    if (isDuplicate (passDoubleMuonHLT, passSingleEleHLT, passDoubleEgHLT, 
		     passSingleMuonHLT, passMuonEgHLT, getDatasetName())) continue;
    //Making the object collections ready!!!
    findObjects();

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
    if (preselMuColl.size() > 0) packLeptons <vhtm::Muon> (preselMuColl, preselLepColl, isMC());
    if (preselElColl.size() > 0) packLeptons <vhtm::Electron> (preselElColl, preselLepColl, isMC());
    std::sort(std::begin(preselLepColl), std::end(preselLepColl), PtComparator<LeptonCand>());

    std::vector<LeptonCand>fakeableLepColl;
    if (fakeableMuColl.size() > 0) packLeptons <vhtm::Muon> (fakeableMuColl, fakeableLepColl, isMC());
    if (fakeableElColl.size() > 0) packLeptons <vhtm::Electron> (fakeableElColl, fakeableLepColl, isMC());
    std::sort(std::begin(fakeableLepColl), std::end(fakeableLepColl), PtComparator<LeptonCand>()); //sorting fakeable lepton candidates

    std::vector<LeptonCand>tightLepColl;
    if (tightMuColl.size() > 0) packLeptons <vhtm::Muon> (tightMuColl, tightLepColl, isMC());
    if (tightElColl.size() > 0) packLeptons <vhtm::Electron> (tightElColl, tightLepColl, isMC());
    std::sort(std::begin(tightLepColl), std::end(tightLepColl), PtComparator<LeptonCand>()); //sorting fakeable lepton candidates

    if (fakeableLepColl.size() < 2) continue;
    AnaUtil::fillHist1DBasic("evtCutFlow", 2);
    AnaUtil::fillHist1DBasic("evtCutFlowWt", 2, MCweight*lumiFac, isMC());

    auto lep1 = fakeableLepColl[0];
    auto lep2 = fakeableLepColl[1];

    if (!(lep1.pt > 25 && lep2.pt > 20)) continue;
    AnaUtil::fillHist1DBasic("evtCutFlow", 3);
    AnaUtil::fillHist1DBasic("evtCutFlowWt", 3, MCweight*lumiFac, isMC());
    
    // ------------------------------- HLT Condition --------------------------------- //
    int  hltCounter   = 0;
    bool isMuMu   {false}; 
    bool isEleEle {false};
    bool isEleMu  {false};
    
    if (lep1.flavour == 1 && lep2.flavour == 1) {
      isMuMu = true;
      if (passDoubleMuonHLT || passSingleMuonHLT) hltCounter++;
    }
    else if (lep1.flavour == 2 && lep2.flavour == 2) {
      isEleEle = true;
      if (passDoubleEgHLT || passSingleEleHLT) hltCounter++;
    }
    else if ((lep1.flavour == 1 && lep2.flavour == 2) || (lep1.flavour == 2 && lep2.flavour == 1)) {
      isEleMu = true;
      if (passSingleMuonHLT || passSingleEleHLT || passMuonEgHLT) hltCounter++;
    }
    
    if (hltCounter == 0) continue;
    AnaUtil::fillHist1DBasic("evtCutFlow", 4);
    AnaUtil::fillHist1DBasic("evtCutFlowWt", 4, MCweight*lumiFac, isMC());
    // ------------------------------------------------------------------------------- //

    // Leptons must be of same sign
    if (lep1.charge * lep2.charge < 0) continue;
    AnaUtil::fillHist1DBasic("evtCutFlow", 5);
    AnaUtil::fillHist1DBasic("evtCutFlowWt", 5, MCweight*lumiFac, isMC());
    
    // no low mass resonance
    if (hasLowMassResonance(preselLepColl)) continue;
    AnaUtil::fillHist1DBasic("evtCutFlow", 6);
    AnaUtil::fillHist1DBasic("evtCutFlowWt", 6, MCweight*lumiFac, isMC());

    // no Z like candidate (Z->ll)
    if (hasZcandidate(preselLepColl)) continue;
    AnaUtil::fillHist1DBasic("evtCutFlow", 7);
    AnaUtil::fillHist1DBasic("evtCutFlowWt", 7, MCweight*lumiFac, isMC());

    // Max 2 tight leptons
    if (tightLepColl.size() > 2) continue;
    AnaUtil::fillHist1DBasic("evtCutFlow", 8);
    AnaUtil::fillHist1DBasic("evtCutFlowWt", 8, MCweight*lumiFac, isMC());

    // tau veto
    if (tauColl.size() > 0) continue;
    AnaUtil::fillHist1DBasic("evtCutFlow", 9);
    AnaUtil::fillHist1DBasic("evtCutFlowWt", 9, MCweight*lumiFac, isMC());

    bool hasOnlyPrompt = isPrompt(lep1) && isPrompt(lep2);
    bool hasNonPrompt  = (hasOnlyPrompt) ? false : true;
    //if (!(isPrompt(lep1) && isPrompt(lep2))) continue;
    AnaUtil::fillHist1DBasic("evtCutFlow", 10);
    AnaUtil::fillHist1DBasic("evtCutFlowWt", 10, MCweight*lumiFac, isMC());

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
    // nomatch : not yet genMatched : prompt or nonprompt -> Not decided yet
    size_t ntmu = 0;
    size_t ntel = 0;
    bool leadIsTightMuon {false};
    bool leadIsTightEle {false};

    for (size_t i = 0; i < 2; ++i) {
      auto indx = fakeableLepColl[i].index;
      auto flav = fakeableLepColl[i].flavour;
      if (flav == 1) {
	for (auto& mu : tightMuColl) 
	  if (indx == mu.index) {
	    ntmu++;
	    if (i == 0) leadIsTightMuon = true;
	  }
      }
      else if (flav == 2) {
	for (auto& el : tightElColl)
	  if (indx == el.index) {
	    ntel++;
	    if (i == 0) leadIsTightEle = true;
	  }
      }
    }
    
    bool isSR_nomatch = ((ntmu + ntel) == 2) ? true : false;
    bool isFR_nomatch = (isSR_nomatch) ? false : true;

    // ---------------------------------------------------------------------------------------------------------------- //
    // --------------------------------------- !!! Weight Factory and Flags !!! --------------------------------------- //
    // ---------------------------------------------------------------------------------------------------------------- //
    // For Signal region
    if (isMC() && isSR_nomatch) {
      if (isEleEle)     MCweight = MCweight 
			  * AnaBase::getIdSF("Tight", lep1.pt, lep1.SCeta, "Electron") 
			  * AnaBase::getIdSF("Tight", lep2.pt, lep2.SCeta, "Electron");
      else if (isMuMu)  MCweight = MCweight 
			  * AnaBase::getIdSF("Medium", lep1.pt, lep1.eta, "Muon") 
			  * AnaBase::getIdSF("Medium", lep2.pt, lep2.eta, "Muon")
			  * AnaBase::getIsoSF("Tight", lep1.pt, lep1.eta, "Muon")
			  * AnaBase::getIsoSF("Tight", lep2.pt, lep2.eta, "Muon");
      else if (isEleMu) {
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
    else if (isMC() && isFR_nomatch) {
      if (isEleEle) {
	if (ntel == 0) MCweight = - MCweight 
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
      else if (isMuMu) { 
	if (ntmu == 0) MCweight = - MCweight 
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
      else if (isEleMu) {
	if (ntel == 0 && ntmu == 0) MCweight = (lep1.flavour == 1) ? 
				    - MCweight 
				    * 0.01 
				    * 0.01 
				    * AnaBase::getIdSF("Loose", lep1.pt, lep1.eta, "Muon") 
				    * AnaBase::getIdSF("Loose", lep2.pt, lep2.SCeta, "Electron")
				    : MCweight
				    * 0.01 
				    * 0.01 
				    * AnaBase::getIdSF("Loose", lep1.pt, lep1.SCeta, "Electron") 
				    * AnaBase::getIdSF("Loose", lep2.pt, lep2.eta, "Muon");
	
	else if (ntel == 1 && ntmu == 0) 
	  MCweight = (leadIsTightEle) ? MCweight 
	    * 0.01 
	    * AnaBase::getIdSF("Tight", lep1.pt, lep1.SCeta, "Electron") 
	    * AnaBase::getIdSF("Loose", lep2.pt, lep2.eta, "Muon") 
	    : MCweight 
	    * 0.01 
	    * AnaBase::getIdSF("Loose", lep1.pt, lep1.eta, "Muon") 
	    * AnaBase::getIdSF("Tight", lep2.pt, lep2.SCeta, "Electron");
	else if (ntel == 0 && ntmu == 1) 
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
    // --------------------------------------------------------------------------------------------------------//
    //                                   Regions and Flags : Very Essential                                    //
    // --------------------------------------------------------------------------------------------------------//
    // isSR : Signal Region, isFR : Fake or Sideband Region, isNPFR : Non-Prompt and Fake Region
    bool isSR   = hasOnlyPrompt && isSR_nomatch; 
    bool isFR   = hasOnlyPrompt && isFR_nomatch;
    bool isNPFR = hasNonPrompt && isFR_nomatch;

    if (dumpIdx < 10) {
      std::cout<<"Event : "<<ev<<"\n";
      std::cout<<"ntel :"<<ntel<<"\t"
	       <<"ntmu :"<<ntmu<<"\t"
	       <<"lep1Flv :"<<lep1.flavour<<"\t"
	       <<"lep2Flv :"<<lep2.flavour<<"\t"
	       <<"leadIsTightMuon :"<<leadIsTightMuon<<"\t"
	       <<"leadIsTightEle :"<<leadIsTightEle<<"\n";
      std::cout<<"isSR :"<<isSR<<"\t"
	       <<"isFR :"<<isFR<<"\t"
	       <<"isNPFR :"<<isNPFR<<"\n"
	       <<"isEleEle :"<<isEleEle<<"\t"
	       <<"isEleMu  :"<<isEleMu<<"\t"
	       <<"isMuMu   :"<<isMuMu<<"\n"
	       <<"isSR_EleEle :"<<(isSR && isEleEle)<<"\t"
	       <<"isSR_EleMu  :"<<(isSR && isEleMu)<<"\t"
	       <<"isSR_MuMu   :"<<(isSR && isMuMu)<<"\n"
	       <<"isFR_EleEle :"<<(isFR && isEleEle)<<"\t"
	       <<"isFR_EleMu  :"<<(isFR && isEleMu)<<"\t"
	       <<"isFR_MuMu   :"<<(isFR && isMuMu)<<"\n"
	       <<"isNPFR_EleEle :"<<(isNPFR && isEleEle)<<"\t"
	       <<"isNPFR_EleMu  :"<<(isNPFR && isEleMu)<<"\t"
	       <<"isNPFR_MuMu   :"<<(isNPFR && isMuMu)<<"\n";
      std::cout<<setprecision(10)<<"MC_weight : "<<MCweight<<"\t"<<"Lumi_weight : "<<MCweight*lumiFac<<"\n"<<"\n";
    }
    dumpIdx++;

    // Channel Flags ... very important
    std::map <std::string, bool> channelFlags {}; 
    // Signal Region :: 
    // 1. both of the leading and sub-leading fakeable leptons are prompt (GenMatched)
    // 2. both of the fakeable leptons are tight
    // 3. channel name specified
    channelFlags.insert(std::make_pair("EleEle_SR", (isSR && isEleEle)));
    channelFlags.insert(std::make_pair("EleMu_SR",  (isSR && isEleMu)));
    channelFlags.insert(std::make_pair("MuMu_SR",   (isSR && isMuMu)));
    // Sideband region (/Fake application region) ::
    // 1. both of the leading and sub-leading fakeable leptons are prompt (GenMatched)
    // 2. at least one of the two leptons fail tightId
    // 3. channel name specified 
    // 4. not necessasry for signal samples
    channelFlags.insert(std::make_pair("EleEle_FR", (isFR && isEleEle && !isSignal())));
    channelFlags.insert(std::make_pair("EleMu_FR",  (isFR && isEleMu  && !isSignal())));
    channelFlags.insert(std::make_pair("MuMu_FR",   (isFR && isMuMu   && !isSignal())));
    // Fake Extrapolated region :: This will be estimated from data driven fake factor estimation
    // 1. At least one of the leading and sub-leading fakeable leptons is not prompt (GenMatched)
    // 2. at least one of the two leptons fail tightId 
    // 3. channel name specified 
    // 4. Only for MC : Non-Prompt & Fake Region for MC only
    // 5. not necessasry for signal samples  
    channelFlags.insert(std::make_pair("EleEle_NPFR", (isNPFR && isEleEle && isMC() && !isSignal())));
    channelFlags.insert(std::make_pair("EleMu_NPFR",  (isNPFR && isEleMu  && isMC() && !isSignal())));
    channelFlags.insert(std::make_pair("MuMu_NPFR",   (isNPFR && isMuMu   && isMC() && !isSignal())));

    // ----------------------------------------------------------------------------------------------------------------------- //
    // ----------------------------------------------------------------------------------------------------------------------- //
    // ----------------------------------------------------------------------------------------------------------------------- //

    // MET Distribution
    AnaUtil::fillHist1D("metPt", met.pt, 100, 0, 500, channelFlags, MCweight);
    if (met.pt < 60) continue;
    // From this point onwards, 
    // only prompt lepton contribution is goint to be considered for the evtCutFlow histogram
    AnaUtil::fillHist1DBasic("evtCutFlow", 11);
    AnaUtil::fillHist1DBasic("evtCutFlowWt", 11, MCweight*lumiFac, isMC());
    // nEvents in SR
    AnaUtil::fillHist1DBasic("evtCutFlow", 12, isSR);
    AnaUtil::fillHist1DBasic("evtCutFlowWt", 12, MCweight*lumiFac, isSR && isMC());

    // To get the yields for different regions
    AnaUtil::fillHist1DBasic("SR_yield", 0, (isSR && isEleEle));
    AnaUtil::fillHist1DBasic("SR_yieldWt", 0, MCweight*lumiFac, (isSR && isEleEle && isMC()));
    AnaUtil::fillHist1DBasic("SR_yield", 1, (isSR && isEleMu));
    AnaUtil::fillHist1DBasic("SR_yieldWt", 1, MCweight*lumiFac, (isSR && isEleMu && isMC()));
    AnaUtil::fillHist1DBasic("SR_yield", 2, (isSR && isMuMu));
    AnaUtil::fillHist1DBasic("SR_yieldWt", 2, MCweight*lumiFac, (isSR && isMuMu && isMC()));

    AnaUtil::fillHist1DBasic("FR_yield", 0, (isFR && isEleEle && !isSignal()));
    AnaUtil::fillHist1DBasic("FR_yieldWt", 0, MCweight*lumiFac, (isFR && isEleEle && isMC() && !isSignal()));
    AnaUtil::fillHist1DBasic("FR_yield", 1, (isFR && isEleMu && !isSignal()));
    AnaUtil::fillHist1DBasic("FR_yieldWt", 1, MCweight*lumiFac, (isFR && isEleMu && isMC() && !isSignal()));
    AnaUtil::fillHist1DBasic("FR_yield", 2, (isFR && isMuMu && !isSignal()));
    AnaUtil::fillHist1DBasic("FR_yieldWt", 2, MCweight*lumiFac, (isFR && isMuMu && isMC() && !isSignal()));

    bool isResolved_WZ = (jetColl.size() >= 3 && fatJetColl.size() == 0 && bJetColl.size()==0) ? true : false;
    bool isBoosted_WZ  = (fatJetColl.size() >= 1 && bJetColl.size()==0 && fatbJetColl.size()==0) ? true : false;

    // Need to save in skimmed tree
    // Channel EleEle : chTag = 1 | EleMu : chTag = 2 | MuMu : chTag = 3 
    float chTag = 0.0; //channel
    if (isSR) {
      if (isEleEle)      chTag = 1.0;
      else if (isEleMu)  chTag = 2.0;
      else if (isMuMu)   chTag = 3.0;
    }
  
    // ----------------------------- Resolved WZ Region ----------------------------- //
    if (isResolved_WZ) {
      AnaUtil::fillHist1DBasic("evtCutFlow",   13,                   isSR);
      AnaUtil::fillHist1DBasic("evtCutFlowWt", 13, MCweight*lumiFac, (isMC() && isSR));

      AnaUtil::fillHist1DBasic("SR_yield",      3,                   (isSR && isEleEle));
      AnaUtil::fillHist1DBasic("SR_yieldWt",    3, MCweight*lumiFac, (isMC() && isSR && isEleEle));
      AnaUtil::fillHist1DBasic("SR_yield",      4,                   (isSR && isEleMu));
      AnaUtil::fillHist1DBasic("SR_yieldWt",    4, MCweight*lumiFac, (isMC() && isSR && isEleMu));
      AnaUtil::fillHist1DBasic("SR_yield",      5,                   (isSR && isMuMu));
      AnaUtil::fillHist1DBasic("SR_yieldWt",    5, MCweight*lumiFac, (isMC() && isSR && isMuMu));

      AnaUtil::fillHist1DBasic("FR_yield",      3,                   (isFR && isEleEle && !isSignal()));
      AnaUtil::fillHist1DBasic("FR_yieldWt",    3, MCweight*lumiFac, (isMC() && isFR && isEleEle && !isSignal()));
      AnaUtil::fillHist1DBasic("FR_yield",      4,                   (isFR && isEleMu && !isSignal()));
      AnaUtil::fillHist1DBasic("FR_yieldWt",    4, MCweight*lumiFac, (isMC() && isFR && isEleMu && !isSignal()));
      AnaUtil::fillHist1DBasic("FR_yield",      5,                   (isFR && isMuMu && !isSignal()));
      AnaUtil::fillHist1DBasic("FR_yieldWt",    5, MCweight*lumiFac, (isMC() && isFR && isMuMu && !isSignal()));

      AnaUtil::fillHist1D("nAk4Jets_Resolved_WZ", jetColl.size(), 10, -0.5, 9.5, channelFlags, MCweight);
      AnaUtil::fillHist1D("ak4Jet1Pt_Resolved_WZ", jetColl[0].pt, 20, 0, 300, channelFlags, MCweight);
      AnaUtil::fillHist1D("ak4Jet2Pt_Resolved_WZ", jetColl[1].pt, 20, 0, 300, channelFlags, MCweight);
      AnaUtil::fillHist1D("ak4Jet3Pt_Resolved_WZ", jetColl[2].pt, 20, 0, 300, channelFlags, MCweight);
      AnaUtil::fillHist1D("ak4Jet4Pt_Resolved_WZ", jetColl[3].pt, 20, 0, 300, channelFlags, MCweight, jetColl.size() > 3);
      
      // jet inv mass
      float jetsInvM = (jetColl.size() > 3) ? (AnaUtil::getP4(jetColl[0]) + AnaUtil::getP4(jetColl[1]) + 
					       AnaUtil::getP4(jetColl[2]) + AnaUtil::getP4(jetColl[3])).M()
	: (AnaUtil::getP4(jetColl[0]) + AnaUtil::getP4(jetColl[1]) + AnaUtil::getP4(jetColl[2])).M();

      AnaUtil::fillHist1D("jetsInvMass_Resolved_WZ", jetsInvM, 100, 0, 500, channelFlags, MCweight);

      // --------------------------- Varibales to be plotted and stored in ntuple ------------------------- //      
      
      if (isSR && skimObj_) {
	TreeVariablesResolved varListR;

	varListR.MCweight               = MCweight;
	varListR.Channel                = chTag;
	varListR.pt_lep1                = lep1.pt;
	varListR.pt_lep2                = lep2.pt;
	varListR.invM_jets              = jetsInvM;
    
	skimObj_->fill(varListR, isResolved_WZ, isBoosted_WZ);    
      }
      histf()->cd(); // Very Very Very Essential
      // -------------------------------------------------------------------------------------------------- //

    }

    // ----------------------------- Boosted WZ Region ----------------------------- //
    if (isBoosted_WZ) {
      AnaUtil::fillHist1DBasic("evtCutFlow",   14,                   isSR);
      AnaUtil::fillHist1DBasic("evtCutFlowWt", 14, MCweight*lumiFac, (isMC() && isSR));

      AnaUtil::fillHist1DBasic("SR_yield",      6,                   (isSR && isEleEle));
      AnaUtil::fillHist1DBasic("SR_yieldWt",    6, MCweight*lumiFac, (isMC() && isSR && isEleEle));
      AnaUtil::fillHist1DBasic("SR_yield",      7,                   (isSR && isEleMu));
      AnaUtil::fillHist1DBasic("SR_yieldWt",    7, MCweight*lumiFac, (isMC() && isSR && isEleMu));
      AnaUtil::fillHist1DBasic("SR_yield",      8,                   (isSR && isMuMu));
      AnaUtil::fillHist1DBasic("SR_yieldWt",    8, MCweight*lumiFac, (isMC() && isSR && isMuMu));

      AnaUtil::fillHist1DBasic("FR_yield",      6,                   (isFR && isEleEle && !isSignal()));
      AnaUtil::fillHist1DBasic("FR_yieldWt",    6, MCweight*lumiFac, (isMC() && isFR && isEleEle && !isSignal()));
      AnaUtil::fillHist1DBasic("FR_yield",      7,                   (isFR && isEleMu && !isSignal()));
      AnaUtil::fillHist1DBasic("FR_yieldWt",    7, MCweight*lumiFac, (isMC() && isFR && isEleMu && !isSignal()));
      AnaUtil::fillHist1DBasic("FR_yield",      8,                   (isFR && isMuMu && !isSignal()));
      AnaUtil::fillHist1DBasic("FR_yieldWt",    8, MCweight*lumiFac, (isMC() && isFR && isMuMu && !isSignal()));

      AnaUtil::fillHist1D("nAk4Jets_Boosted_WZ", jetColl_ak8Cleaned.size(), 10, -0.5, 9.5,channelFlags, MCweight);
      AnaUtil::fillHist1D("ak4Jet1Pt_Boosted_WZ", jetColl_ak8Cleaned[0].pt, 50, 0, 1000, channelFlags, MCweight, jetColl_ak8Cleaned.size() >= 1);
      AnaUtil::fillHist1D("ak4Jet2Pt_Boosted_WZ", jetColl_ak8Cleaned[1].pt, 50, 0, 1000, channelFlags, MCweight, jetColl_ak8Cleaned.size() >= 2);
      AnaUtil::fillHist1D("nAk4Jets_has1FatJet_Boosted_WZ", jetColl_ak8Cleaned.size(), 10, -0.5, 9.5, channelFlags, MCweight, fatJetColl.size() == 1);
      AnaUtil::fillHist1D("nAk4Jets_has2OrMoreFatJet_Boosted_WZ", jetColl_ak8Cleaned.size(), 10, -0.5, 9.5, channelFlags, MCweight, fatJetColl.size() >= 2);

      // --------------------------- Varibales to be plotted and stored in ntuple ------------------------- //      
      if (isSR && skimObj_) {
	TreeVariablesBoosted  varListB;
	
	varListB.MCweight               = MCweight;
	varListB.Channel                = chTag;
	varListB.pt_lep1                = lep1.pt;
	varListB.pt_lep2                = lep2.pt;

	skimObj_->fill(varListB, isResolved_WZ, isBoosted_WZ);
      }
      histf()->cd(); // Very Very Very Essential
      // -------------------------------------------------------------------------------------------------- //
    } 
    if (!isMC()) selEvLog() << evt.run << " " << evt.lumis << " " << evt.event << std::endl;
  } // Event loop ends
}

bool MultiLeptonMVAna::isDuplicate(bool passDoubleMuonHLT, bool passDoubleEgHLT, bool passMuonEgHLT,
				   bool passSingleMuonHLT, bool passSingleEleHLT, std::string dataset) {
  if (isMC()) return false;
  else {
    std::map<std::string, bool>triggerOptions;
    std::map<std::string, bool>::iterator it;
    triggerOptions["DoubleMuon"]     =  passDoubleMuonHLT;
    triggerOptions["DoubleEG"]       =  passDoubleEgHLT && !passDoubleMuonHLT;
    triggerOptions["MuonEG"]         =  passMuonEgHLT && !(passDoubleMuonHLT || passDoubleEgHLT);
    triggerOptions["SingleMuon"]     =  passSingleMuonHLT && !(passDoubleMuonHLT || passDoubleEgHLT || passMuonEgHLT);
    triggerOptions["SingleElectron"] =  passSingleEleHLT && !(passDoubleMuonHLT || passDoubleEgHLT || passMuonEgHLT || passSingleMuonHLT);
    
    it = triggerOptions.find(dataset.c_str());
    if (it->second) return true;
  }
  return false;
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

bool MultiLeptonMVAna::isPrompt(LeptonCand lep) {
  if (!isMC()) return true;
  else {
    if (lep.genFlv == 1  || // prompt : final state 
	lep.genFlv == 15 || // from tau
	lep.genFlv == 22)   // from photon conversion (for ele only)
      return true;
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
      "isPrompt                            : ",
      "met > 60                            : ",
      "is SR                               : ",
      "isResolved_WZ                       : ",
      "isBoosted_WZ                        : "
      };
  vector<string> yieldLabels {
    "EleEle [All]",
      "EleMu  [All]",
      "MuMu   [All]",
      "EleEle [Resolved_WZ]",
      "EleMu  [Resolved_WZ]",
      "MuMu   [Resolved_WZ]",
      "EleEle [Boosted_WZ]",
      "EleMu  [Boosted_WZ]",
      "MuMu   [Boosted_WZ]"
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
