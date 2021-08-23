#include "configana.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <numeric>
#include <sstream>
#include <fstream>

#include "MultiLeptonMVAna.h"
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
using AnaUtil::bookHist1D;
using AnaUtil::fillHist1D;
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

  if (_createMVATree)
    skimObj_ = std::make_unique<MVASkim>(_mvaInputFile);
  else if (_readMVA)
    mvaObj_ = std::make_unique<MVAnalysis>(_MVAnetwork, _MVAxmlFile);

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
  bookHist1D("evtCutFlow", "Event CutFlow", 15, -0.5, 14.5);
  if (isMC()) {
    bookHist1D("EventWtSum", "Event Weight Sum", 1, -0.5, 0.5);
    bookHist1D("evtCutFlowWt", "Event CutFlow (Weighted)", 15, -0.5, 14.5);
  }
  bookHist1D("SR_yield", "Yield in Signal Region", 9, -0.5, 8.5);
  if (isMC()) bookHist1D("SR_yieldWt", "Yield in Signal Region Weighted", 9, -0.5, 8.5);
  bookHist1D("SB_yield", "Yield in Fake Extrapolated Region", 9, -0.5, 8.5);
  if (isMC()) bookHist1D("SB_yieldWt", "Yield in Fake Extrapolated Region Weighted", 9, -0.5, 8.5);


  std::vector<std::string>RegionFlags_ = {"SignalRegion","SidebandRegion"};
  std::vector<std::string>ChannelFlags_ = {"EleEle","EleMu", "MuMu"};

  bookHist1D("TestHist", "working on hist booking", 50, 0, 300, ChannelFlags_, RegionFlags_);
  bookHist1D("MetPt", "Missing E_{T} (GeV)", 50, 0, 300, ChannelFlags_, RegionFlags_);

  std::vector<std::string>ResolvedFlags_ = AnaUtil::mergeToList(RegionFlags_, "IsResolved");
  bookHist1D("NoAk4Jets", "nAk4 Jets", 10, 0, 10, ChannelFlags_, ResolvedFlags_);
  bookHist1D("Ak4Jet1Pt", "jet1 p_{T} (GeV)", 50, 0, 300, ChannelFlags_, ResolvedFlags_);
  bookHist1D("Ak4Jet2Pt", "jet2 p_{T} (GeV)", 50, 0, 300, ChannelFlags_, ResolvedFlags_);
  bookHist1D("Ak4Jet3Pt", "jet3 p_{T} (GeV)", 50, 0, 300, ChannelFlags_, ResolvedFlags_);
  bookHist1D("Ak4Jet4Pt", "jet4 p_{T} (GeV)", 50, 0, 300, ChannelFlags_, ResolvedFlags_);
  bookHist1D("JetsInvMass", "jets inv mass (GeV)", 100, 0, 500, ChannelFlags_, ResolvedFlags_);

  std::vector<std::string>BoostedFlags_ = AnaUtil::mergeToList(RegionFlags_, "IsBoosted");
  bookHist1D("NoAk4Jets", "No. of ak4 jets", 10, 0, 10, ChannelFlags_, BoostedFlags_);
  bookHist1D("Ak4Jet1Pt", "Leading ak4 jet p_{T} (GeV)", 50, 0, 300, ChannelFlags_, BoostedFlags_);
  bookHist1D("Ak4Jet2Pt", "Sub-leading ak4 jet p_{T} (GeV)", 50, 0, 300, ChannelFlags_, BoostedFlags_);
  bookHist1D("NoAk4JetsHas1FatJet", "No. of ak4 jets with one ak8", 10, 0, 10, ChannelFlags_, BoostedFlags_);
  bookHist1D("NoAk4JetsHas2orMoreFatJet", "No. of ak4 jets with ak8", 10, 0, 10, ChannelFlags_, BoostedFlags_);

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
    AnaUtil::fillHist1D("EventWtSum", 0, evtWeightSum_);
  }
  //--------------------------------------------------//
  //--------------------Event Loop--------------------//
  //--------------------------------------------------//
  int ev = 0;
  int dumpIdx = 0;
  while (treeReader()->Next()) {
    ev++;
    clearLists(); // reset analysis related lists for each event
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
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 0, MCweight*lumiFac);

    if (evt.nGoodPV < 1) continue;
    AnaUtil::fillHist1D("evtCutFlow", 1);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 1, MCweight*lumiFac);

    // Here one will only get to know if any of the <Type> triggers give true or false
    // If you want to know the HLT_path also, you can modify the isTriggered function in such a way
    // so that it returns the HLT_paths as well as the values which are giving true values
    bool passDoubleMuonHLT = AnaBase::isTriggered(getDoubleMuonHLTpaths(),  AnaBase::getHLTscores(getDoubleMuonHLTptrs()));
    bool passSingleEleHLT  = AnaBase::isTriggered(getSingleElectronHLTpaths(), AnaBase::getHLTscores(getSingleElectronHLTptrs()));
    bool passDoubleEgHLT   = AnaBase::isTriggered(getDoubleEgHLTpaths(), AnaBase::getHLTscores(getDoubleEgHLTptrs()));
    bool passSingleMuonHLT = AnaBase::isTriggered(getSingleMuonHLTpaths(), AnaBase::getHLTscores(getSingleMuonHLTptrs()));
    bool passMuonEgHLT     = AnaBase::isTriggered(getMuonEgHLTpaths(), AnaBase::getHLTscores(getMuonEgHLTptrs()));

    // Duplicate Event removal for data
    if (!isMC() && AnaBase::isDuplicate(passDoubleMuonHLT, passSingleEleHLT, passDoubleEgHLT, 
					passSingleMuonHLT, passMuonEgHLT, getDatasetName())) continue;
    // Prepare object collections
    findObjects();

    histf()->cd(); //required

    // Access Selected Objects
    const auto& preselMuColl       = getPreSelMuList();
    const auto& fakeableMuColl     = getFakeableMuList();
    const auto& tightMuColl        = getTightMuList();

    const auto& preselElColl       = getPreSelEleList();
    const auto& fakeableElColl     = getFakeableEleList();
    const auto& tightElColl        = getTightEleList();

    const auto& jetColl            = getCleanJetList();
    const auto& bJetColl           = getBJetList();
    const auto& jetColl_ak8Cleaned = getAk8CleanJetList();

    const auto& fatJetColl         = getCleanFatJetList();
    const auto& fatbJetColl        = getBTaggedFatJetList();

    const auto& tauColl            = getLepCleanTauList();
    const vhtm::MET& met           = getMETList().at(0);


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
    AnaUtil::fillHist1D("evtCutFlow", 2);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 2, MCweight*lumiFac);

    auto lep1 = fakeableLepColl[0];
    auto lep2 = fakeableLepColl[1];

    if (!(lep1.pt > 25 && lep2.pt > 20)) continue;
    AnaUtil::fillHist1D("evtCutFlow", 3);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 3, MCweight*lumiFac);
    
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
    AnaUtil::fillHist1D("evtCutFlow", 4);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 4, MCweight * lumiFac);

    // Leptons must be of same sign
    if (lep1.charge * lep2.charge < 0) continue;
    AnaUtil::fillHist1D("evtCutFlow", 5);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 5, MCweight * lumiFac);
    
    // no low mass resonance
    if (hasLowMassResonance(preselLepColl)) continue;
    AnaUtil::fillHist1D("evtCutFlow", 6);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 6, MCweight * lumiFac);

    // no Z like candidate (Z->ll)
    if (hasZcandidate(preselLepColl)) continue;
    AnaUtil::fillHist1D("evtCutFlow", 7);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 7, MCweight * lumiFac);

    // Max 2 tight leptons
    if (tightLepColl.size() > 2) continue;
    AnaUtil::fillHist1D("evtCutFlow", 8);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 8, MCweight * lumiFac);

    // tau veto
    if (tauColl.size() > 0) continue;
    AnaUtil::fillHist1D("evtCutFlow", 9);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 9, MCweight * lumiFac);

    bool hasOnlyPrompt = isMC() ? (isPrompt(lep1) && isPrompt(lep2)) : true;
    bool hasNonPrompt  = hasOnlyPrompt ? false : true;
    AnaUtil::fillHist1D("evtCutFlow", 10);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 10, MCweight * lumiFac);

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

    // beforeGenMatch : not yet genMatched : prompt or nonprompt -> Not decided yet
    size_t ntmu{0}, ntel{0};
    bool leadIsTightMuon{false}, leadIsTightEle {false};

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
    
    bool isSR_beforeGenMatch = ((ntmu + ntel) == 2) ? true : false;
    bool isSB_beforeGenMatch = (isSR_beforeGenMatch) ? false : true;

    // ---------------------------------------------------------------------------------------------------------------- //
    // --------------------------------------- !!! Weight Factory and Flags !!! --------------------------------------- //
    // ---------------------------------------------------------------------------------------------------------------- //
    // For Signal region
    if (isMC()) {
      if (isSR_beforeGenMatch) {
	if (isEleEle)     MCweight *= SFHandler().getIdSF("Tight", lep1.pt, lep1.SCeta, "Electron") 
			            * SFHandler().getIdSF("Tight", lep2.pt, lep2.SCeta, "Electron");
	else if (isMuMu)  MCweight *= SFHandler().getIdSF("Medium", lep1.pt, lep1.eta, "Muon") 
			            * SFHandler().getIdSF("Medium", lep2.pt, lep2.eta, "Muon")
			            * SFHandler().getIsoSF("Tight", lep1.pt, lep1.eta, "Muon")
			            * SFHandler().getIsoSF("Tight", lep2.pt, lep2.eta, "Muon");
	else if (isEleMu) {
	  if (leadIsTightMuon) MCweight *= SFHandler().getIdSF("Medium", lep1.pt, lep1.eta, "Muon")
				         * SFHandler().getIsoSF("Tight", lep1.pt, lep1.eta, "Muon")
				         * SFHandler().getIdSF("Tight", lep2.pt, lep2.SCeta, "Electron");
	  else MCweight *= SFHandler().getIdSF("Tight", lep1.pt, lep1.SCeta, "Electron") 
		         * SFHandler().getIdSF("Medium", lep2.pt, lep2.eta, "Muon");
	}
      }
      // For Sideband region
      else if (isSB_beforeGenMatch) {
	if (isEleEle) {
	  if (ntel == 0) MCweight = - MCweight
			   * 0.01
			   * 0.01
			   * SFHandler().getIdSF("Loose", lep1.pt, lep1.SCeta, "Electron")
			   * SFHandler().getIdSF("Loose", lep2.pt, lep2.SCeta, "Electron");
	  else MCweight = (leadIsTightEle) ?
		 MCweight
		 * 0.01
		 * SFHandler().getIdSF("Tight", lep1.pt, lep1.SCeta, "Electron")
		 * SFHandler().getIdSF("Loose", lep2.pt, lep2.SCeta, "Electron")
		 : MCweight
		 * 0.01
		 * SFHandler().getIdSF("Loose", lep1.pt, lep1.SCeta, "Electron")
		 * SFHandler().getIdSF("Tight", lep2.pt, lep2.SCeta, "Electron");
	}
	else if (isMuMu) {
	  if (ntmu == 0) MCweight = - MCweight
			   * 0.01
			   * 0.01
			   * SFHandler().getIdSF("Loose", lep1.pt, lep1.eta, "Muon")
			   * SFHandler().getIdSF("Loose", lep2.pt, lep2.eta, "Muon");
	  else MCweight = (leadIsTightMuon) ? 
		 MCweight
		 * 0.01
		 * SFHandler().getIdSF("Medium", lep1.pt, lep1.eta, "Muon")
		 * SFHandler().getIsoSF("Tight", lep1.pt, lep1.eta, "Muon")
		 * SFHandler().getIdSF("Loose", lep2.pt, lep2.eta, "Muon")
		 : MCweight
		 * 0.01
		 * SFHandler().getIdSF("Loose", lep1.pt, lep1.eta, "Muon")
		 * SFHandler().getIdSF("Medium", lep2.pt, lep2.eta, "Muon")
		 * SFHandler().getIdSF("Tight", lep2.pt, lep2.eta, "Muon");
	}
	else if (isEleMu) {
	  if (ntel == 0 && ntmu == 0) MCweight = (lep1.flavour == 1) ?
					- MCweight
					* 0.01
					* 0.01
					* SFHandler().getIdSF("Loose", lep1.pt, lep1.eta, "Muon")
					* SFHandler().getIdSF("Loose", lep2.pt, lep2.SCeta, "Electron")
					: MCweight
					* 0.01
					* 0.01
					* SFHandler().getIdSF("Loose", lep1.pt, lep1.SCeta, "Electron")
					* SFHandler().getIdSF("Loose", lep2.pt, lep2.eta, "Muon");
	  
	  else if (ntel == 1 && ntmu == 0)
	    MCweight = (leadIsTightEle) ? MCweight
	      * 0.01
	      * SFHandler().getIdSF("Tight", lep1.pt, lep1.SCeta, "Electron")
	      * SFHandler().getIdSF("Loose", lep2.pt, lep2.eta, "Muon")
	      : MCweight
	      * 0.01
	      * SFHandler().getIdSF("Loose", lep1.pt, lep1.eta, "Muon")
	      * SFHandler().getIdSF("Tight", lep2.pt, lep2.SCeta, "Electron");
	  else if (ntel == 0 && ntmu == 1)
	    MCweight = (leadIsTightMuon) ? MCweight
	      * 0.01
	      * SFHandler().getIdSF("Medium", lep1.pt, lep1.eta, "Muon")
	      * SFHandler().getIsoSF("Tight", lep1.pt, lep1.eta, "Muon")
	      * SFHandler().getIdSF("Loose", lep2.pt, lep2.SCeta, "Electron")
	      : MCweight
	      * 0.01
	      * SFHandler().getIdSF("Loose", lep1.pt, lep1.SCeta, "Electron")
	      * SFHandler().getIdSF("Medium", lep2.pt, lep2.eta, "Muon")
	      * SFHandler().getIsoSF("Tight", lep2.pt, lep2.eta, "Muon");
	}
      }
    }
    // --------------------------------------------------------------------------------------------------------//
    //                                   Regions and Flags : Very Essential                                    //
    // --------------------------------------------------------------------------------------------------------//
    //    SR : Both Leptons are prompt & tight          |     SB : Prompt but atleast one fails tight Id       //
    //                                                  |        DataSB has both prompt and non-prompt         //
    //                                                  |            FFx(DataSB-McSB) = FFxNPSB                //
    //                                                  |                        |                             //
    //---------------------------------------------------------------------------V-----------------------------//
    //     FR : Both tight but non-prompt               |   NPSB : Non-prompt and atleast one fails tight Id   //
    //     aka MCfake : Fake from MC simulation         |   aka MC-Closure region ( NPSB from MC x FF )        //
    //                  (using MC only)                 |                                                      //
    //             [data-driven : Fake Band] <----------|-------------- FFxNPSB : Data driven                  //
    //---------------------------------------------------------------------------------------------------------//
    bool isSR          = hasOnlyPrompt && isSR_beforeGenMatch;
    bool isSB          = hasOnlyPrompt && isSB_beforeGenMatch;
    bool isMCclosure   = hasNonPrompt && isSB_beforeGenMatch; 
    bool isMCfake      = hasNonPrompt && isSR_beforeGenMatch; 

    if (dumpIdx < 10) {
      cout << "------------------------------------------------------------" << endl;
      cout << "Event : " << ev << endl;
      cout << "  ntel  ntmu lep1Flv lep2Flv leadIsTightMuon leadIsTightEle" << endl;
      cout << setw(6)  << ntel
           << setw(6)  << ntmu
	   << setw(8)  << lep1.flavour
	   << setw(8)  << lep2.flavour
	   << setw(16) << leadIsTightMuon
	   << setw(15) << leadIsTightEle
	   << endl << endl;
      cout << "  isSR  isSB MCclosure  MCFake"
           << endl;
      cout << setw(6) << isSR 
	   << setw(6) << isSB
	   << setw(10) << isMCclosure
	   << setw(8) << isMCfake
	   << endl << endl;
      cout << "EleEle EleMu  MuMu"
	   << endl;
      cout << setw(6) << isEleEle
	   << setw(6) << isEleMu
	   << setw(6) << isMuMu
	   << endl << endl;
      if (isMC()) {
	cout << "   MC_weight Lumi_weight"
	     << endl;
	cout << setprecision(5) << setw(12) << MCweight
	     << setw(12) << MCweight * lumiFac 
	     << endl;
      }
    }
    dumpIdx++;

    // Region Flags
    // One can add other rehions also e.g. MC_fake region, MC_closure region etc
    std::map <std::string, bool> regionFlags {
        {"SignalRegion", isSR},
	{"SidebandRegion", (isSB && !isSignal())}
    };
    // Channel Flags ... very important
    std::map <std::string, bool> channelFlags { 
        {"EleEle", isEleEle},
	{"EleMu",  isEleMu},
	{"MuMu",   isMuMu}
    };

    // ----------------------------------------------------------------------------------------------------------------------- //
    // ----------------------------------------------------------------------------------------------------------------------- //
    // ----------------------------------------------------------------------------------------------------------------------- //

    AnaUtil::fillHist1D("TestHist", met.pt, channelFlags, regionFlags, MCweight);

    // MET Distribution
    if (met.pt < 40) continue;
    // From this point onwards, 
    //dumpEvent(evt.event);

    // only prompt lepton contribution is goint to be considered for the evtCutFlow histogram
    AnaUtil::fillHist1D("evtCutFlow", 11);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 11, MCweight * lumiFac);

    // nEvents in SR
    if (isSR) {
      AnaUtil::fillHist1D("evtCutFlow", 12);
      if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 12, MCweight * lumiFac);

      // To get the yields for different regions
      if (isEleEle) {
	AnaUtil::fillHist1D("SR_yield", 0);
	if (isMC()) AnaUtil::fillHist1D("SR_yieldWt", 0, MCweight * lumiFac);
      }
      else if (isEleMu) {
	AnaUtil::fillHist1D("SR_yield", 1);
	if (isMC()) AnaUtil::fillHist1D("SR_yieldWt", 1, MCweight*lumiFac);
      }
      else if (isMuMu) {
	AnaUtil::fillHist1D("SR_yield", 2);
	if (isMC()) AnaUtil::fillHist1D("SR_yieldWt", 2, MCweight * lumiFac);
      }
    }
    if (isSB && !isSignal()) {
      if (isEleEle) {
	AnaUtil::fillHist1D("SB_yield", 0);
        if (isMC()) AnaUtil::fillHist1D("SB_yieldWt", 0, MCweight * lumiFac);
      }
      else if (isEleMu) {
	AnaUtil::fillHist1D("SB_yield", 1);
	if (isMC()) AnaUtil::fillHist1D("SB_yieldWt", 1, MCweight * lumiFac);
      }
      else if (isMuMu) {
	AnaUtil::fillHist1D("SB_yield", 2);
	if (isMC()) AnaUtil::fillHist1D("SB_yieldWt", 2, MCweight * lumiFac);
      }
    }

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
      std::map <std::string, bool> ResolvedFlags = AnaUtil::combineMaps(regionFlags, {{"IsResolved", isResolved_WZ}});

      if (isSR) {
	AnaUtil::fillHist1D("evtCutFlow", 13);
	if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 13, MCweight*lumiFac);

	if (isEleEle) {
	  AnaUtil::fillHist1D("SR_yield", 3);
	  if (isMC()) AnaUtil::fillHist1D("SR_yieldWt", 3, MCweight*lumiFac);
	}
	if (isEleMu) {
	  AnaUtil::fillHist1D("SR_yield", 4);
	  if (isMC()) AnaUtil::fillHist1D("SR_yieldWt", 4, MCweight * lumiFac);
	}
	if (isMuMu) {
	  AnaUtil::fillHist1D("SR_yield", 5);
	  if (isMC()) AnaUtil::fillHist1D("SR_yieldWt", 5, MCweight * lumiFac);
	}
      }
      if (isSB  && !isSignal()) {
	if (isEleEle) {
	  AnaUtil::fillHist1D("SB_yield", 3);
	  if (isMC()) AnaUtil::fillHist1D("SB_yieldWt", 3, MCweight * lumiFac);
	}
	if (isEleMu) {
	  AnaUtil::fillHist1D("SB_yield", 4);
	  if (isMC()) AnaUtil::fillHist1D("SB_yieldWt", 4, MCweight * lumiFac);
	}
	if (isMuMu) {
	  AnaUtil::fillHist1D("SB_yield", 5);
	  if (isMC()) AnaUtil::fillHist1D("SB_yieldWt", 5, MCweight * lumiFac);
	}
      }

      AnaUtil::fillHist1D("MetPt", met.pt, channelFlags, ResolvedFlags, MCweight);
      AnaUtil::fillHist1D("NoAk4Jets", jetColl.size(), channelFlags, ResolvedFlags, MCweight);
      AnaUtil::fillHist1D("Ak4Jet1Pt", jetColl[0].pt, channelFlags, ResolvedFlags, MCweight);
      AnaUtil::fillHist1D("Ak4Jet2Pt", jetColl[1].pt, channelFlags, ResolvedFlags, MCweight);
      AnaUtil::fillHist1D("Ak4Jet3Pt", jetColl[2].pt, channelFlags, ResolvedFlags, MCweight);
      if (jetColl.size() > 3) AnaUtil::fillHist1D("Ak4Jet4Pt", jetColl[3].pt, channelFlags, ResolvedFlags, MCweight);
      
      // jet inv mass
      float jetsInvM = (jetColl.size() > 3) ? (AnaUtil::getP4(jetColl[0]) + AnaUtil::getP4(jetColl[1]) + 
					       AnaUtil::getP4(jetColl[2]) + AnaUtil::getP4(jetColl[3])).M()
	: (AnaUtil::getP4(jetColl[0]) + AnaUtil::getP4(jetColl[1]) + AnaUtil::getP4(jetColl[2])).M();

      AnaUtil::fillHist1D("JetsInvMass", jetsInvM, channelFlags, ResolvedFlags, MCweight);

      // --------------------------- Varibales to be plotted and stored in ntuple ------------------------- //      
      
      if (isSR && skimObj_) {
	TreeVariablesResolved varList;

	varList.MCweight               = MCweight;
	varList.Channel                = chTag;
	varList.pt_lep1                = lep1.pt;
	varList.pt_lep2                = lep2.pt;
	varList.invM_jets              = jetsInvM;
    
	skimObj_->fill(varList);    
      }
      histf()->cd(); // Very Very Very Essential
    }

    // ----------------------------- Boosted WZ Region ----------------------------- //
    if (isBoosted_WZ) {
      std::map <std::string, bool> BoostedFlags = AnaUtil::combineMaps(regionFlags, {{"IsBoosted", isBoosted_WZ}});
      if (isSR) {
	AnaUtil::fillHist1D("evtCutFlow", 14);
	if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 14, MCweight*lumiFac);

	if (isEleEle) {
	  AnaUtil::fillHist1D("SR_yield", 6);
	  if (isMC()) AnaUtil::fillHist1D("SR_yieldWt", 6, MCweight*lumiFac);
	}
	else if (isEleMu) {
	  AnaUtil::fillHist1D("SR_yield", 7);
	  if (isMC()) AnaUtil::fillHist1D("SR_yieldWt", 7, MCweight * lumiFac);
	}
	else if (isMuMu) {
	  AnaUtil::fillHist1D("SR_yield", 8);
	  if (isMC()) AnaUtil::fillHist1D("SR_yieldWt", 8, MCweight * lumiFac);
	}
      }
      if (isSB  && !isSignal()) {
	if (isEleEle) {
	  AnaUtil::fillHist1D("SB_yield", 6);
	  if (isMC()) AnaUtil::fillHist1D("SB_yieldWt", 6, MCweight * lumiFac);
	}
	else if (isEleMu) {
	  AnaUtil::fillHist1D("SB_yield", 7);
	  if (isMC()) AnaUtil::fillHist1D("SB_yieldWt", 7, MCweight * lumiFac);
	}
	else if (isMuMu) {
	  AnaUtil::fillHist1D("SB_yield", 8);
	  if (isMC()) AnaUtil::fillHist1D("SB_yieldWt", 8, MCweight * lumiFac);
	}
      }
      AnaUtil::fillHist1D("NoAk4Jets", jetColl_ak8Cleaned.size(), channelFlags, BoostedFlags, MCweight);
      if (jetColl_ak8Cleaned.size() >= 1) 
	AnaUtil::fillHist1D("Ak4Jet1Pt", jetColl_ak8Cleaned[0].pt, channelFlags, BoostedFlags, MCweight);
      if (jetColl_ak8Cleaned.size() >= 2) 
	AnaUtil::fillHist1D("Ak4Jet2Pt", jetColl_ak8Cleaned[1].pt, channelFlags, BoostedFlags, MCweight);
      if (fatJetColl.size() == 1) 
	AnaUtil::fillHist1D("NoAk4JetsHas1FatJet", jetColl_ak8Cleaned.size(), channelFlags, BoostedFlags, MCweight);
      if (fatJetColl.size() >= 2) 
	AnaUtil::fillHist1D("NoAk4JetsHas2orMoreFatJet", jetColl_ak8Cleaned.size(), channelFlags, BoostedFlags, MCweight);

      // --------------------------- Varibales to be plotted and stored in ntuple ------------------------- //      
      if (isSR && skimObj_) {
	TreeVariablesBoosted  varList;
	
	varList.MCweight               = MCweight;
	varList.Channel                = chTag;
	varList.pt_lep1                = lep1.pt;
	varList.pt_lep2                = lep2.pt;

	skimObj_->fill(varList);
      }
      histf()->cd(); // Very Very Very Essential
      // -------------------------------------------------------------------------------------------------- //
    } 
    if (!isMC()) selEvLog() << evt.run << " " << evt.lumis << " " << evt.event << std::endl;
  } // Event loop ends
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

bool MultiLeptonMVAna::isPrompt(const LeptonCand& lep) {
  if (lep.genFlv == 1  || // prompt : final state 
      lep.genFlv == 15 || // from tau
      lep.genFlv == 22)   // from photon conversion (for ele only)
    return true;
  return false;
}

void MultiLeptonMVAna::endJob() {
  PhysicsObjSelector::endJob();
  
  histf()->cd();
  vector<string> evLabels {
    "Events processed                                : ",
      "has GoodPV                                    : ",
      "at least 2 fakeable leptons                   : ",
      "lep1pt grthan 25 and lep2pt grthan 20         : ",
      "pass HLT                                      : ",
      "2 or more fakeable leptons                    : ",
      "low mass resonance veto                       : ",
      "Z mass resonance veto                         : ",
      "has max 2 tight leptons                       : ",
      "tau veto                                      : ",
      "isPrompt                                      : ",
      "met above 40                                  : ",
      "is SR                                         : ",
      "isResolved_WZ                                 : ",
      "isBoosted_WZ                                  : "
      };
  AnaUtil::SetEvtCutFlowBinLabels("evtCutFlow", evLabels);
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
  if (!isSignal()) {
    AnaUtil::showYield("SB_yield",   yieldLabels, "Fake Extrapolation in Signal Region [unweighted]", "Yield");
  }
  if (isMC()) {
    cout << endl
         << "evtWeightSum: " << setw(10) << setprecision(0) << evtWeightSum_ << endl
         << "      lumiWt: " << setw(10) << setprecision(5) << lumiFac
         << endl;
    //AnaUtil::showEfficiency("evtCutFlowWt", evLabels, "Event Selection (Weighted)", "Events");
    AnaUtil::showYield("SR_yieldWt",   yieldLabels, "Prompt Contribution in Signal Region [Lumi weighted]", "Yield");
    if (!isSignal()) {
      AnaUtil::showYield("SB_yieldWt",   yieldLabels, "Fake Extrapolation in Signal Region [Lumi weighted]", "Yield");
    }
  }
}

void MultiLeptonMVAna::closeFiles() {
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
bool MultiLeptonMVAna::readJob(const string& jobFile, int& nFiles)
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
void MultiLeptonMVAna::printJob(ostream& os) const
{
  AnaBase::printJob(os);
  os << endl;
  
  os << "   useEventList: " << std::boolalpha << useEventList_ << endl
     << " dumpEventCount: " << dumpEventCount_ << endl
     << "   dumpEventMax: " << dumpEventCount_ << endl;
}
