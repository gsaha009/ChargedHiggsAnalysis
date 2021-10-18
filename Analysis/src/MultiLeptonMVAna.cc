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
  HistBooker hbooker;
  histf()->cd();
  hbooker.bookHistograms(isMC());
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
  op.verbose = (verbosity()) ? true : false;
  op.usesbit = true;  // Crucial
  op.printselected = false;
  
  int nPrint = std::max(100000, nEvents()/10000);  
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
    lumiFac = lumiWt(evtWeightSum_, op.verbose);
    cout <<setprecision(5)<< ">>> evtWeightSum: "<<evtWeightSum_
	 <<"\t"<<">>> lumiWt: "<<lumiFac
	 <<" <<<<"<<endl;
    histf()->cd();
    fillHist1D("EventWtSum", 0, evtWeightSum_);
  }
  //--------------------------------------------------//
  //--------------------Event Loop--------------------//
  //--------------------------------------------------//
  int ev = 0;
  int dumpIdx = 0;
  while (treeReader()->Next()) {
    ev++;
    clearLists(); // reset analysis related lists for each event

    // Get Event info
    if (!findEventInfo()) { 
      cerr << ">>> Event Info Not Found! ev = " 
	   << endl;
      continue;
    }

    // Get Event
    const vhtm::Event& evt = getEventList().at(0); 
    if (op.verbose) std::cout << " ---------------------------- Event No. " << evt.event << " ---------------------------- "<< std::endl;

    // Get GenParticle info and LHE info
    if (isMC() && readGenInfo()) {
      if (!findGenPartInfo()) {
	cerr << ">>> Gen Info Not Found, ev = " 
	     << evt.event << ", go to the next event!" 
	     << endl;
	continue;
      }
      if (!findLHEPartInfo()) {
	cerr << ">>> LHE Info Not Found, ev = " 
	     << evt.event << ", go to the next event!" 
	     << endl;
	continue;
      }
    }

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
    float totalWeight = evWt * puWt;

    if (ev == 1 || op.verbose) {
      cout << "Weights INFO ===> \n";
      cout << ">>> eventNo: "  << setw(8) << ev << endl
           << setprecision(3)
	   << "     eventWt: " << setw(8) << evWt << endl
	   << "    pileUpWt: " << setw(8) << puWt << endl
	   << "      lumiWt: " << setw(8) << lumiFac << endl
	   << "      bTagWt: " << setw(8) << evt.btagWeight_CSVV2 << endl
	   << "     totalWt: " << setw(8) << totalWeight
	   << endl;
    }
    
    fillHist1D("evtCutFlow", 0);
    if (isMC()) fillHist1D("evtCutFlowWt", 0, totalWeight);
    
    if (evt.nGoodPV < 1) continue;
    fillHist1D("evtCutFlow", 1);
    if (isMC()) fillHist1D("evtCutFlowWt", 1, totalWeight);

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
					passSingleMuonHLT, passMuonEgHLT, getDatasetName(), getEra(), op.verbose)) continue;

    fillHist1D("evtCutFlow", 2);
    if (isMC()) fillHist1D("evtCutFlowWt", 2, totalWeight);
    
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
    fillHist1D("evtCutFlow", 3);
    if (isMC()) fillHist1D("evtCutFlowWt", 3, totalWeight);

    auto lep1 = fakeableLepColl[0];
    auto lep2 = fakeableLepColl[1];

    if (!(lep1.pt > 25 && lep2.pt > 20)) continue;
    fillHist1D("evtCutFlow", 4);
    if (isMC()) fillHist1D("evtCutFlowWt", 4, totalWeight);
    
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
    else std::cout<<"wrong lepton flavour\n";
    
    if (hltCounter == 0) continue;
    fillHist1D("evtCutFlow", 5);
    if (isMC()) fillHist1D("evtCutFlowWt", 5, totalWeight);

    // Leptons must be of same sign
    if (lep1.charge * lep2.charge < 0) continue;
    fillHist1D("evtCutFlow", 6);
    if (isMC()) fillHist1D("evtCutFlowWt", 6, totalWeight);
    
    // no low mass resonance
    if (hasLowMassResonance(preselLepColl)) continue;
    fillHist1D("evtCutFlow", 7);
    if (isMC()) fillHist1D("evtCutFlowWt", 7, totalWeight);

    // no Z like candidate (Z->ll)
    if (hasZcandidate(preselLepColl)) continue;
    fillHist1D("evtCutFlow", 8);
    if (isMC()) fillHist1D("evtCutFlowWt", 8, totalWeight);

    // Max 2 tight leptons
    if (tightLepColl.size() > 2) continue;
    fillHist1D("evtCutFlow", 9);
    if (isMC()) fillHist1D("evtCutFlowWt", 9, totalWeight);

    // tau veto
    if (tauColl.size() > 0) continue;
    fillHist1D("evtCutFlow", 10);
    if (isMC()) fillHist1D("evtCutFlowWt", 10, totalWeight);


    // MET Distribution
    if (met.pt < 40) continue;
    fillHist1D("evtCutFlow", 11);
    if (isMC()) fillHist1D("evtCutFlowWt", 11, totalWeight);

    bool hasOnlyPrompt = isMC() ? (isPrompt(lep1) && isPrompt(lep2)) : true;
    bool hasNonPrompt  = hasOnlyPrompt ? false : true;

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
    if (op.verbose) std::cout << "Weights and ScaleFactors >>>------> \n";

    if (isSR_beforeGenMatch) {
      // both of the fakeable leptons are tight
      if (op.verbose) std::cout << " --- SR :: PromptNotCheckedYet --- \n";
      if (isEleEle) { // both leptons are electrons
	double lep1_IdSF = (isMC()) ? SFHandler().getIdSF("Tight", lep1.pt, lep1.SCeta, "Electron") : 1.0;
	double lep2_IdSF = (isMC()) ? SFHandler().getIdSF("Tight", lep2.pt, lep2.SCeta, "Electron") : 1.0;
	totalWeight *= lep1_IdSF * lep2_IdSF;
	if (op.verbose) std::cout << "EleEle\n"
				  << "lep1_IdSF"   << setw(6) << lep1_IdSF << "\n"
				  << "lep2_IdSF"   << setw(6) << lep2_IdSF << "\n"
				  << "totalWeight" << setw(6) << totalWeight << "\n";
      }
      else if (isMuMu) { // both leptons are muons
	double lep1_IdSF  = (isMC()) ? SFHandler().getIdSF("Medium", lep1.pt, lep1.eta, "Muon") : 1.0;
	double lep2_IdSF  = (isMC()) ? SFHandler().getIdSF("Medium", lep2.pt, lep2.eta, "Muon") : 1.0;
	double lep1_IsoSF = (isMC()) ? SFHandler().getIsoSF("Tight", lep1.pt, lep1.eta, "Muon") : 1.0;
	double lep2_IsoSF = (isMC()) ? SFHandler().getIsoSF("Tight", lep2.pt, lep2.eta, "Muon") : 1.0; 
	totalWeight *= lep1_IdSF * lep2_IdSF * lep1_IsoSF * lep2_IsoSF;
	if (op.verbose) std::cout << "MuMu\n"
				  << "lep1_IdSF"   << setw(6) << lep1_IdSF << "\n"
				  << "lep2_IdSF"   << setw(6) << lep2_IdSF << "\n"
				  << "lep1_IsoSF"  << setw(6) << lep1_IsoSF << "\n"
				  << "lep2_IsoSF"  << setw(6) << lep2_IsoSF << "\n"
				  << "totalWeight" << setw(6) << totalWeight << "\n";	
      }
      else if (isEleMu) { // one is electron and other is muon
	if (leadIsTightMuon) { // leading lepton is tight Muon
	  double lep1_IdSF  = (isMC()) ? SFHandler().getIdSF("Medium", lep1.pt, lep1.eta, "Muon") : 1.0;
	  double lep1_IsoSF = (isMC()) ? SFHandler().getIsoSF("Tight", lep1.pt, lep1.eta, "Muon") : 1.0;
	  double lep2_IdSF  = (isMC()) ? SFHandler().getIdSF("Tight", lep2.pt, lep2.SCeta, "Electron") : 1.0;
	  totalWeight *= lep1_IdSF * lep1_IsoSF * lep2_IdSF;
	  if (op.verbose) std::cout << "MuEle\n"
				    << "lep1_IdSF"   << setw(6) << lep1_IdSF << "\n"
				    << "lep2_IdSF"   << setw(6) << lep2_IdSF << "\n"
				    << "lep1_IsoSF"  << setw(6) << lep1_IsoSF << "\n"
				    << "totalWeight" << setw(6) << totalWeight << "\n";
	}
	else if (leadIsTightEle) { // leading lepton is tight Electron
	  double lep1_IdSF  = (isMC()) ? SFHandler().getIdSF("Tight", lep1.pt, lep1.SCeta, "Electron") : 1.0;
	  double lep2_IdSF  = (isMC()) ? SFHandler().getIdSF("Medium", lep2.pt, lep2.eta, "Muon") : 1.0;
	  double lep2_IsoSF = (isMC()) ? SFHandler().getIsoSF("Tight", lep2.pt, lep2.eta, "Muon") : 1.0;
	  totalWeight *= lep1_IdSF * lep2_IdSF * lep2_IsoSF;
	  if (op.verbose) std::cout << "EleMu\n"
				    << "lep1_IdSF"   << setw(6) << lep1_IdSF   << "\n"
				    << "lep2_IdSF"   << setw(6) << lep2_IdSF   << "\n"
				    << "lep2_IsoSF"  << setw(6) << lep2_IsoSF  << "\n"
				    << "totalWeight" << setw(6) << totalWeight << "\n";
	}
	else {
	  std::cerr << "In SR : {isEleMu}, no tight lepton is found"
		    << std::endl;
	}
      }
    }
    
    // At least one fakeable lepton fails the tight criteria 
    else if (isSB_beforeGenMatch) {
      if (op.verbose) std::cout << " --- SB :: PromptNotCheckedYet --- \n";
      if (isEleEle) { // both leptons are electrons
	if (op.verbose) std::cout << " --- EleEle --- \n";
	if (ntel == 0) { // both electrons fail tight id
	  double lep1_FF = SFHandler().getFF(lep1.pt, lep1.eta,"Electron");
	  double lep2_FF = SFHandler().getFF(lep2.pt, lep2.eta,"Electron");
	  double lep1_IdSF = (isMC()) ? SFHandler().getIdSF("Loose", lep1.pt, lep1.SCeta, "Electron") : 1.0;
	  double lep2_IdSF = (isMC()) ? SFHandler().getIdSF("Loose", lep2.pt, lep2.SCeta, "Electron") : 1.0;
	  totalWeight *= -1.0 * lep1_FF * lep2_FF * lep1_IdSF * lep2_IdSF;
	  if (op.verbose) std::cout << "EleEle\n"
				    << "Both fail tight Id \n"
				    << "lep1_IdSF"   << setw(6) << lep1_IdSF << "\n"
				    << "lep2_IdSF"   << setw(6) << lep2_IdSF << "\n"
				    << "lep1_FF "    << setw(6) << lep1_FF   << "\n"
				    << "lep2_FF "    << setw(6) << lep2_FF   << "\n"
				    << "totalWeight" << setw(6) << totalWeight << "\n";
	}
	
	else { // only one electron fails tight id
	  if (leadIsTightEle) { // leading fakeable electron is tight and the sub-leading electron fails tight id
	    double lep1_IdSF = (isMC()) ? SFHandler().getIdSF("Tight", lep1.pt, lep1.SCeta, "Electron") : 1.0;
	    double lep2_IdSF = (isMC()) ? SFHandler().getIdSF("Loose", lep2.pt, lep2.SCeta, "Electron") : 1.0;
	    double lep2_FF   = SFHandler().getFF(lep2.pt, lep2.eta, "Electron");
	    totalWeight *= lep1_IdSF * lep2_IdSF * lep2_FF;
	    if (op.verbose) std::cout << "EleEle\n"
				      << "lead electron is tight | sub-lead electron fails tight \n"
				      << "lep1_IdSF"   << setw(6) << lep1_IdSF << "\n"
				      << "lep2_IdSF"   << setw(6) << lep2_IdSF << "\n"
				      << "lep2_FF "    << setw(6) << lep2_FF   << "\n"
				      << "totalWeight" << setw(6) << totalWeight << "\n";
	  }
	  else { // leading electron is the fakeable one and the sub-leading is the tight electron
	    double lep1_FF   = SFHandler().getFF(lep1.pt, lep1.eta, "Electron");
	    double lep1_IdSF = (isMC()) ? SFHandler().getIdSF("Loose", lep1.pt, lep1.SCeta, "Electron") : 1.0;
	    double lep2_IdSF = (isMC()) ? SFHandler().getIdSF("Tight", lep2.pt, lep2.SCeta, "Electron") : 1.0;
	    totalWeight *= lep1_FF * lep1_IdSF * lep2_IdSF;
	    if (op.verbose) std::cout << "EleEle\n"
				      << "lead electron fails tight | sub-lead electron is tight \n"
				      << "lep1_IdSF"   << setw(6) << lep1_IdSF << "\n"
				      << "lep2_IdSF"   << setw(6) << lep2_IdSF << "\n"
				      << "lep1_FF "    << setw(6) << lep1_FF   << "\n"
				      << "totalWeight" << setw(6) << totalWeight << "\n";
	  }
	}
      }

      else if (isMuMu) { // both leptons are muons
	if (op.verbose) std::cout << " --- MuMu --- \n";
	if (ntmu == 0) { // both muons fail tight id
	  double lep1_FF = SFHandler().getFF(lep1.pt, std::fabs(lep1.eta),"Muon");
	  double lep2_FF = SFHandler().getFF(lep2.pt, std::fabs(lep2.eta),"Muon");
	  double lep1_IdSF  = (isMC()) ? SFHandler().getIdSF("Loose", lep1.pt, lep1.eta, "Muon") : 1.0;
	  double lep2_IdSF  = (isMC()) ? SFHandler().getIdSF("Loose", lep2.pt, lep2.eta, "Muon") : 1.0;
	  double lep1_IsoSF = (isMC()) ? SFHandler().getIsoSF("Tight", lep1.pt, lep1.eta, "Muon") : 1.0;
	  double lep2_IsoSF = (isMC()) ? SFHandler().getIsoSF("Tight", lep2.pt, lep2.eta, "Muon") : 1.0;
	  totalWeight = -1.0 * lep1_FF * lep2_FF * lep1_IdSF * lep2_IdSF * lep1_IsoSF * lep2_IsoSF;
	  if (op.verbose) std::cout << "MuMu\n"
				    << "Both fail tight Id \n"
				    << "lep1_IdSF"   << setw(6) << lep1_IdSF << "\n"
				    << "lep2_IdSF"   << setw(6) << lep2_IdSF << "\n"
				    << "lep1_IsoSF"  << setw(6) << lep1_IsoSF << "\n"
				    << "lep2_IsoSF"  << setw(6) << lep2_IsoSF << "\n"
				    << "lep1_FF "    << setw(6) << lep1_FF   << "\n"
				    << "lep2_FF "    << setw(6) << lep2_FF   << "\n"
				    << "totalWeight" << setw(6) << totalWeight << "\n";
	}

	else { // only one muon passes tight id
	  if (leadIsTightMuon) { // leading fakeable muon is tight and the sub-leading muon fails tight id
	    double lep1_IdSF  = (isMC()) ? SFHandler().getIdSF("Medium", lep1.pt, lep1.eta, "Muon") : 1.0;
	    double lep1_IsoSF = (isMC()) ? SFHandler().getIsoSF("Tight", lep1.pt, lep1.eta, "Muon") : 1.0;
	    double lep2_IdSF  = (isMC()) ? SFHandler().getIdSF("Loose", lep2.pt, lep2.eta, "Muon") : 1.0;
	    double lep2_IsoSF = (isMC()) ? SFHandler().getIsoSF("Tight", lep2.pt, lep2.eta, "Muon") : 1.0;
	    double lep2_FF    = SFHandler().getFF(lep2.pt, std::fabs(lep2.eta),"Muon");
	    totalWeight *= lep2_FF * lep1_IdSF * lep1_IsoSF * lep2_IdSF * lep2_IsoSF;
	    if (op.verbose) std::cout << "MuMu\n"
				      << "leading muon is tight | sub leading fails tight id \n"
				      << "lep1_IdSF"   << setw(6) << lep1_IdSF   << "\n"
				      << "lep1_IsoSF"  << setw(6) << lep1_IsoSF  << "\n"
				      << "lep2_IdSF"   << setw(6) << lep2_IdSF   << "\n"
				      << "lep2_IsoSF"  << setw(6) << lep2_IsoSF  << "\n"
				      << "lep2_FF "    << setw(6) << lep2_FF     << "\n"
				      << "totalWeight" << setw(6) << totalWeight << "\n";
	  }
	  else { // leading muon is the fakeable one and the sub-leading is the tight muon
	    double lep1_FF    = SFHandler().getFF(lep1.pt, std::fabs(lep1.eta),"Muon");
	    double lep1_IdSF  = (isMC()) ? SFHandler().getIdSF("Loose", lep1.pt, lep1.eta, "Muon") : 1.0;
	    double lep1_IsoSF = (isMC()) ? SFHandler().getIsoSF("Tight", lep1.pt, lep1.eta, "Muon") : 1.0;
	    double lep2_IdSF  = (isMC()) ? SFHandler().getIdSF("Medium", lep2.pt, lep2.eta, "Muon") : 1.0;
	    double lep2_IsoSF = (isMC()) ? SFHandler().getIsoSF("Tight", lep2.pt, lep2.eta, "Muon") : 1.0;
	    totalWeight *= lep1_FF * lep1_IdSF * lep1_IsoSF * lep2_IdSF * lep2_IsoSF;
	    if (op.verbose) std::cout << "MuMu\n"
				      << "leading muon fails tight id | sub leading is tight \n"
				      << "lep1_FF"     << setw(6) << lep1_FF     << "\n"
				      << "lep1_IdSF"   << setw(6) << lep1_IdSF   << "\n"
				      << "lep1_IsoSF"  << setw(6) << lep1_IsoSF  << "\n"
				      << "lep2_IdSF"   << setw(6) << lep2_IdSF   << "\n"
				      << "lep2_IsoSF"  << setw(6) << lep2_IsoSF  << "\n"
				      << "totalWeight" << setw(6) << totalWeight << "\n";
	  }
	}
      }
      
      else if (isEleMu) { // one of the electron or muon fails tight id
	if (op.verbose) std::cout << " --- EleMu --- \n";
	if (ntel == 0 && ntmu == 0) { // both fail tight ID
	  if (lep1.flavour == 1 && lep2.flavour == 2) { // leading one is muon {flavour = 1} and sub-leading is electron {flavour = 2}
	    double lep1_FF = SFHandler().getFF(lep1.pt, lep1.eta,"Muon");
	    double lep2_FF = SFHandler().getFF(lep2.pt, lep2.eta,"Electron");
	    double lep1_IdSF  = (isMC()) ? SFHandler().getIdSF("Loose", lep1.pt, lep1.eta, "Muon") : 1.0;
	    double lep1_IsoSF = (isMC()) ? SFHandler().getIsoSF("Tight", lep1.pt, lep1.eta, "Muon") : 1.0;
	    double lep2_IdSF  = (isMC()) ? SFHandler().getIdSF("Loose", lep2.pt, lep2.SCeta, "Electron") : 1.0;
	    totalWeight *= -1.0 * lep1_FF * lep2_FF * lep1_IdSF * lep1_IsoSF * lep2_IdSF;
	    if (op.verbose) std::cout << "MuEle\n"
				      << "both fails tight id \n"
				      << "lep1_FF"     << setw(6) << lep1_FF     << "\n"
				      << "lep2_FF"     << setw(6) << lep2_FF     << "\n"
				      << "lep1_IdSF"   << setw(6) << lep1_IdSF   << "\n"
				      << "lep1_IsoSF"  << setw(6) << lep1_IsoSF  << "\n"
				      << "lep2_IdSF"   << setw(6) << lep2_IdSF   << "\n"
				      << "totalWeight" << setw(6) << totalWeight << "\n";
	  }
	  else if (lep1.flavour == 2 && lep2.flavour == 1) {
	    double lep1_FF = SFHandler().getFF(lep1.pt, lep1.eta, "Electron");
	    double lep2_FF = SFHandler().getFF(lep2.pt, lep2.eta, "Muon");
	    double lep1_IdSF  = (isMC()) ? SFHandler().getIdSF("Loose", lep1.pt, lep1.SCeta, "Electron") : 1.0;
	    double lep2_IdSF  = (isMC()) ? SFHandler().getIdSF("Loose", lep2.pt, lep2.eta, "Muon") : 1.0;
	    double lep2_IsoSF = (isMC()) ? SFHandler().getIsoSF("Tight", lep2.pt, lep2.eta, "Muon") : 1.0;
	    totalWeight *= -1.0 * lep1_FF * lep2_FF * lep1_IdSF * lep2_IdSF * lep2_IsoSF;
	    if (op.verbose) std::cout << "EleMu\n"
				      << "both fails tight id \n"
				      << "lep1_FF"     << setw(6) << lep1_FF     << "\n"
				      << "lep2_FF"     << setw(6) << lep2_FF     << "\n"
				      << "lep1_IdSF"   << setw(6) << lep1_IdSF   << "\n"
				      << "lep2_IdSF"   << setw(6) << lep2_IdSF   << "\n"
				      << "lep2_IsoSF"  << setw(6) << lep2_IsoSF  << "\n"
				      << "totalWeight" << setw(6) << totalWeight << "\n";
	  }
	}

	else if (ntel == 1 && ntmu == 0) { // electron is tight but muon is loose
	  if (leadIsTightEle) { // leading fakeable electron is tight and the sub-leading muon fails tight id  
	    double lep1_IdSF = (isMC()) ? SFHandler().getIdSF("Loose", lep1.pt, lep1.eta, "Muon") : 1.0;
	    double lep1_IsoSF = (isMC()) ? SFHandler().getIsoSF("Tight", lep1.pt, lep1.eta, "Muon") : 1.0;
	    double lep2_IdSF = (isMC()) ? SFHandler().getIdSF("Tight", lep2.pt, lep2.SCeta, "Electron") : 1.0;
	    double lep2_FF   = SFHandler().getFF(lep2.pt, std::fabs(lep2.eta),"Muon");
	    totalWeight *= lep1_IdSF * lep1_IsoSF * lep2_IdSF * lep2_FF;
	    if (op.verbose) std::cout << "EleMu\n"
				      << "lead electron is tight | sub-leading muon fails tight id \n"
				      << "lep2_FF"     << setw(6) << lep2_FF     << "\n"
				      << "lep1_IdSF"   << setw(6) << lep1_IdSF   << "\n"
				      << "lep1_IsoSF"  << setw(6) << lep1_IsoSF  << "\n"
				      << "lep2_IdSF"   << setw(6) << lep2_IdSF   << "\n"
				      << "totalWeight" << setw(6) << totalWeight << "\n";
	  }
	  else { // leading muon is the fakeable one and the sub-leading is the tight electron
	    double lep1_IdSF = (isMC()) ? SFHandler().getIdSF("Loose", lep1.pt, lep1.eta, "Muon") : 1.0;
	    double lep1_IsoSF = (isMC()) ? SFHandler().getIsoSF("Tight", lep1.pt, lep1.eta, "Muon") : 1.0;
	    double lep2_IdSF = (isMC()) ? SFHandler().getIdSF("Tight", lep2.pt, lep2.SCeta, "Electron") : 1.0;
	    double lep1_FF   = SFHandler().getFF(lep1.pt, lep1.eta,"Muon");
	    totalWeight *= lep1_IdSF * lep1_IsoSF * lep2_IdSF * lep1_FF;
	    if (op.verbose) std::cout << "MuEle\n"
				      << "lead muon fails tight | sub-leading electron is tight \n"
				      << "lep1_FF"     << setw(6) << lep1_FF     << "\n"
				      << "lep1_IdSF"   << setw(6) << lep1_IdSF   << "\n"
				      << "lep1_IsoSF"  << setw(6) << lep1_IsoSF  << "\n"
				      << "lep2_IdSF"   << setw(6) << lep2_IdSF   << "\n"
				      << "totalWeight" << setw(6) << totalWeight << "\n";
	  }
	}

	else if (ntel == 0 && ntmu == 1) { // electron is losse but muon is tight
	  if (leadIsTightMuon) { // leading fakeable muon is tight and the sub-leading electron fails tight id 
	    double lep1_IdSF  = (isMC()) ? SFHandler().getIdSF("Medium", lep1.pt, lep1.eta, "Muon") : 1.0;
	    double lep1_IsoSF = (isMC()) ? SFHandler().getIsoSF("Tight", lep1.pt, lep1.eta, "Muon") : 1.0;
	    double lep2_IdSF  = (isMC()) ? SFHandler().getIdSF("Loose", lep2.pt, lep2.SCeta, "Electron") : 1.0;
	    double lep2_FF    = SFHandler().getFF(lep2.pt, lep2.eta, "Electron");
	    totalWeight *= lep1_IdSF * lep2_IdSF * lep1_IsoSF * lep2_FF;
	    if (op.verbose) std::cout << "MuEle\n"
				      << "lead muon is tight | sub-leading electron fails tight id \n"
				      << "lep1_IdSF"   << setw(6) << lep1_IdSF   << "\n"
				      << "lep1_IsoSF"  << setw(6) << lep1_IsoSF  << "\n"
				      << "lep2_IdSF"   << setw(6) << lep2_IdSF   << "\n"
				      << "lep2_FF"     << setw(6) << lep2_FF     << "\n"
				      << "totalWeight" << setw(6) << totalWeight << "\n";
	  }
	  else { // leading electron is the fakeable one and the sub-leading is the tight muon
	    double lep1_IdSF  = (isMC()) ? SFHandler().getIdSF("Loose", lep1.pt, lep1.SCeta, "Electron") : 1.0;
	    double lep2_IdSF  = (isMC()) ? SFHandler().getIdSF("Medium", lep2.pt, lep2.eta, "Muon") : 1.0;
	    double lep2_IsoSF = (isMC()) ? SFHandler().getIsoSF("Tight", lep2.pt, lep2.eta, "Muon") : 1.0;
	    double lep1_FF    = SFHandler().getFF(lep1.pt, std::fabs(lep1.SCeta),"Electron");
	    totalWeight *= lep1_IdSF * lep2_IdSF * lep2_IsoSF * lep1_FF;
	    if (op.verbose) std::cout << "EleMu\n"
				      << "lead electron fails tight id | sub-leading muon is tight \n"
				      << "lep1_IdSF"   << setw(6) << lep1_IdSF   << "\n"
				      << "lep2_IdSF"   << setw(6) << lep2_IdSF   << "\n"
				      << "lep2_IsoSF"  << setw(6) << lep2_IsoSF  << "\n"
				      << "lep1_FF"     << setw(6) << lep1_FF     << "\n"
				      << "totalWeight" << setw(6) << totalWeight << "\n";
	  }
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
    
    if (dumpIdx < 10 || op.verbose) {
      cout << "------------------------------------------------------------" << endl;
      cout << "Event : " << evt.event << endl;
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
	cout << setprecision(5) << setw(12) << totalWeight
	     << setw(12) << totalWeight * lumiFac 
	     << endl;
      }
    }
    if (op.verbose || dumpIdx == 0) dumpEvent(evt.event);
    dumpIdx++;
    
    // Region Flags
    // One can add other rehions also e.g. MC_fake region, MC_closure region etc
    std::map <std::string, bool> regionFlags {
       {"_SR_", isSR},
       {"_SB_", (isSB && !isSignal())}
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

    // nEvents in SR
    if (isSR) {
      fillHist1D("evtCutFlow", 12);
      if (isMC()) fillHist1D("evtCutFlowWt", 12, totalWeight);
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
      else               chTag = -1.0;
    }

    TLorentzVector lep1p4 = AnaUtil::getP4(lep1);
    TLorentzVector lep2p4 = AnaUtil::getP4(lep2);
    TLorentzVector metp4;
    metp4.SetPtEtaPhiE(met.pt, 0.0, met.phi, met.pt);
    float dr_l1l2   = lep1p4.DeltaR(lep2p4);
    float dphi_l1l2 = TVector2::Phi_mpi_pi(lep1.phi - lep2.phi);
    float deta_l1l2 = lep1.eta - lep2.eta;
    float invM_l1l2 = (lep1p4+lep2p4).M();
    // ----------------------------- Resolved WZ Region ----------------------------- //
    if (isResolved_WZ) {
      std::map <std::string, bool> ResolvedFlags = AnaUtil::combineMaps(regionFlags, {{"IsResolved_", isResolved_WZ}});

      if (isSR) {
	fillHist1D("evtCutFlow", 13);
	if (isMC()) fillHist1D("evtCutFlowWt", 13, totalWeight);
      }

      fillHist1D("MetPt", met.pt, channelFlags, ResolvedFlags, totalWeight);
      fillHist1D("MetPhi", met.phi, channelFlags, ResolvedFlags, totalWeight);
      fillHist1D("NoAk4Jets", jetColl.size(), channelFlags, ResolvedFlags, totalWeight);
      //int nLoopJets = (jetColl.size() >= 4) ? 4 : jetColl.size();
      float HT = 0.0;
      TLorentzVector jetsp4;
      jetsp4.SetPtEtaPhiM(0.,0.,0.,0.);
      for (size_t i=0; i < jetColl.size(); ++i) {
	auto& jet1 = jetColl[i];
	std::string jetPt_hname  = "Ak4Jet"+std::to_string(i+1)+"Pt";
	fillHist1D (jetPt_hname.c_str(), jet1.pt, channelFlags, ResolvedFlags, totalWeight);
	std::string jetEta_hname  = "Ak4Jet"+std::to_string(i+1)+"Eta";
	fillHist1D (jetEta_hname.c_str(), jet1.eta, channelFlags, ResolvedFlags, totalWeight);
	std::string jetPhi_hname  = "Ak4Jet"+std::to_string(i+1)+"Phi";
	fillHist1D (jetPhi_hname.c_str(), jet1.phi, channelFlags, ResolvedFlags, totalWeight);
	std::string jetBtag_hname  = "Ak4Jet"+std::to_string(i+1)+"btagDeepFlvB";
	fillHist1D (jetBtag_hname.c_str(), jet1.btagDeepFlavB, channelFlags, ResolvedFlags, totalWeight);

	TLorentzVector jet1p4 = AnaUtil::getP4(jet1);
	for (size_t j=i+1; j < jetColl.size(); ++j) {
	  auto& jet2 = jetColl[j];
	  TLorentzVector jet2p4 = AnaUtil::getP4(jet2);
	  std::string invM_hname = "InvM_jet"+std::to_string(i+1)+"_jet"+std::to_string(j+1);
	  std::string dR_hname   = "DR_jet"+std::to_string(i+1)+"_jet"+std::to_string(j+1);
	  std::string dPhi_hname = "DPhi_jet"+std::to_string(i+1)+"_jet"+std::to_string(j+1);
	  std::string dEta_hname = "DEta_jet"+std::to_string(i+1)+"_jet"+std::to_string(j+1);
	  fillHist1D (invM_hname.c_str(), (jet1p4+jet2p4).M(), channelFlags, ResolvedFlags, totalWeight);
	  fillHist1D (dR_hname.c_str(), jet1p4.DeltaR(jet2p4), channelFlags, ResolvedFlags, totalWeight);
	  fillHist1D (dPhi_hname.c_str(), TVector2::Phi_mpi_pi(jet1.phi - jet2.phi), channelFlags, ResolvedFlags, totalWeight);
	  fillHist1D (dEta_hname.c_str(), (jet1.eta - jet2.eta), channelFlags, ResolvedFlags, totalWeight);
	}	
	HT += jet1.pt;
	jetsp4 += AnaUtil::getP4(jet1);

	float lep1jetDR   = lep1p4.DeltaR(jet1p4);
	float lep1jetDPhi = TVector2::Phi_mpi_pi(lep1.phi - jet1.phi);
	float lep1jetDEta = lep1.eta - jet1.eta;
	float lep2jetDR   = lep2p4.DeltaR(jet1p4);
	float lep2jetDPhi = TVector2::Phi_mpi_pi(lep2.phi - jet1.phi);
	float lep2jetDEta = lep2.eta - jet1.eta;
	string lep1jetdr_hname = "DR_lep1jet"+std::to_string(i+1);
	fillHist1D(lep1jetdr_hname.c_str(), lep1jetDR, channelFlags, ResolvedFlags, totalWeight);
	string lep2jetdr_hname = "DR_lep2jet"+std::to_string(i+1);
	fillHist1D(lep2jetdr_hname.c_str(), lep2jetDR, channelFlags, ResolvedFlags, totalWeight);
	string lep1jetdphi_hname = "DPhi_lep1jet"+std::to_string(i+1);
	fillHist1D(lep1jetdphi_hname.c_str(), lep1jetDPhi, channelFlags, ResolvedFlags, totalWeight);
	string lep2jetdphi_hname = "DPhi_lep2jet"+std::to_string(i+1);
	fillHist1D(lep2jetdphi_hname.c_str(), lep2jetDPhi, channelFlags, ResolvedFlags, totalWeight);
	string lep1jetdeta_hname = "DEta_lep1jet"+std::to_string(i+1);
	fillHist1D(lep1jetdeta_hname.c_str(), lep1jetDEta, channelFlags, ResolvedFlags, totalWeight);
	string lep2jetdeta_hname = "DEta_lep2jet"+std::to_string(i+1);
	fillHist1D(lep2jetdeta_hname.c_str(), lep2jetDEta, channelFlags, ResolvedFlags, totalWeight);
      }
  
      // jet inv mass
      float jetsInvM = (jetColl.size() > 3) ? (AnaUtil::getP4(jetColl[0]) + AnaUtil::getP4(jetColl[1]) + 
					       AnaUtil::getP4(jetColl[2]) + AnaUtil::getP4(jetColl[3])).M()
	: (AnaUtil::getP4(jetColl[0]) + AnaUtil::getP4(jetColl[1]) + AnaUtil::getP4(jetColl[2])).M();

      fillHist1D("JetsInvMass", jetsInvM, channelFlags, ResolvedFlags, totalWeight);
      fillHist1D("HT", HT, channelFlags, ResolvedFlags, totalWeight);
      fillHist1D("HT_vectSum", jetsp4.Pt(), channelFlags, ResolvedFlags, totalWeight);

      // Lepton histograms
      fillHist1D("Lep1pt",lep1.pt,channelFlags, ResolvedFlags, totalWeight);
      fillHist1D("Lep2pt",lep2.pt,channelFlags, ResolvedFlags, totalWeight);
      fillHist1D("DiLepPt", (lep1p4+lep2p4).Pt(), channelFlags, ResolvedFlags, totalWeight);
      fillHist1D("Lep1eta",lep1.eta,channelFlags, ResolvedFlags, totalWeight);
      fillHist1D("Lep2eta",lep2.eta,channelFlags, ResolvedFlags, totalWeight);
      fillHist1D("Lep1phi",lep1.phi,channelFlags, ResolvedFlags, totalWeight);
      fillHist1D("Lep2phi",lep2.phi,channelFlags, ResolvedFlags, totalWeight);

      float ST = HT+lep1.pt+lep2.pt+met.pt;
      fillHist1D("ST",ST,channelFlags, ResolvedFlags, totalWeight);
      fillHist1D("DR_lep1lep2", dr_l1l2, channelFlags, ResolvedFlags, totalWeight);
      fillHist1D("DPhi_lep1lep2", dphi_l1l2, channelFlags, ResolvedFlags, totalWeight);
      fillHist1D("DEta_lep1lep2", deta_l1l2, channelFlags, ResolvedFlags, totalWeight);
      fillHist1D("InvM_l1l2", invM_l1l2, channelFlags, ResolvedFlags, totalWeight);

      // --------------------------- Varibales to be plotted and stored in ntuple ------------------------- //      
      
      if (isSR && skimObj_) {
	TreeVariablesResolved varList;

	varList.totalWeight               = totalWeight;
	varList.Channel                = chTag;
	// lepton1
	varList.px_lep1                = lep1p4.Px();
	varList.py_lep1                = lep1p4.Py();
	varList.pz_lep1                = lep1p4.Pz();
	varList.E_lep1                 = lep1p4.E();
	varList.pt_lep1                = lep1.pt;
	varList.eta_lep1               = lep1.eta;
	varList.phi_lep1               = lep1.phi;
	varList.charge_lep1            = lep1.charge;
	varList.pdgid_lep1             = lep1.pdgId;
	// lepton2
	varList.px_lep2                = lep2p4.Px();
	varList.py_lep2                = lep2p4.Py();
	varList.pz_lep2                = lep2p4.Pz();
	varList.E_lep2                 = lep2p4.E();
	varList.pt_lep2                = lep2.pt;
	varList.eta_lep2               = lep2.eta;
	varList.phi_lep2               = lep2.phi;
	varList.charge_lep2            = lep2.charge;
	varList.pdgid_lep2             = lep2.pdgId;
	// jet1
	varList.px_jet1                = AnaUtil::getP4(jetColl[0]).Px();
	varList.py_jet1                = AnaUtil::getP4(jetColl[0]).Py();
	varList.pz_jet1                = AnaUtil::getP4(jetColl[0]).Pz();
	varList.E_jet1                 = AnaUtil::getP4(jetColl[0]).E();
	varList.pt_jet1                = jetColl[0].pt;
	varList.eta_jet1               = jetColl[0].eta;
	varList.phi_jet1               = jetColl[0].phi;
	varList.btag_jet1              = jetColl[0].btagDeepFlavB;
	// jet2
	varList.px_jet2                = AnaUtil::getP4(jetColl[1]).Px();
	varList.py_jet2                = AnaUtil::getP4(jetColl[1]).Py();
	varList.pz_jet2                = AnaUtil::getP4(jetColl[1]).Pz();
	varList.E_jet2                 = AnaUtil::getP4(jetColl[1]).E();
	varList.pt_jet2                = jetColl[1].pt;
	varList.eta_jet2               = jetColl[1].eta;
	varList.phi_jet2               = jetColl[1].phi;
	varList.btag_jet2              = jetColl[1].btagDeepFlavB;
	// jet3
	varList.px_jet3                = AnaUtil::getP4(jetColl[2]).Px();
	varList.py_jet3                = AnaUtil::getP4(jetColl[2]).Py();
	varList.pz_jet3                = AnaUtil::getP4(jetColl[2]).Pz();
	varList.E_jet3                 = AnaUtil::getP4(jetColl[2]).E();
	varList.pt_jet3                = jetColl[2].pt;
	varList.eta_jet3               = jetColl[2].eta;
	varList.phi_jet3               = jetColl[2].phi;
	varList.btag_jet3              = jetColl[2].btagDeepFlavB;
	// jet4
	varList.px_jet4                = (jetColl.size() > 3) ? AnaUtil::getP4(jetColl[3]).Px() : 0.0;
	varList.py_jet4                = (jetColl.size() > 3) ? AnaUtil::getP4(jetColl[3]).Py() : 0.0;
	varList.pz_jet4                = (jetColl.size() > 3) ? AnaUtil::getP4(jetColl[3]).Pz() : 0.0;
	varList.E_jet4                 = (jetColl.size() > 3) ? AnaUtil::getP4(jetColl[3]).E(): 0.0;
	varList.pt_jet4                = (jetColl.size() > 3) ? jetColl[3].pt : 0.0;
	varList.eta_jet4               = (jetColl.size() > 3) ? jetColl[3].eta : 0.0;
	varList.phi_jet4               = (jetColl.size() > 3) ? jetColl[3].phi : 0.0;
	varList.btag_jet4              = (jetColl.size() > 3) ? jetColl[3].btagDeepFlavB : 0.0;
	/*
	// di-jet
	varList.dR_jet1jet2            = (AnaUtil::getP4(jetColl[0])).DeltaR(AnaUtil::getP4(jetColl[1]));
	varList.dR_jet2jet3            = (AnaUtil::getP4(jetColl[0])).DeltaR(AnaUtil::getP4(jetColl[1]));
	varList.dR_jet1jet3            = (AnaUtil::getP4(jetColl[0])).DeltaR(AnaUtil::getP4(jetColl[1]));
	varList.dR_jet2jet4            = (jetColl.size() > 3) ? (AnaUtil::getP4(jetColl[1])).DeltaR(AnaUtil::getP4(jetColl[3])) : 0.0;
	varList.dR_jet3jet4            = (jetColl.size() > 3) ? (AnaUtil::getP4(jetColl[2])).DeltaR(AnaUtil::getP4(jetColl[3])) : 0.0;
	varList.dPhi_jet1jet2          = std::fabs(TVector2::Phi_mpi_pi(jetColl[0].phi - jetColl[1].phi));
	varList.dPhi_jet2jet3          = std::fabs(TVector2::Phi_mpi_pi(jetColl[1].phi - jetColl[2].phi));
	varList.dPhi_jet1jet3          = std::fabs(TVector2::Phi_mpi_pi(jetColl[0].phi - jetColl[2].phi));
	varList.dPhi_jet2jet4          = (jetColl.size() > 3) ? std::fabs(TVector2::Phi_mpi_pi(jetColl[1].phi - jetColl[3].phi)) : 0.0;
	varList.dPhi_jet3jet4          = (jetColl.size() > 3) ? std::fabs(TVector2::Phi_mpi_pi(jetColl[2].phi - jetColl[3].phi)) : 0.0;
	varList.invM_jet1jet2          = (AnaUtil::getP4(jetColl[0]) + AnaUtil::getP4(jetColl[1])).M();
	varList.invM_jet2jet3          = (AnaUtil::getP4(jetColl[1]) + AnaUtil::getP4(jetColl[2])).M();
	varList.invM_jet1jet3          = (AnaUtil::getP4(jetColl[0]) + AnaUtil::getP4(jetColl[2])).M();
	varList.invM_jet2jet4          = (jetColl.size() > 3) ? (AnaUtil::getP4(jetColl[1]) + AnaUtil::getP4(jetColl[3])).M() : 0.0;
	varList.invM_jet3jet4          = (jetColl.size() > 3) ? (AnaUtil::getP4(jetColl[2]) + AnaUtil::getP4(jetColl[3])).M() : 0.0;
	// lep-jet
	*/
	// met
	varList.px_met                 = metp4.Px();
	varList.py_met                 = metp4.Py();
	varList.pz_met                 = 0.0;
	varList.E_met                  = met.pt;
	varList.pt_met                 = met.pt;
	varList.eta_met                = 0.0;
	varList.phi_met                = met.phi;

	skimObj_->fill(varList);    
      }
      histf()->cd(); // Very Very Very Essential
    }

    // ----------------------------- Boosted WZ Region ----------------------------- //
    if (isBoosted_WZ) {
      std::map <std::string, bool> BoostedFlags = AnaUtil::combineMaps(regionFlags, {{"IsBoosted_", isBoosted_WZ}});
      if (isSR) {
	fillHist1D("evtCutFlow", 14);
	if (isMC()) fillHist1D("evtCutFlowWt", 14, totalWeight);
      }

      fillHist1D("MetPt", met.pt, channelFlags, BoostedFlags, totalWeight);
      fillHist1D("MetPhi", met.phi, channelFlags, BoostedFlags, totalWeight);

      fillHist1D("NoAk4Jets", jetColl_ak8Cleaned.size(), channelFlags, BoostedFlags, totalWeight);

      fillHist1D("Lep1pt",lep1.pt,channelFlags, BoostedFlags, totalWeight);
      fillHist1D("Lep2pt",lep2.pt,channelFlags, BoostedFlags, totalWeight);
      fillHist1D("DiLepPt", (lep1p4+lep2p4).Pt(), channelFlags, BoostedFlags, totalWeight);
      fillHist1D("Lep1eta",lep1.eta,channelFlags, BoostedFlags, totalWeight);
      fillHist1D("Lep2eta",lep2.eta,channelFlags, BoostedFlags, totalWeight);
      fillHist1D("Lep1phi",lep1.phi,channelFlags, BoostedFlags, totalWeight);
      fillHist1D("Lep2phi",lep2.phi,channelFlags, BoostedFlags, totalWeight);

      fillHist1D("DR_lep1lep2", dr_l1l2, channelFlags, BoostedFlags, totalWeight);
      fillHist1D("DPhi_lep1lep2", dphi_l1l2, channelFlags, BoostedFlags, totalWeight);
      fillHist1D("DEta_lep1lep2", deta_l1l2, channelFlags, BoostedFlags, totalWeight);
      fillHist1D("InvM_l1l2", invM_l1l2, channelFlags, BoostedFlags, totalWeight);

      fillHist1D("Ak8Jet1Pt", fatJetColl[0].pt, channelFlags, BoostedFlags, totalWeight);
      fillHist1D("Ak8Jet1Eta", fatJetColl[0].eta, channelFlags, BoostedFlags, totalWeight);
      fillHist1D("Ak8Jet1Phi", fatJetColl[0].phi, channelFlags, BoostedFlags, totalWeight);
      if (fatJetColl.size() > 1) {
	fillHist1D("Ak8Jet2Pt", fatJetColl[1].pt, channelFlags, BoostedFlags, totalWeight);
	fillHist1D("Ak8Jet2Eta", fatJetColl[1].eta, channelFlags, BoostedFlags, totalWeight);
	fillHist1D("Ak8Jet2Phi", fatJetColl[1].phi, channelFlags, BoostedFlags, totalWeight);
      }
      if (jetColl_ak8Cleaned.size() >= 1) fillHist1D("Ak4Jet1Pt", jetColl_ak8Cleaned[0].pt, channelFlags, BoostedFlags, totalWeight);
      if (jetColl_ak8Cleaned.size() >= 2) fillHist1D("Ak4Jet2Pt", jetColl_ak8Cleaned[1].pt, channelFlags, BoostedFlags, totalWeight);
      if (fatJetColl.size() == 1) fillHist1D("NoAk4JetsHas1FatJet", jetColl_ak8Cleaned.size(), channelFlags, BoostedFlags, totalWeight);
      if (fatJetColl.size() >= 2) fillHist1D("NoAk4JetsHas2orMoreFatJet", jetColl_ak8Cleaned.size(), channelFlags, BoostedFlags, totalWeight);

      // --------------------------- Varibales to be plotted and stored in ntuple ------------------------- //      
      if (isSR && skimObj_) {
	TreeVariablesBoosted  varList;
	
	varList.totalWeight            = totalWeight;
	varList.Channel                = chTag;
	varList.pt_lep1                = lep1.pt;
	varList.pt_lep2                = lep2.pt;

	skimObj_->fill(varList);
      }
      histf()->cd(); // Very Very Very Essential
      // -------------------------------------------------------------------------------------------------- //
    }
    if (!isMC()) selEvLog() << evt.run << " " << evt.lumis << " " << evt.event << std::endl;
  } //Event loop ends
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
  //std::cout<<lep.genFlv<<std::endl;
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
    "Events processed",
      "has GoodPV",
      "Data Masking",
      "nfakeable leptons >= 2",
      "lep1pt > 25 and lep2pt > 20 GeV",
      "pass HLT",
      "leptons having same charge",
      "low mass resonance veto",
      "Z mass resonance veto",
      "max 2 tight leptons",
      "tau veto",
      "met > 40 GeV",
      "is SR",
      "isResolved_WZ",
      "isBoosted_WZ"
      };
  
  AnaUtil::SetEvtCutFlowBinLabels("evtCutFlow", evLabels);
  AnaUtil::showEfficiency("evtCutFlow", evLabels, "Event Selection (Unweighted)");  
  if (isMC()) {
    cout << endl
         << "evtWeightSum: " << setw(10) << setprecision(0) << evtWeightSum_ << endl
         << "      lumiWt: " << setw(10) << setprecision(5) << lumiFac
         << endl;
    AnaUtil::showEfficiency("evtCutFlowWt", evLabels, "Event Selection (Weighted)", "Events");
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
