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
  
  new TH1D("evtCutFlow", "Event CutFlow", 14, -0.5, 13.5);
  if (isMC()) new TH1D("evtCutFlowWt", "Event CutFlow (Weighted)", 14, -0.5, 13.5);
  new TH1D("evYield_nJetBinned", "", 3, 2.5, 5.5);
  if (isMC()) new TH1D("evYieldWt_nJetBinned", "(Weighted)", 3, 2.5, 5.5);
  new TH1D("evYield_lepFlvBinned", "", 3, -0.5, 2.5);
  if (isMC()) new TH1D("evYieldWt_lepFlvBinned", "(Weighted)", 3, -0.5, 2.5);

  new TH1D("nMuons", "nTightIsoMuons", 10, -0.5, 9.5);
  new TH1D("nElectrons", "nTightIsoElectrons", 10, -0.5, 9.5);
  new TH1D("nLeptons_pc", "nTightIsoLeptons", 10, -0.5, 9.5);
  new TH1D("nJets_", "nTightJetsWithLeptonVeto", 10, -0.5, 9.5);
  new TH1D("nJets_pc", "nTightJetsWithLeptonVeto:isoLepCleaned", 10, -0.5, 9.5);
  new TH1D("nFatJets", "", 10, -0.5, 9.5);

  new TH1D ("fjpt", "", 1000, 0., 2000.);
  new TH1D ("fjeta", "", 100, -5., 5.);
  new TH1D ("fjmass", "", 1000, 0, 2000);
  new TH1D ("fjmassSD", "", 1000, 0, 2000);
  new TH1D ("fjWtag", "", 100, 0, 100);
  new TH1D ("fjZtag", "", 100, 0, 100);
  new TH1D ("fjTtag", "", 100, 0, 100);
  new TH1D ("fjWtagMD", "", 100, 0, 100);
  new TH1D ("fjZtagMD", "", 100, 0, 100);
  new TH1D ("fjTtagMD", "", 100, 0, 100);
  new TH1D ("fjtau1", "", 100, 0, 1);
  new TH1D ("fjtau2", "", 100, 0, 1);
  new TH1D ("fjtau3", "", 100, 0, 1);
  new TH1D ("fjtau3ovtau2", "", 100, 0, 1);
  new TH1D ("fjtau2ovtau1", "", 100, 0, 1);
  new TH1D ("nAk4in", "", 10, -0.5, 9.5);
  new TH1D ("nAk4out", "", 10, -0.5, 9.5);
  new TH2D ("nAk8Ak4in", "", 10, 0.0, 10.0, 10, 0.0, 10.0);
  new TH2D ("nAk8Ak4out", "", 10, 0.0, 10.0, 10, 0.0, 10.0);
  new TH2D ("Ak8ptAk8mSD", "", 400, 200.0, 1000.0, 200, 0.0, 400.0);
  new TH2D ("Ak8mSDt3_t2", "", 200, 0.0, 400.0, 100, 0.0, 1.0);
  new TH2D ("Ak8mSDt2_t1", "", 200, 0.0, 400.0, 100, 0.0, 1.0);
  new TH1D ("fJLepDR", "", 100, -0.5, 4.5);
  new TH1D ("ak4Idx", "", 50, -1.5, 48.5);
  new TH1D ("ak8subJetIdx1", "", 50, -1.5, 48.5);
  new TH1D ("ak8subJetIdx2", "", 50, -1.5, 48.5);

  //nFJ = 1
  new TH1D ("fj1pt_1", "", 1000, 0., 2000.);
  new TH1D ("fj1mSD_1", "", 500, 0., 1000.);
  new TH1D ("fjInvM_1", "", 500, 0., 1000.);
  new TH1D ("nAk4_1", "", 10, -0.5, 9.5);
  new TH1D ("nAk4in_1", "", 10, -0.5, 9.5);
  new TH1D ("ak4inPt_1", "", 1000, 0., 2000.);
  new TH1D ("nAk4out_1", "", 10, -0.5, 9.5);
  new TH1D ("ak4outPt_1", "", 1000, 0., 2000.);
  new TH1D ("ak4invM_1", "", 500, 0., 1000.);
  new TH1D ("ak4inInvM_1", "", 500, 0., 1000.);

  //nFJ = 2
  new TH1D ("fj1pt_2", "", 1000, 0., 2000.);
  new TH1D ("fj2pt_2", "", 1000, 0., 2000.);
  new TH1D ("fj1mSD_2", "", 500, 0., 1000.);
  new TH1D ("fj2mSD_2", "", 500, 0., 1000.);
  new TH1D ("fjInvM_2", "", 500, 0., 1000.);
  new TH1D ("nAk4_2", "", 10, -0.5, 9.5);
  new TH1D ("nAk4in_2", "", 10, -0.5, 9.5);
  new TH1D ("ak4inPt_2", "", 1000, 0., 2000.);
  new TH1D ("nAk4out_2", "", 10, -0.5, 9.5);
  new TH1D ("ak4outPt_2", "", 1000, 0., 2000.);
  new TH1D ("ak4invM_2", "", 500, 0., 1000.);
  new TH1D ("ak4inInvM_2", "", 500, 0., 1000.);

  //nFJ = 3
  new TH1D ("fj1pt_3", "", 1000, 0., 2000.);
  new TH1D ("fj2pt_3", "", 1000, 0., 2000.);
  new TH1D ("fj3pt_3", "", 1000, 0., 2000.);
  new TH1D ("fj1mSD_3", "", 500, 0., 1000.);
  new TH1D ("fj2mSD_3", "", 500, 0., 1000.);
  new TH1D ("fj3mSD_3", "", 500, 0., 1000.);
  new TH1D ("fjInvM_3", "", 500, 0., 1000.);
  new TH1D ("nAk4_3", "", 10, -0.5, 9.5);
  new TH1D ("nAk4in_3", "", 10, -0.5, 9.5);
  new TH1D ("ak4inPt_3", "", 1000, 0., 2000.);
  new TH1D ("nAk4out_3", "", 10, -0.5, 9.5);
  new TH1D ("ak4outPt_3", "", 1000, 0., 2000.);
  new TH1D ("ak4invM_3", "", 500, 0., 1000.);
  new TH1D ("ak4inInvM_3", "", 500, 0., 1000.);

  //nFJ > 3
  new TH1D ("fj1pt_4", "", 1000, 0., 2000.);
  new TH1D ("fj2pt_4", "", 1000, 0., 2000.);
  new TH1D ("fj3pt_4", "", 1000, 0., 2000.);
  new TH1D ("fj4pt_4", "", 1000, 0., 2000.);
  new TH1D ("fj1mSD_4", "", 500, 0., 1000.);
  new TH1D ("fj2mSD_4", "", 500, 0., 1000.);
  new TH1D ("fj3mSD_4", "", 500, 0., 1000.);
  new TH1D ("fj4mSD_4", "", 500, 0., 1000.);
  new TH1D ("fjInvM_4", "", 500, 0., 1000.);
  new TH1D ("nAk4_4", "", 10, -0.5, 9.5);
  new TH1D ("nAk4in_4", "", 10, -0.5, 9.5);
  new TH1D ("ak4inPt_4", "", 1000, 0., 2000.);
  new TH1D ("nAk4out_4", "", 10, -0.5, 9.5);
  new TH1D ("ak4outPt_4", "", 1000, 0., 2000.);
  new TH1D ("ak4invM_4", "", 500, 0., 1000.);
  new TH1D ("ak4inInvM_4", "", 500, 0., 1000.);


  new TH1D("totCharge_pc", "", 5, -2.5, 2.5);
  new TH1D("l1pt", "", 200, 0., 400.);
  new TH1D("l2pt", "", 200, 0., 400.);
  new TH1D("jet1pt_pc", "", 500, 0., 1000.);
  new TH1D("jet2pt_pc", "", 500, 0., 1000.);
  new TH1D("hLepPt", "", 300, 0., 300.);
  new TH1D("h2LepPt", "", 300, 0., 300.);
  new TH1D("MET", "Missing Transverse Energy", 300, 0., 300.);
  new TH1D("lT_pc", "", 300, 0., 300.);
  new TH1D("met_pc", "Missing Transverse Energy", 300, 0., 300.);
  new TH1D("nbJets_pc", "", 5, -0.5, 4.5);

  new TH1D("LT", "", 300, 0., 600.);
  new TH1D("LTMET", "", 300, 0., 600.);
  new TH1D("LTMETvec", "", 300, 0., 600.);
  new TH1D("LTMETfrac", "", 100, 0., 10.);
  new TH1D("LTMETfracInv", "", 100, 0., 1.);
  new TH1D("w1w2mT", "", 500, 0., 1000.);
  new TH1D("w1w2mT_set1", "", 500, 0., 1000.);
  new TH1D("w1w2mT_set2", "", 500, 0., 1000.);
  new TH1D("w1w2mTD", "", 500, 0., 1000.);
  new TH1D("w1w2mToverLT", "", 100, 0., 10.);
  new TH1D("l1l2InvM", "", 200, 0., 200.);
  new TH1D("lepMetAngle", "", 100, 0., 5.);
  new TH1D("jetLepAngle", "", 100, 0., 5.);
  new TH1D("jetMetAngle", "", 100, 0., 5.);
  new TH1D("j1j2MetAngle", "", 100, 0., 5.);

  new TH1D("l1l2DR", "", 100, 0.0, 5.0);
  new TH1D("l1l2DEta", "", 100, -0.5, 4.0);
  new TH1D("l1l2DPhi", "", 100, -0.5, 4.5);

  new TH1D("j1l1DR", "", 100, -0.5, 4.5);
  new TH1D("j1l1DEta", "", 100, -0.5, 4.0);
  new TH1D("j1l1DPhi", "", 100, -0.5, 4.5);

  new TH1D("j1l2DR", "", 100, -0.5, 4.5);
  new TH1D("j1l2DEta", "", 100, -0.5, 4.0);
  new TH1D("j1l2DPhi", "", 100, -0.5, 4.5);

  new TH1D("j2l1DR", "", 100, -0.5, 4.5);
  new TH1D("j2l1DEta", "", 100, -0.5, 4.0);
  new TH1D("j2l1DPhi", "", 100, -0.5, 4.5);

  new TH1D("j2l2DR", "", 100, -0.5, 4.5);
  new TH1D("j2l2DEta", "", 100, -0.5, 4.0);
  new TH1D("j2l2DPhi", "", 100, -0.5, 4.5);

  new TH1D("j3l1DR", "", 100, -0.5, 4.5);
  new TH1D("j3l1DEta", "", 100, -0.5, 4.0);
  new TH1D("j3l1DPhi", "", 100, -0.5, 4.5);

  new TH1D("j3l2DR", "", 100, -0.5, 4.5);
  new TH1D("j3l2DEta", "", 100, -0.5, 4.0);
  new TH1D("j3l2DPhi", "", 100, -0.5, 4.5);

  new TH1D("maxJLDR", "", 100, -0.5, 4.5);
  new TH1D("maxJLDPhi", "", 100, -0.5, 4.5);
  new TH1D("minJLDR", "", 100, -0.5, 4.5);
  new TH1D("minJLDPhi", "", 100, -0.5, 4.5);

  new TH1D("j1METDPhi", "", 100, -0.5, 4.5);
  new TH1D("j2METDPhi", "", 100, -0.5, 4.5);
  new TH1D("j3METDPhi", "", 100, -0.5, 4.5);
  new TH1D("l1METDPhi", "", 100, -0.5, 4.5);
  new TH1D("l2METDPhi", "", 100, -0.5, 4.5);

  new TH1D("minLepMetDPhi", "", 100, -0.5, 4.5);
  new TH1D("maxLepMetDPhi", "", 100, -0.5, 4.5);
  new TH1D("minJetMetDPhi", "", 100, -0.5, 4.5);
  new TH1D("maxJetMetDPhi", "", 100, -0.5, 4.5);

  new TH1D("jet1pt", "", 500, 0., 1000.);
  new TH1D("jet2pt", "", 500, 0., 1000.);
  new TH1D("jet3pt", "", 500, 0., 1000.);
  new TH1D("jet4pt", "", 500, 0., 1000.);
  new TH1D("jet5pt", "", 500, 0., 1000.);
  new TH1D("jet6pt", "", 500, 0., 1000.);
  new TH1D("HT", "", 400, 0., 800.);
  new TH1D("HTvec", "", 400, 0., 800.);
  new TH1D("HTMET", "", 400, 0., 800.);
  new TH1D("HTMETvec", "", 300, 0., 600.);
  new TH1D("HTfrac", "", 100, 0., 10.);
  new TH1D("HTMETfrac", "", 100, 0., 10.);

  new TH1D("jInvM", "", 500, 0., 1000.);
  new TH1D("j1j2InvM", "", 500, 0., 1000.);
  new TH1D("zp1", "", 100, 0., 10.);
  new TH1D("zp2", "", 100, 0., 10.);
  new TH1D("zpMax", "", 100, 0., 10.);
  new TH1D("zpMin", "", 100, 0., 10.);


  new TH1D("j1j2DR", "", 100, -0.5, 4.5);
  new TH1D("j1j2DEta", "", 100, -0.5, 4.0);
  new TH1D("j1j2DPhi", "", 100, -0.5, 4.5);
  new TH1D("j1j3DR", "", 100, -0.5, 4.5);
  new TH1D("j1j3DEta", "", 100, -0.5, 4.0);
  new TH1D("j1j3DPhi", "", 100, -0.5, 4.5);
  new TH1D("j2j3DR", "", 100, -0.5, 4.5);
  new TH1D("j2j3DEta", "", 100, -0.5, 4.0);
  new TH1D("j2j3DPhi", "", 100, -0.5, 4.5);

  new TH1D("minJJDR", "", 100, -0.5, 4.5);
  new TH1D("minJJDPhi", "", 100, -0.5, 4.5);
  new TH1D("maxJJDR", "", 100, -0.5, 4.5);
  new TH1D("maxJJDPhi", "", 100, -0.5, 4.5);


  new TH1D("METovSQHT", "", 500, 0., 50.);
  new TH1D("METovSQLT", "", 500, 0., 50.);
  new TH1D("ST", "", 500, 0., 1000.);
  new TH1D("LTovSqrtST", "", 500, 0., 100.);
  new TH1D("HTovSqrtST", "", 500, 0., 100.);
  new TH1D("STvec", "", 300, 0., 600.);
  new TH1D("STMETvec", "", 300, 0., 600.);

  new TH1D("l1l2InvM_soCsoF_pc", "", 100, 0., 200.);
  new TH1D("l1l2InvM_sCsF", "", 100, 0., 200.);
  new TH1D("l1l2InvM_oCsF", "", 100, 0., 200.);
  new TH1D("ZcandProb", "", 2, -0.5, 1.5);

  new TH2D("LepVsJetDR", "DeltaR between L1L2 vs. J1J2", 100, -0.5, 4.5, 100, -0.5, 4.5);
  new TH2D("jdrVsInvM", "DeltaR between J1J2 vs. J1J2InvM", 100, -0.5, 4.5, 250, 0.0, 500.0);

  new TH1D("bTagFac", "", 100, 0.0, 1.0);
  new TH1D("mvaOut", "", 100, -1., 1.);
  new TH1D("mvaOut_1", "", 100, -1., 1.);

  new TH1D("ST_PBDT", "", 500, 0., 1000.);
  new TH1D("jInvM_PBDT", "", 500, 0., 1000.);
  new TH1D("j1j2InvM_PBDT", "", 500, 0., 1000.);
  new TH1D("w1w2mT_PBDT", "", 500, 0., 1000.);
  new TH1D("w1w2mToverLT_PBDT", "", 100, 0., 10.);
  new TH1D("LTMET_PBDT", "", 300, 0., 600.);
  new TH1D("nJets_PBDT", "", 10, -0.5, 9.5);

  new TH1D("evYield_nBjets", "", 3, -0.5, 2.5);

  histf()->cd();
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
    float puWt = 1.0; //for Data
    float evWt = 1.0; //for Data
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

    // Duplicate event removal for data
    if (!isMC()) {
#ifdef SKIP_DUPLICATE_IN_MEMORY
      std::ostringstream mkey;
      mkey << evt.run << "-" << evt.lumis << "-" << evt.event;
      std::string evs {mkey.str()};
      if (chain()->GetTreeNumber() > skipDuplicateFileindex_ && eventIdStore_.find(evs) != eventIdStore_.end()) {
	// cout << "DuplicateAll: " << evs << endl;
	continue;
      }
      eventIdStore_.insert({evs, 1});
#endif
#ifdef SKIP_DUPLICATE_FROM_FILE
      if (skipDuplicate_) { 
	std::ostringstream mkey;
	mkey << evt.run << "-" << evt.lumis << "-" << evt.event;
	std::string evs {mkey.str()};
	if (eventIdStore_.find(evs) != eventIdStore_.end()) {
	  // cout << "Duplicate Event: " << evs << endl;
	  continue;
	}
      }
      evLog() << evt.run << " " << evt.lumis << " " << evt.event << std::endl;
#endif
    }

    // HLT Selection
    std::vector <bool> hltScores = HLT_Scores;
    bool HLT_pass {false};
    for (auto&& hlt : hltScores){
      if (hlt) HLT_pass = true;
    }
    if (!HLT_pass) continue;

    //Making the object collections ready!!!
    findObjects();

    //if (ev < 100) dumpEverything (ev, fLog());

    histf()->cd(); //required
    histf()->cd("TMVAnalysis");

    //Access Selected Objects
    const auto& looseIsoEleList = getTightIsoEleList();
    int nEle = looseIsoEleList.size();
    const auto& mediumIsoMuList = getTightIsoMuList();
    int nMu = mediumIsoMuList.size();
    //const auto& tightJetListLV = getTightLepVetoJetList();
    const auto& tightJetListLV = getCleanJetList();
    int nJet = tightJetListLV.size();
    //const auto& fatJetList = getFatJetList();
    //int nfJet = fatJetList.size();
    //    const auto& isoTauList = getIsoTauList();
    //    int nTau = isoTauList.size();
    const vhtm::MET& met = getMETList().at(0);

    //P A C K I N G  L E P T O N S in LepCandList_
    std::vector<LeptonCand>LepCandList;
    if (nMu > 0) packLeptons<vhtm::Muon>(mediumIsoMuList, LepCandList);
    if (nEle > 0)  packLeptons<vhtm::Electron>(looseIsoEleList, LepCandList);
    AnaUtil::fillHist1D("nLeptonCand", LepCandList.size(), allWt);

    std::sort(std::begin(LepCandList), std::end(LepCandList), PtComparator<LeptonCand>()); //sorting lepton candidates

    AnaUtil::fillHist1D ("nMuons", nMu, allWt);
    AnaUtil::fillHist1D ("nElectrons", nEle, allWt);
    AnaUtil::fillHist1D ("nLeptons_pc", LepCandList.size(), allWt);
    //AnaUtil::fillHist1D ("nFatJets", nfJet, allWt);

    if (LepCandList.size() < 2) continue;
    AnaUtil::fillHist1D("evtCutFlow", 2);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 2, allWt);

    AnaUtil::fillHist1D ("l1pt", LepCandList[0].pt, allWt);
    if (LepCandList[0].pt < 25) continue;
    AnaUtil::fillHist1D("evtCutFlow", 3);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 3, allWt);

    AnaUtil::fillHist1D ("l2pt", LepCandList[1].pt, allWt);
    if (LepCandList[1].pt < 15) continue;
    AnaUtil::fillHist1D("evtCutFlow", 4);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 4, allWt);
    
    auto& lep1 = LepCandList[0];
    auto& lep2 = LepCandList[1];
    TLorentzVector lep1p4 = AnaUtil::getP4(lep1);
    TLorentzVector lep2p4 = AnaUtil::getP4(lep2);

    //AnaUtil::fillHist1D ("l1l2InvM", l1l2InvM, allWt);

    float totCharge = lep1.charge + lep2.charge;
    AnaUtil::fillHist1D ("totCharge_pc", totCharge, allWt);
    if (totCharge == 0.0) continue;
    AnaUtil::fillHist1D("evtCutFlow", 5);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 5, allWt);

    if (hasZcandidate(LepCandList, allWt)) continue;
    AnaUtil::fillHist1D("evtCutFlow", 6);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 6, allWt);

    if (met.pt < 50) continue;
    float hT = 0.0;
    TLorentzVector jSumP4;
    for (int j = 0; j < nJet; ++j){
      auto& jet = tightJetListLV[j];
      if (j < 4) {
	jSumP4 += AnaUtil::getP4(jet);
	hT += jet.pt;
      }
      if (j == 0) AnaUtil::fillHist1D ("jet1pt", jet.pt, allWt);
      else if (j == 1) AnaUtil::fillHist1D ("jet2pt", jet.pt, allWt);
      else if (j == 2) AnaUtil::fillHist1D ("jet3pt", jet.pt, allWt);
      else if (j == 3) AnaUtil::fillHist1D ("jet4pt", jet.pt, allWt);
      else if (j == 4) AnaUtil::fillHist1D ("jet5pt", jet.pt, allWt);
      else if (j == 5) AnaUtil::fillHist1D ("jet6pt", jet.pt, allWt);
    }
    /*

    if (nfJet > 0) 
      for (auto& j: fatJetList) {
	AnaUtil::fillHist1D ("fjpt", j.pt, allWt);
	AnaUtil::fillHist1D ("fjeta", j.eta, allWt);
	//	if (j.pt < 200) continue;
	//	if (std::fabs(j.eta) > 2.4) continue;
	AnaUtil::fillHist1D ("fjmass", j.mass, allWt);
	AnaUtil::fillHist1D ("fjmassSD", j.softDropMass, allWt);
	AnaUtil::fillHist1D ("fjWtag", j.deepTag_WvsQCD, allWt);
	AnaUtil::fillHist1D ("fjZtag", j.deepTag_ZvsQCD, allWt);
	AnaUtil::fillHist1D ("fjTtag", j.deepTag_TvsQCD, allWt);
	AnaUtil::fillHist1D ("fjWtagMD", j.deepTagMD_WvsQCD, allWt);
	AnaUtil::fillHist1D ("fjZtagMD", j.deepTagMD_ZvsQCD, allWt);
	AnaUtil::fillHist1D ("fjTtagMD", j.deepTagMD_TvsQCD, allWt);
	AnaUtil::fillHist1D ("fjtau1", j.tau1, allWt);
	AnaUtil::fillHist1D ("fjtau2", j.tau2, allWt);
	AnaUtil::fillHist1D ("fjtau3", j.tau3, allWt);
	AnaUtil::fillHist1D ("fjtau3ovtau2", j.tau3/j.tau2, allWt);
	AnaUtil::fillHist1D ("fjtau2ovtau1", j.tau2/j.tau1, allWt);
	AnaUtil::fillHist2D ("Ak8ptAk8mSD", j.pt, j.softDropMass, allWt);
	AnaUtil::fillHist2D ("Ak8mSDt3_t2", j.softDropMass, (j.tau3/j.tau2), allWt);
	AnaUtil::fillHist2D ("Ak8mSDt2_t1", j.softDropMass, (j.tau2/j.tau1), allWt);
	for (auto& l: LepCandList) {
	  AnaUtil::fillHist1D ("fJLepDR", AnaUtil::getP4(l).DeltaR(AnaUtil::getP4(j)), allWt);
	}
	AnaUtil::fillHist1D ("ak8subJetIdx1", j.subJetIdx1, allWt);
	AnaUtil::fillHist1D ("ak8subJetIdx2", j.subJetIdx2, allWt);
      }

    if (fatJetList.size() > 0) {
      int nAk4out = 0;
      int nAk4in = 0;
      TLorentzVector Ak4inP4;
      for (auto& ak4: tightJetListLV) {
	bool isSubJet {false};
	AnaUtil::fillHist1D ("ak4Idx", ak4.index, allWt);
	for (auto& ak8: fatJetList) {
	  if ((ak4.index == ak8.subJetIdx1) || (ak4.index == ak8.subJetIdx2)) {
	    isSubJet = true;
	    nAk4in++;
	    Ak4inP4 += AnaUtil::getP4 (ak4);
	    if (fatJetList.size() == 1) AnaUtil::fillHist1D ("ak4inPt_1", ak4.pt, allWt);
	    else if (fatJetList.size() == 2) AnaUtil::fillHist1D ("ak4inPt_2", ak4.pt, allWt);
	    else if (fatJetList.size() == 3) AnaUtil::fillHist1D ("ak4inPt_3", ak4.pt, allWt);
	    else if (fatJetList.size() > 3) AnaUtil::fillHist1D ("ak4inPt_4", ak4.pt, allWt);
	    break;
	  }
	}
	if (isSubJet) continue;
	nAk4out++;
	if (fatJetList.size() == 1) AnaUtil::fillHist1D ("ak4outPt_1", ak4.pt, allWt);
	else if (fatJetList.size() == 2) AnaUtil::fillHist1D ("ak4outPt_2", ak4.pt, allWt);
	else if (fatJetList.size() == 3) AnaUtil::fillHist1D ("ak4outPt_3", ak4.pt, allWt);
	else if (fatJetList.size() > 3) AnaUtil::fillHist1D ("ak4outPt_4", ak4.pt, allWt);
      }
      
      AnaUtil::fillHist1D ("nAk4in", nAk4in, allWt);
      AnaUtil::fillHist1D ("nAk4out", nAk4out, allWt);
      AnaUtil::fillHist2D ("nAk8Ak4in", fatJetList.size(), nAk4in, allWt);
      AnaUtil::fillHist2D ("nAk8Ak4out", fatJetList.size(), nAk4out, allWt);     
      
      if (fatJetList.size() == 1) {
	AnaUtil::fillHist1D ("fj1pt_1", fatJetList[0].pt, allWt);
	AnaUtil::fillHist1D ("fj1mSD_1", fatJetList[0].softDropMass, allWt);
	AnaUtil::fillHist1D ("fjInvM_1", AnaUtil::getP4(fatJetList[0]).M(), allWt);
	AnaUtil::fillHist1D ("nAk4_1", tightJetListLV.size(), allWt);
	AnaUtil::fillHist1D ("nAk4in_1", nAk4in, allWt);
	AnaUtil::fillHist1D ("nAk4out_1", nAk4out, allWt);
	AnaUtil::fillHist1D ("ak4invM_1", jSumP4.M(), allWt);
	AnaUtil::fillHist1D ("ak4inInvM_1", Ak4inP4.M(), allWt);
      }
      else if (fatJetList.size() == 2) {
	AnaUtil::fillHist1D ("fj1pt_2", fatJetList[0].pt, allWt);
	AnaUtil::fillHist1D ("fj2pt_2", fatJetList[1].pt, allWt);
	AnaUtil::fillHist1D ("fj1mSD_2", fatJetList[0].softDropMass, allWt);
	AnaUtil::fillHist1D ("fj2mSD_2", fatJetList[1].softDropMass, allWt);
	AnaUtil::fillHist1D ("fjInvM_2", (AnaUtil::getP4(fatJetList[0])+AnaUtil::getP4(fatJetList[1])).M(), allWt);
	AnaUtil::fillHist1D ("nAk4_2", tightJetListLV.size(), allWt);
	AnaUtil::fillHist1D ("nAk4in_2", nAk4in, allWt);
	AnaUtil::fillHist1D ("nAk4out_2", nAk4out, allWt);
	AnaUtil::fillHist1D ("ak4invM_2", jSumP4.M(), allWt);
	AnaUtil::fillHist1D ("ak4inInvM_2", Ak4inP4.M(), allWt);
      }
      else if (fatJetList.size() == 3) {
	AnaUtil::fillHist1D ("fj1pt_3", fatJetList[0].pt, allWt);
	AnaUtil::fillHist1D ("fj2pt_3", fatJetList[1].pt, allWt);
	AnaUtil::fillHist1D ("fj3pt_3", fatJetList[2].pt, allWt);
	AnaUtil::fillHist1D ("fj1mSD_3", fatJetList[0].softDropMass, allWt);
	AnaUtil::fillHist1D ("fj2mSD_3", fatJetList[1].softDropMass, allWt);
	AnaUtil::fillHist1D ("fj3mSD_3", fatJetList[2].softDropMass, allWt);
	AnaUtil::fillHist1D ("fjInvM_3", (AnaUtil::getP4(fatJetList[0])+AnaUtil::getP4(fatJetList[1])+AnaUtil::getP4(fatJetList[2])).M(), allWt);
	AnaUtil::fillHist1D ("nAk4_3", tightJetListLV.size(), allWt);
	AnaUtil::fillHist1D ("nAk4in_3", nAk4in, allWt);
	AnaUtil::fillHist1D ("nAk4out_3", nAk4out, allWt);
	AnaUtil::fillHist1D ("ak4invM_3", jSumP4.M(), allWt);
	AnaUtil::fillHist1D ("ak4inInvM_3", Ak4inP4.M(), allWt);
      }
      else if (fatJetList.size() > 3) {
	AnaUtil::fillHist1D ("fj1pt_4", fatJetList[0].pt, allWt);
	AnaUtil::fillHist1D ("fj2pt_4", fatJetList[1].pt, allWt);
	AnaUtil::fillHist1D ("fj3pt_4", fatJetList[2].pt, allWt);
	AnaUtil::fillHist1D ("fj4pt_4", fatJetList[3].pt, allWt);
	AnaUtil::fillHist1D ("fj1mSD_4", fatJetList[0].softDropMass, allWt);
	AnaUtil::fillHist1D ("fj2mSD_4", fatJetList[1].softDropMass, allWt);
	AnaUtil::fillHist1D ("fj3mSD_4", fatJetList[2].softDropMass, allWt);
	AnaUtil::fillHist1D ("fj4mSD_4", fatJetList[3].softDropMass, allWt);
	AnaUtil::fillHist1D ("fjInvM_4", (AnaUtil::getP4(fatJetList[0])+AnaUtil::getP4(fatJetList[1])+AnaUtil::getP4(fatJetList[2])+AnaUtil::getP4(fatJetList[3])).M(), allWt);
	AnaUtil::fillHist1D ("nAk4_4", tightJetListLV.size(), allWt);
	AnaUtil::fillHist1D ("nAk4in_4", nAk4in, allWt);
	AnaUtil::fillHist1D ("nAk4out_4", nAk4out, allWt);
	AnaUtil::fillHist1D ("ak4invM_4", jSumP4.M(), allWt);
	AnaUtil::fillHist1D ("ak4inInvM_4", Ak4inP4.M(), allWt);
      }

      ///////////////////////
      if (Ak4inIdx.size() > 0) {
	TLorentzVector Ak4inP4;
	for (size_t i = 0; i < Ak4inIdx.size(); ++i) {
	  Ak4inP4 += AnaUtil::getP4 (tightJetListLV.at(Ak4inIdx.at(i)));
	}
	
	if (fatJetList.size() == 1) AnaUtil::fillHist1D ("ak4inInvM_1", Ak4inP4.M(), allWt);
	else if (fatJetList.size() == 2) AnaUtil::fillHist1D ("ak4inInvM_2", Ak4inP4.M(), allWt);
	else if (fatJetList.size() == 3) AnaUtil::fillHist1D ("ak4inInvM_3", Ak4inP4.M(), allWt);
	else if (fatJetList.size() > 3) AnaUtil::fillHist1D ("ak4inInvM_4", Ak4inP4.M(), allWt);
	
	}/////////////////
    }
    
    */
    //AnaUtil::fillHist1D ("nJets_", nJet_, allWt);
    AnaUtil::fillHist1D ("nJets_pc", nJet, allWt);
    if (nJet < 3) continue;
    AnaUtil::fillHist1D("evtCutFlow", 7);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 7, allWt);

    TLorentzVector metp4;
    metp4.SetPtEtaPhiE(met.pt, 0.0, met.phi, met.pt);



    AnaUtil::fillHist1D ("jet1pt_pc", tightJetListLV[0].pt, allWt);
    if (tightJetListLV[0].pt < 60) continue;
    AnaUtil::fillHist1D("evtCutFlow", 8);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 8, allWt);

    AnaUtil::fillHist1D ("jet2pt_pc", tightJetListLV[1].pt, allWt);
    if (tightJetListLV[1].pt < 35) continue;
    AnaUtil::fillHist1D("evtCutFlow", 9);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 9, allWt);

    float lT = lep1.pt + lep2.pt;
    AnaUtil::fillHist1D("lT_pc", lT, allWt);
    //if (lT < 65) continue;
    AnaUtil::fillHist1D("evtCutFlow", 10);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 10, allWt);

    AnaUtil::fillHist1D ("met_pc", met.pt, allWt);
    //    if (met.pt < 50) continue;
    AnaUtil::fillHist1D("evtCutFlow", 11);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 11, allWt);

    int nbJets = 0;
    for (auto& j: tightJetListLV){
      float bTag = j.btagCSVV2;
      AnaUtil::fillHist1D ("bTagFac", bTag, allWt);
      if (bTag < AnaUtil::cutValue(jetCutMap(), "btagFactor")) continue;
      nbJets++;
    }
    AnaUtil::fillHist1D ("nbJets_pc", nbJets, allWt);

    if (nbJets == 0) AnaUtil::fillHist1D ("evYield_nBjets", 0, allWt);
    if (nbJets == 1) AnaUtil::fillHist1D ("evYield_nBjets", 1, allWt);
    if (nbJets >= 2) AnaUtil::fillHist1D ("evYield_nBjets", 2, allWt);

    float lepDR = lep1p4.DeltaR(lep2p4);
    float lepDEta = std::fabs(lep1.eta - lep2.eta);
    float lepDPhi = std::fabs(TVector2::Phi_mpi_pi(lep1.phi - lep2.phi));
    float mT  = Calculate_TotalMT(lep1p4, lep2p4, metp4);
    float mTD = Calculate_TotalMT_Doubt(lep1p4, lep2p4, metp4);
    float mT_set1 = Calculate_MT((lep1p4+lep2p4), metp4);
    float mT_set2 = (lep1p4+lep2p4+metp4).Mt();
    float mTovLt = mT/lT;
    float LTMET = lT + met.pt;
    float LTMETvec = (lep1p4 + lep2p4 + metp4).Pt();
    float LTMETfrac = LTMET/LTMETvec;
    float LTMETfracInv = 1.0/LTMETfrac;
    float l1l2InvM = (lep1p4 + lep2p4).M();

    TVector3 lepUnitVec = (lep1p4+lep2p4).Vect().Unit();
    TVector3 lepUnitVecOrtho = lepUnitVec.Orthogonal();
    TVector3 metUnitVec = metp4.Vect().Unit();
    TVector3 j1j2UnitVecOrtho = (AnaUtil::getP4(tightJetListLV[0]) + AnaUtil::getP4(tightJetListLV[1])).Vect().Unit().Orthogonal();
    TVector3 jetUnitVec = jSumP4.Vect().Unit();
    TVector3 jetUnitVecOrtho = jetUnitVec.Orthogonal();
    double ang = lepUnitVecOrtho.Angle(metUnitVec);
    double j1j2MetAng = j1j2UnitVecOrtho.Angle(metUnitVec);
    double jetMetAng = jetUnitVecOrtho.Angle(metUnitVec);
    double jetLepAng = jetUnitVecOrtho.Angle(lepUnitVecOrtho);
    


    std::vector<float>jlDR, jlDPhi;
    float j1l1DR = (AnaUtil::getP4(tightJetListLV[0])).DeltaR(lep1p4);
    float j1l1DPhi = std::fabs(TVector2::Phi_mpi_pi(tightJetListLV[0].phi-lep1.phi));
    float j1l1DEta = std::fabs(tightJetListLV[0].eta - lep1.eta);
    jlDR.push_back(j1l1DR);
    jlDPhi.push_back(j1l1DPhi);

    float j1l2DR = (AnaUtil::getP4(tightJetListLV[0])).DeltaR(lep2p4);
    float j1l2DPhi = std::fabs(TVector2::Phi_mpi_pi(tightJetListLV[0].phi-lep2.phi));
    float j1l2DEta = std::fabs(tightJetListLV[0].eta - lep2.eta);
    jlDR.push_back(j1l2DR);
    jlDPhi.push_back(j1l2DPhi);

    float j2l1DR = (AnaUtil::getP4(tightJetListLV[1])).DeltaR(lep1p4);
    float j2l1DPhi = std::fabs(TVector2::Phi_mpi_pi(tightJetListLV[1].phi-lep1.phi));
    float j2l1DEta = std::fabs(tightJetListLV[1].eta - lep1.eta);
    jlDR.push_back(j2l1DR);
    jlDPhi.push_back(j2l1DPhi);

    float j2l2DR = (AnaUtil::getP4(tightJetListLV[1])).DeltaR(lep2p4);
    float j2l2DPhi = std::fabs(TVector2::Phi_mpi_pi(tightJetListLV[1].phi-lep2.phi));
    float j2l2DEta = std::fabs(tightJetListLV[1].eta - lep2.eta);
    jlDR.push_back(j2l2DR);
    jlDPhi.push_back(j2l2DPhi);

    float j3l1DR = (AnaUtil::getP4(tightJetListLV[2])).DeltaR(lep1p4);
    float j3l1DPhi = std::fabs(TVector2::Phi_mpi_pi(tightJetListLV[2].phi-lep1.phi));
    float j3l1DEta = std::fabs(tightJetListLV[2].eta - lep1.eta);
    jlDR.push_back(j3l1DR);
    jlDPhi.push_back(j3l1DPhi);

    float j3l2DR = (AnaUtil::getP4(tightJetListLV[2])).DeltaR(lep2p4);
    float j3l2DPhi = std::fabs(TVector2::Phi_mpi_pi(tightJetListLV[2].phi-lep2.phi));
    float j3l2DEta = std::fabs(tightJetListLV[2].eta - lep2.eta);
    jlDR.push_back(j3l2DR);
    jlDPhi.push_back(j3l2DPhi);

    std::sort(jlDR.begin(), jlDR.end());
    std::sort(jlDPhi.begin(), jlDPhi.end());
    float maxJLDR = jlDR[5];
    float minJLDR = jlDR[0];
    float maxJLDPhi = jlDPhi[5];
    float minJLDPhi = jlDPhi[0];

    float j1j2InvM = (AnaUtil::getP4(tightJetListLV[0]) + AnaUtil::getP4(tightJetListLV[1])).M();
    //AN-16/428
    //zeppenfeld variable
    float zp1 = std::fabs(lep1p4.Eta() - 0.5*(tightJetListLV[0].eta + tightJetListLV[1].eta))/std::fabs(tightJetListLV[0].eta - tightJetListLV[1].eta);
    float zp2 = std::fabs(lep2p4.Eta() - 0.5*(tightJetListLV[0].eta + tightJetListLV[1].eta))/std::fabs(tightJetListLV[0].eta - tightJetListLV[1].eta);
    float zpMax = (zp1 >= zp2) ? zp1 : zp2;
    float zpMin = (zp1 >= zp2) ? zp2 : zp1;

    float lep1MetDPhi = std::fabs(TVector2::Phi_mpi_pi(lep1.phi - met.phi));
    float lep2MetDPhi = std::fabs(TVector2::Phi_mpi_pi(lep2.phi - met.phi));
    std::vector <float> JMDphi;
    float jet1MetDPhi = std::fabs(TVector2::Phi_mpi_pi(tightJetListLV[0].phi - met.phi));
    float jet2MetDPhi = std::fabs(TVector2::Phi_mpi_pi(tightJetListLV[1].phi - met.phi));
    float jet3MetDPhi = std::fabs(TVector2::Phi_mpi_pi(tightJetListLV[2].phi - met.phi));
    JMDphi.push_back(jet1MetDPhi);
    JMDphi.push_back(jet2MetDPhi);
    JMDphi.push_back(jet3MetDPhi);
    std::sort(JMDphi.begin(), JMDphi.end());
    
    float minLepMetDPhi = (lep1MetDPhi >= lep2MetDPhi) ? lep2MetDPhi : lep1MetDPhi;
    float maxLepMetDPhi = (lep1MetDPhi >= lep2MetDPhi) ? lep1MetDPhi : lep2MetDPhi;
    //float minJetMetDPhi = (jet1MetDPhi >= jet2MetDPhi) ? jet2MetDPhi : jet1MetDPhi;
    //float maxJetMetDPhi = (jet1MetDPhi >= jet2MetDPhi) ? jet1MetDPhi : jet2MetDPhi;
    float minJetMetDPhi = JMDphi[0];
    float maxJetMetDPhi = JMDphi[2];





    float hTMET = hT+met.pt;
    float hTvec = jSumP4.Pt();
    float hTMETvec = (jSumP4+metp4).Pt();
    float hTfrac = hT/hTvec;
    float hTMETfrac = hTMETvec/hTvec;
    float jInvM = jSumP4.M();
    std::vector<float>jjDR, jjDPhi;
    float j1j2DR = (AnaUtil::getP4(tightJetListLV[0])).DeltaR(AnaUtil::getP4(tightJetListLV[1]));
    jjDR.push_back(j1j2DR);
    float j1j2DEta = std::fabs(tightJetListLV[0].eta - tightJetListLV[1].eta);
    float j1j2DPhi = std::fabs(TVector2::Phi_mpi_pi(tightJetListLV[0].phi-tightJetListLV[1].phi));
    jjDPhi.push_back(j1j2DPhi);
    float j1j3DR = (AnaUtil::getP4(tightJetListLV[0])).DeltaR(AnaUtil::getP4(tightJetListLV[2]));
    jjDR.push_back(j1j3DR);
    float j1j3DEta = std::fabs(tightJetListLV[0].eta - tightJetListLV[2].eta);
    float j1j3DPhi = std::fabs(TVector2::Phi_mpi_pi(tightJetListLV[0].phi-tightJetListLV[2].phi));
    jjDPhi.push_back(j1j3DPhi);
    float j2j3DR = (AnaUtil::getP4(tightJetListLV[1])).DeltaR(AnaUtil::getP4(tightJetListLV[2]));
    jjDR.push_back(j2j3DR);
    float j2j3DEta = std::fabs(tightJetListLV[1].eta - tightJetListLV[2].eta);
    float j2j3DPhi = std::fabs(TVector2::Phi_mpi_pi(tightJetListLV[1].phi-tightJetListLV[2].phi));
    jjDPhi.push_back(j2j3DPhi);
    std::sort(jjDR.begin(), jjDR.end());
    std::sort(jjDPhi.begin(), jjDPhi.end());

    float minJJDR   = jjDR[0];
    float maxJJDR   = jjDR[2];
    float minJJDPhi = jjDPhi[0];
    float maxJJDPhi = jjDPhi[2];

    float METovSQHT  = met.pt/TMath::Sqrt(hT);
    float METovSQLT  = met.pt/TMath::Sqrt(lT);
    float ST         = hT+lT;
    float LTovSqrtST = lT/TMath::Sqrt(ST);
    float HTovSqrtST = hT/TMath::Sqrt(ST);
    float STvec      = (lep2p4+lep1p4+jSumP4).Pt();
    float STMETvec   = (lep2p4+lep1p4+jSumP4+metp4).Pt();




    AnaUtil::fillHist1D ("nLeptons", LepCandList.size(), allWt);
    AnaUtil::fillHist1D ("nJets", nJet, allWt);
    AnaUtil::fillHist1D ("MET", met.pt, allWt);
    AnaUtil::fillHist1D ("hLepPt", lep1.pt, allWt);
    AnaUtil::fillHist1D ("h2LepPt", lep2.pt, allWt);
    AnaUtil::fillHist1D ("LT", lT, allWt);
    AnaUtil::fillHist1D ("w1w2mT", mT, allWt);
    AnaUtil::fillHist1D ("w1w2mT_set1", mT_set1, allWt);
    AnaUtil::fillHist1D ("w1w2mT_set2", mT_set2, allWt);
    AnaUtil::fillHist1D ("lepMetAngle", ang, allWt);
    AnaUtil::fillHist1D ("j1j2MetAngle", j1j2MetAng, allWt);
    AnaUtil::fillHist1D ("jetMetAngle", jetMetAng, allWt);
    AnaUtil::fillHist1D ("jetLepAngle", jetLepAng, allWt);

    AnaUtil::fillHist1D ("w1w2mTD", mTD, allWt);
    AnaUtil::fillHist1D ("w1w2mToverLT", mTovLt, allWt);
    AnaUtil::fillHist1D ("LTMET", LTMET, allWt);
    AnaUtil::fillHist1D ("LTMETvec", LTMETvec, allWt);
    AnaUtil::fillHist1D ("LTMETfrac", LTMETfrac, allWt);
    AnaUtil::fillHist1D ("LTMETfracInv", LTMETfracInv, allWt);
    AnaUtil::fillHist1D ("l1l2InvM", l1l2InvM, allWt);
    AnaUtil::fillHist1D ("jInvM", jInvM, allWt);

    AnaUtil::fillHist1D ("l1l2DR", lepDR, allWt);
    AnaUtil::fillHist1D ("l1l2DEta", lepDEta, allWt);
    AnaUtil::fillHist1D ("l1l2DPhi", lepDPhi, allWt);

    AnaUtil::fillHist2D ("LepVsJetDR", lepDR, j1j2DR, allWt);
    AnaUtil::fillHist2D ("jdrVsInvM", j1j2DR, j1j2InvM, allWt);


    AnaUtil::fillHist1D ("j1l1DR", j1l1DR, allWt);
    AnaUtil::fillHist1D ("j1l1DEta", j1l1DEta, allWt);
    AnaUtil::fillHist1D ("j1l1DPhi", j1l1DPhi, allWt);

    AnaUtil::fillHist1D ("j1l2DR", j1l2DR, allWt);
    AnaUtil::fillHist1D ("j1l2DEta", j1l2DEta, allWt);
    AnaUtil::fillHist1D ("j1l2DPhi", j1l2DPhi, allWt);

    AnaUtil::fillHist1D ("j2l1DR", j2l1DR, allWt);
    AnaUtil::fillHist1D ("j2l1DEta", j2l1DEta, allWt);
    AnaUtil::fillHist1D ("j2l1DPhi", j2l1DPhi, allWt);

    AnaUtil::fillHist1D ("j2l2DR", j2l2DR, allWt);
    AnaUtil::fillHist1D ("j2l2DEta", j2l2DEta, allWt);
    AnaUtil::fillHist1D ("j2l2DPhi", j2l2DPhi, allWt);

    AnaUtil::fillHist1D ("j3l1DR", j3l1DR, allWt);
    AnaUtil::fillHist1D ("j3l1DEta", j3l1DEta, allWt);
    AnaUtil::fillHist1D ("j3l1DPhi", j3l1DPhi, allWt);

    AnaUtil::fillHist1D ("j3l2DR", j3l2DR, allWt);
    AnaUtil::fillHist1D ("j3l2DEta", j3l2DEta, allWt);
    AnaUtil::fillHist1D ("j3l2DPhi", j3l2DPhi, allWt);

    AnaUtil::fillHist1D ("maxJLDR", maxJLDR, allWt);
    AnaUtil::fillHist1D ("minJLDR", minJLDR, allWt);
    AnaUtil::fillHist1D ("maxJLDPhi", maxJLDPhi, allWt);
    AnaUtil::fillHist1D ("minJLDPhi", minJLDPhi, allWt);

    AnaUtil::fillHist1D ("l1METDPhi", lep1MetDPhi, allWt);
    AnaUtil::fillHist1D ("l2METDPhi", lep2MetDPhi, allWt);
    AnaUtil::fillHist1D ("j1METDPhi", jet1MetDPhi, allWt);
    AnaUtil::fillHist1D ("j2METDPhi", jet2MetDPhi, allWt);
    AnaUtil::fillHist1D ("j3METDPhi", jet3MetDPhi, allWt);

    AnaUtil::fillHist1D ("minLepMetDPhi", minLepMetDPhi, allWt);
    AnaUtil::fillHist1D ("maxLepMetDPhi", maxLepMetDPhi, allWt);
    AnaUtil::fillHist1D ("minJetMetDPhi", minJetMetDPhi, allWt);
    AnaUtil::fillHist1D ("maxJetMetDPhi", maxJetMetDPhi, allWt);

    AnaUtil::fillHist1D ("j1j2InvM", j1j2InvM, allWt);
    AnaUtil::fillHist1D ("zp1", zp1, allWt);
    AnaUtil::fillHist1D ("zp2", zp2, allWt);
    AnaUtil::fillHist1D ("zpMax", zpMax, allWt);
    AnaUtil::fillHist1D ("zpMin", zpMin, allWt);
    AnaUtil::fillHist1D ("HT", hT, allWt);
    AnaUtil::fillHist1D ("HTvec", hTvec, allWt);
    AnaUtil::fillHist1D ("HTMET", hTMET, allWt);
    AnaUtil::fillHist1D ("HTMETvec", hTMETvec, allWt);
    AnaUtil::fillHist1D ("HTfrac", hTfrac, allWt);
    AnaUtil::fillHist1D ("HTMETfrac", hTMETfrac, allWt);

    AnaUtil::fillHist1D ("j1j2DR", j1j2DR, allWt);
    AnaUtil::fillHist1D ("j1j2DEta", j1j2DEta, allWt);
    AnaUtil::fillHist1D ("j1j2DPhi", j1j2DPhi, allWt);

    AnaUtil::fillHist1D ("j1j3DR", j1j3DR, allWt);
    AnaUtil::fillHist1D ("j1j3DEta", j1j3DEta, allWt);
    AnaUtil::fillHist1D ("j1j3DPhi", j1j3DPhi, allWt);

    AnaUtil::fillHist1D ("j2j3DR", j2j3DR, allWt);
    AnaUtil::fillHist1D ("j2j3DEta", j2j3DEta, allWt);
    AnaUtil::fillHist1D ("j2j3DPhi", j2j3DPhi, allWt);

    AnaUtil::fillHist1D ("minJJDR", minJJDR, allWt);
    AnaUtil::fillHist1D ("minJJDPhi", minJJDPhi, allWt);
    AnaUtil::fillHist1D ("maxJJDR", maxJJDR, allWt);
    AnaUtil::fillHist1D ("maxJJDPhi", maxJJDPhi, allWt);


    AnaUtil::fillHist1D ("totCharge", totCharge, allWt);    
    AnaUtil::fillHist1D ("METovSQHT", METovSQHT, allWt);
    AnaUtil::fillHist1D ("METovSQLT", METovSQLT, allWt);
    AnaUtil::fillHist1D ("ST", ST, allWt);
    AnaUtil::fillHist1D ("LTovSqrtST", LTovSqrtST, allWt);
    AnaUtil::fillHist1D ("HTovSqrtST", HTovSqrtST, allWt);
    AnaUtil::fillHist1D ("STvec", STvec, allWt);
    AnaUtil::fillHist1D ("STMETvec", STMETvec, allWt);


    //------------------End of PreSelection--------------------//

    histf()->cd();

    //---------------------MVA Skimming---------------------//    
   
    if (skimObj_) {

      TreeVariables varList;
      
      varList.puwt       = puWt;
      varList.evwt       = evWt;
      varList.btagwt     = evt.btagWeight_CSVV2;
      varList.nLep       = LepCandList.size();
      varList.nJet       = tightJetListLV.size();
      varList.met        = met.pt;
      varList.jet1pt     = tightJetListLV[0].pt;
      varList.jet2pt     = tightJetListLV[1].pt;
      varList.HT         = hT;
      varList.HTvec      = hTvec;
      varList.HTMET      = hTMET; //new
      varList.HTMETvec   = hTMETvec; //new
      varList.HTfrac     = hTfrac; //new
      varList.HTMETfrac  = hTMETfrac;//new
      varList.jInvM      = jInvM;

      varList.j1j2DR     = j1j2DR;
      varList.j1j2DEta   = j1j2DEta;
      varList.j1j2DPhi   = j1j2DPhi;

      varList.j1j3DR     = j1j3DR;
      varList.j1j3DEta   = j1j3DEta;
      varList.j1j3DPhi   = j1j3DPhi;

      varList.j2j3DR     = j2j3DR;
      varList.j2j3DEta   = j2j3DEta;
      varList.j2j3DPhi   = j2j3DPhi;

      varList.minJJDR    = minJJDR;
      varList.minJJDPhi  = minJJDPhi;
      varList.maxJJDR    = maxJJDR;
      varList.maxJJDPhi  = maxJJDPhi;


      varList.lep1pt     = lep1.pt;
      varList.lep2pt     = lep2.pt;
      varList.LT         = lT;
      varList.LTMET      = LTMET;
      varList.LTMETvec   = LTMETvec;
      varList.LTMETfrac    = LTMETfrac; //new
      varList.LTMETfracInv = LTMETfracInv; //new
      varList.MT         = mT;
      varList.MTovLT     = mTovLt;
      varList.l1l2InvM   = l1l2InvM;
      varList.l1l2DR     = lepDR;
      varList.l1l2DEta   = lepDEta;
      varList.l1l2DPhi   = lepDPhi;

      varList.l1MetDPhi  = lep1MetDPhi;
      varList.l2MetDPhi  = lep2MetDPhi;
      varList.j1MetDPhi  = jet1MetDPhi;
      varList.j2MetDPhi  = jet2MetDPhi;
      varList.j3MetDPhi  = jet3MetDPhi; //new

      varList.j1l1DR     = j1l1DR;
      varList.j1l1DEta   = j1l1DEta;
      varList.j1l1DPhi   = j1l1DPhi;

      varList.j1l2DR     = j1l2DR;
      varList.j1l2DEta   = j1l2DEta;
      varList.j1l2DPhi   = j1l2DPhi;

      varList.j2l1DR     = j2l1DR;
      varList.j2l1DEta   = j2l1DEta;
      varList.j2l1DPhi   = j2l1DPhi;

      varList.j2l2DR     = j2l2DR;
      varList.j2l2DEta   = j2l2DEta;
      varList.j2l2DPhi   = j2l2DPhi;

      //new
      varList.j3l1DR     = j3l1DR;
      varList.j3l1DEta   = j3l1DEta;
      varList.j3l1DPhi   = j3l1DPhi;

      varList.j3l2DR     = j3l2DR;
      varList.j3l2DEta   = j3l2DEta;
      varList.j3l2DPhi   = j3l2DPhi;

      varList.maxJLDR   = maxJLDR;
      varList.minJLDR   = minJLDR;
      varList.maxJLDPhi = maxJLDPhi;
      varList.minJLDPhi = minJLDPhi;
      //new

      varList.minLepMetDPhi = minLepMetDPhi;
      varList.maxLepMetDPhi = maxLepMetDPhi;
      varList.minJetMetDPhi = minJetMetDPhi;
      varList.maxJetMetDPhi = maxJetMetDPhi;
      
      //new
      varList.METovSQHT   = METovSQHT;
      varList.METovSQLT   = METovSQLT;
      varList.ST          = ST;
      varList.LTovSqrtST  = LTovSqrtST;
      varList.HTovSqrtST  = HTovSqrtST;
      varList.STvec       = STvec;
      varList.STMETvec    = STMETvec;

      skimObj_->fill(varList);
    }
    
    //-------------------------MVA Evaluation------------------------//    
    
    histf()->cd();     
    double mvaOut = -999.9;
    
    if (_readMVA) {
      InputVariables varlist;
      
      //      varlist.met        = met.pt;
      varlist.HT         = hT;
      //varlist.HTvec      = hTvec;
      //      varlist.HTMETfrac  = hTMETfrac;//new
      varlist.j1j2DR     = j1j2DR;
      varlist.j1j3DR     = j1j3DR;
      //varlist.minJJDR    = minJJDR;
      varlist.maxJJDR    = maxJJDR;
      //varlist.minJJDPhi  = minJJDPhi;
      varlist.LT         = lT;
      //varlist.LTMETvec   = LTMETvec;
      //varlist.LTMETfracInv  = LTMETfracInv; //new
      varlist.l1l2InvM   = l1l2InvM;
      varlist.l1l2DR     = lepDR;
      varlist.l1MetDPhi  = lep1MetDPhi;
      varlist.j1MetDPhi  = jet1MetDPhi;
      varlist.j2MetDPhi  = jet2MetDPhi;
      // varlist.j3MetDPhi  = jet3MetDPhi; //new
      varlist.j1l1DPhi   = j1l1DPhi;
      varlist.minJLDPhi  = minJLDPhi;
      varlist.maxLepMetDPhi = maxLepMetDPhi;
      varlist.minJetMetDPhi = minJetMetDPhi;

      mvaOut = mvaObj_ -> evaluate(_MVAnetwork, varlist);
    }
    
    histf()->cd();
    histf()->cd("TMVAnalysis");
    
    AnaUtil::fillHist1D ("mvaOut", mvaOut, allWt);    
    
    if (nbJets > 0) continue;
    AnaUtil::fillHist1D("evtCutFlow", 12);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 12, allWt);

    AnaUtil::fillHist1D ("mvaOut_1", mvaOut, allWt);    

    if (mvaOut < 0.11) continue;
    AnaUtil::fillHist1D("evtCutFlow", 13);
    if (isMC()) AnaUtil::fillHist1D("evtCutFlowWt", 13, allWt);

    //PostBDT Session=====================>>>
    
    AnaUtil::fillHist1D ("ST_PBDT", ST, allWt);
    AnaUtil::fillHist1D ("jInvM_PBDT", jInvM, allWt);
    AnaUtil::fillHist1D ("j1j2InvM_PBDT", j1j2InvM, allWt);
    AnaUtil::fillHist1D ("w1w2mT_PBDT", mT, allWt);
    AnaUtil::fillHist1D ("w1w2mToverLT_PBDT", mTovLt, allWt);
    AnaUtil::fillHist1D ("LTMET_PBDT", LTMET, allWt);
    AnaUtil::fillHist1D ("nJets_PBDT", nJet, allWt);
    

    if (nJet == 3) {
      AnaUtil::fillHist1D ("evYield_nJetBinned", 3); 
      if (isMC()) AnaUtil::fillHist1D ("evYieldWt_nJetBinned", 3, allWt);
      AnaUtil::fillHist1D ("jInvM_nJ3", jInvM, allWt); 
    }
    else if (nJet == 4) {
      AnaUtil::fillHist1D ("evYield_nJetBinned", 4); 
      if (isMC()) AnaUtil::fillHist1D ("evYieldWt_nJetBinned", 4, allWt);
      AnaUtil::fillHist1D ("jInvM_nJ4", jInvM, allWt); 
    }
    else if (nJet > 4) {
      AnaUtil::fillHist1D ("evYield_nJetBinned", 5); 
      if (isMC()) AnaUtil::fillHist1D ("evYieldWt_nJetBinned", 5, allWt);
      AnaUtil::fillHist1D ("jInvM_nJ4+", jInvM, allWt); 
    }

    int l1flv = lep1.flavour;
    int l2flv = lep2.flavour;
    if (l1flv == 1 && l2flv == 1) { //mumu
      AnaUtil::fillHist1D ("evYield_lepFlvBinned", 0); 
      if (isMC()) AnaUtil::fillHist1D ("evYieldWt_lepFlvBinned", 0, allWt);
    }
    else if (l1flv == 2 && l2flv == 2) { //ee
      AnaUtil::fillHist1D ("evYield_lepFlvBinned", 1); 
      if (isMC()) AnaUtil::fillHist1D ("evYieldWt_lepFlvBinned", 1, allWt);
    }
    else if (l1flv != l2flv) { //emu
      AnaUtil::fillHist1D ("evYield_lepFlvBinned", 2); 
      if (isMC()) AnaUtil::fillHist1D ("evYieldWt_lepFlvBinned", 2, allWt);
    }

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
      "2) nLeptons >= 2    : ",
      "3) lep1pT > 25 GeV  : ",
      "4) lep2pT > 15 GeV  : ",
      "5) Same Charge      : ",
      "6) no Z Cand        : ",
      "7) nJet >= 3        : ",
      "8) jet1Pt > 60 GeV  : ",
      "9) jet2Pt > 35 GeV  : ",
      "10) lT > 65 GeV(NA) : ",
      "11) MET > 50 GeV    : ",
      "12)no bJet(PostBDT) : ",
      "13)BDT >= 0.11      : "
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

#ifdef SKIP_DUPLICATE_FROM_FILE
  if (skipDuplicate_ && !eventFilelist_.empty()) {
    eventIdStore_.clear();
    for (const auto& f: eventFilelist_) {
      cout << ">>> Reading file: " << f << endl;
      ifstream fin(f, std::ios::in);
      if (!fin) {
	cerr << "Input file: " << f << " could not be opened!" << endl;
	continue;
      }
      char buf[BUF_SIZE];
      while (fin.getline(buf, BUF_SIZE, '\n')) {  // Pops off the newline character
	string line(buf);
	if (line.empty()) continue;   
	
	// enable '#' and '//' style comments
	if (line.substr(0,1) == "#" || line.substr(0,2) == "//") continue;
    
	// Split the line into words
	vector<string> tokens;
	AnaUtil::tokenize(line, tokens);
	assert(tokens.size() > 2);
	string key = tokens.at(0) + "-" + tokens.at(1) + "-" + tokens.at(2);
	eventIdStore_.insert({key, 1});
      }
      // Close the file
      fin.close();
    }
    cout << ">>> Total events present: " << eventIdStore_.size() << endl;
  }  
#endif

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
