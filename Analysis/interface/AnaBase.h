#ifndef __AnaBase__hh
#define __AnaBase__hh

#define NEL(x) (sizeof((x))/sizeof((x)[0]))

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <map>
#include <algorithm>

#include "TLorentzVector.h"
#include "TVector.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TVector3.h"
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include "PhysicsObjects.h"
#include "AnaUtil.h"

using uint = unsigned int;
using ulong = unsigned long;

// REQUIRED, most probably to solve circular dependence problem!!!
class TChain;
class TFile;
using std::ofstream;
using std::ifstream;

template <class T>
class PtComparator {
public:
  bool operator()(const T &a, const T &b) const {
    return a.pt > b.pt;
  }
};

template <class T>
class PtComparatorTL {
public:
  bool operator()(const T &a, const T &b) const {
    return a.Pt() > b.Pt();
  }
};

template <class T>
class PtComparatorLep {
public:
  bool operator()(const T &a, const T &b) const {
    return a.lPt > b.lPt;
  }
};

template <class T>
class MassComparator {
public:
  bool operator()(const T &a, const T &b) const {
    TLorentzVector l1,l2;
    l1.SetPtEtaPhiE(a.pt, a.eta, a.phi, a.energy);
    l2.SetPtEtaPhiE(b.pt, b.eta, b.phi, b.energy);
    return l1.M() > l2.M();
  }
};

//Specially for ZCandidate sorting
template <class T>
class MassDiffComparator {
public:
  bool operator()(const T &a, const T &b) const {
    return a.massDiff < b.massDiff;
  }
};

typedef struct  
{
  bool verbose;
  bool usesbit;
  bool printselected;
} Options;

class AnaBase {
    
public:

  AnaBase();
  virtual ~AnaBase();
    
  virtual void eventLoop() = 0;  // the main analysis
  virtual bool beginJob();
  virtual void endJob();
  virtual bool readJob(const std::string& jobFile, int& nFiles);
  virtual void printJob(std::ostream& os=std::cout) const;
  virtual bool init();

  virtual bool openFiles();
  virtual void closeFiles();
  virtual void closeHistFile();

  int setInputFile(const std::string& fname);
  bool branchFound(const std::string& b);
  int getEntries() const;
  //  bool readPileUpHist(bool verbose=false);
  void clearEvent();
  void enableBranches();
  int getEntry(int lflag) const;

  std::unique_ptr<TFile>& histf() {return histf_;}
  //std::unique_ptr<TFile>& fakehistf() {return fakehistf_;}
  TTreeReader* treeReader() {return treeReader_;}
  TTreeReader* treeReaderRun() {return treeReaderRun_;}
  TChain* chain() {return chain_;}
  TChain* chainRun() {return chainRun_;}
  int nEvents() const {return nEvents_;}
  int firstEvent() const {return firstEvt_;}
  int lastEvent() const {return lastEvt_;}

  ofstream& fLog() {return fLog_;}
  ofstream& evLog() {return evLog_;}
  ofstream& selEvLog() {return selEvLog_;}

  //  const ofstream& fLog() const {return fLog_;}
  //  const ofstream& evLog() const {return evLog_;}
  //  const ofstream& selEvLog() const {return selEvLog_;}

  bool readGenInfo() const {return readGenInfo_;}
  bool isMC() const {return isMC_;}
  bool isSignal() const {return isSignal_;}
  int logOption() const {return logOption_;}
  bool useTrigger() const {return useTrigger_;}
  bool useLumiWt() const {return useLumiWt_;}
  double lumiWt(double evtWeightSum=-1, bool verbose=false) const;
  bool usePUWt() const {return usePUWt_;}
  bool useTrueNInt() const {return useTrueNInt_;}
  int getEra() const {return era_;}
  std::string getDatasetName() {return dataset_;}
  bool openScaleFactorRootFiles();
  double getIdSF(std::string IdType, float pt, float eta, std::string Flav) const;
  double getIsoSF(std::string IsoType, float pt, float eta, std::string Flav) const;

  const std::map<std::string, double>& lumiWtMap() const {return AnaUtil::cutMap(hmap_, "lumiWtList");}
  const std::map<std::string, double>& vtxCutMap() const {return AnaUtil::cutMap(hmap_, "vtxCutList");}
  const std::map<std::string, double>& muonCutMap() const {return AnaUtil::cutMap(hmap_, "muonCutList");}
  const std::map<std::string, double>& photonCutMap() const {return AnaUtil::cutMap(hmap_, "photonCutList");}
  const std::map<std::string, double>& packedPFCandidateCutMap() const {return AnaUtil::cutMap(hmap_, "packedPFCandidateCutList");}
  const std::map<std::string, double>& electronCutMap() const {return AnaUtil::cutMap(hmap_, "electronCutList");}
  const std::map<std::string, double>& tauCutMap() const {return AnaUtil::cutMap(hmap_, "tauCutList");}
  const std::map<std::string, double>& jetCutMap() const {return AnaUtil::cutMap(hmap_, "jetCutList");}
  const std::map<std::string, double>& fatJetCutMap() const {return AnaUtil::cutMap(hmap_, "fatJetCutList");}
  const std::map<std::string, double>& evselCutMap() const {return AnaUtil::cutMap(hmap_, "evselCutList");}

  const std::unordered_map<std::string, int>& eventIdMap() const {return eventIdMap_;}
  int bunchCrossing() const {return bunchCrossing_;}

  const std::vector< std::string > getSingleMuonHLTpaths() const {return singleMuonHltPathList_;}
  const std::vector< std::string > getDoubleMuonHLTpaths() const {return doubleMuonHltPathList_;}
  const std::vector< std::string > getSingleElectronHLTpaths() const {return singleElectronHltPathList_;}
  const std::vector< std::string > getDoubleEgHLTpaths() const {return doubleEgHltPathList_;}
  const std::vector< std::string > getMuonEgHLTpaths() const {return muonEgHltPathList_;}

  std::vector< TTreeReaderValue<bool>* > getDoubleMuonHLTptrs() {return doubleMuonHltPtrList_;}
  std::vector< TTreeReaderValue<bool>* > getSingleMuonHLTptrs() {return singleMuonHltPtrList_;}
  std::vector< TTreeReaderValue<bool>* > getDoubleEgHLTptrs() {return doubleEgHltPtrList_;}
  std::vector< TTreeReaderValue<bool>* > getSingleElectronHLTptrs() {return singleElectronHltPtrList_;}
  std::vector< TTreeReaderValue<bool>* > getMuonEgHLTptrs() {return muonEgHltPtrList_;}

private:
  //  std::unique_ptr<TChain> chain_;      // chain contains a list of root files containing the same tree
  std::unique_ptr<TFile> histf_;       // The output file with histograms
  //std::unique_ptr<TFile> fakehistf_;       // The output file with histograms
  TChain* chain_; 
  TChain* chainRun_;  
  std::vector<std::string> brList_;
  
  int nEvents_;
  ofstream fLog_;   
  ofstream evLog_;   
  ofstream selEvLog_;   

  bool isMC_ {false};
  bool isSignal_ {false};
  bool readGenInfo_ {false};
  bool readTrigObject_ {false};
  bool readPFObject_ {false};
  std::vector<std::string> fileList_;
  std::vector<std::string> singleMuonHltPathList_;
  std::vector<std::string> singleElectronHltPathList_;
  std::vector<std::string> doubleMuonHltPathList_;
  std::vector<std::string> doubleEgHltPathList_;
  std::vector<std::string> muonEgHltPathList_;

  std::vector< TTreeReaderValue<bool>* >doubleMuonHltPtrList_;
  std::vector< TTreeReaderValue<bool>* >singleMuonHltPtrList_;
  std::vector< TTreeReaderValue<bool>* >doubleEgHltPtrList_;
  std::vector< TTreeReaderValue<bool>* >singleElectronHltPtrList_;
  std::vector< TTreeReaderValue<bool>* >muonEgHltPtrList_;

  int logOption_ {0};
  bool useTrigger_ {false};
  bool useLumiWt_ {false};
  bool usePUWt_ {false};
  std::string puHistFile_ {"./reweightFunctionFall11.root"};
  std::string puHistogram_ {"pileup"};
  bool useTrueNInt_ {true};
  int bunchCrossing_ {25};

  std::string evtWtSum_;
  std::string dataset_ {"bla"};
  std::string histFile_ {"default.root"};
  //std::string fakehistFile_ {"fakedefault.root"};
  std::string logFile_ {"default.out"};
  std::string evFile_ {"events.out"};
  std::string selEvFile_ {"selected_events.out"};
  int maxEvt_ {0};
  int nFiles_ {0};
  int firstEvt_ {-1};
  int lastEvt_ {-1};
  int era_ {0};

  std::map<std::string, std::map<std::string, double>> hmap_;
  std::unordered_map<std::string, int> eventIdMap_;

  TTreeReader* treeReader_;   //!pointer to the analyzed TTree::Events 
  TTreeReader* treeReaderRun_;   //!pointer to the analyzed TTree::Runs 
  // SF related
  std::string muonIdSFRootFile_ {"default.root"};
  std::string looseMuonIdSFhistName_ {"hist"};
  TH2D* looseMuonIdSFhist_ {nullptr};
  std::string medMuonIdSFhistName_ {"hist"};
  TH2D* medMuonIdSFhist_ {nullptr};
  std::string tightMuonIdSFhistName_ {"hist"};
  TH2D* tightMuonIdSFhist_ {nullptr};

  std::string electronLooseIdSFRootFile_ {"default.root"};
  std::string looseEleIdSFhistName_ {"hist"};
  TH2F* looseEleIdSFhist_ {nullptr};
  std::string electronTightIdSFRootFile_ {"default.root"};
  std::string tightEleIdSFhistName_ {"hist"};
  TH2F* tightEleIdSFhist_ {nullptr};

  std::string muonTightIsoSFRootFile_ {"default.root"};
  std::string tightMuIsoSFhistName_ {"hist"};
  TH2D* tightMuIsoSFhist_ {nullptr};

 public:
  // Required Branches

  //----------------------------------TTree::Runs-----------------------------------//
  TTreeReaderValue< double >* genEventSumw_;
  TTreeReaderValue< double >* genEventSumw2_;

  //----------------------------------TTree::Events---------------------------------//
  TTreeReaderValue< unsigned int >* run_;
  TTreeReaderValue< unsigned long long>* event_;
  TTreeReaderValue< unsigned int >* lumis_;

  // weights
  TTreeReaderValue< float >* PU_Weight;
  TTreeReaderValue< float >* PU_WeightUp;
  TTreeReaderValue< float >* PU_WeightDown;
  TTreeReaderValue< float >* genEvWt;
  TTreeReaderValue< float >* btagWeight_CSVV2;
  TTreeReaderValue< float >* btagWeight_CMVA;

  // primary vertex
  TTreeReaderValue< float >* PV_ndf;
  TTreeReaderValue< float >* PV_xPos;
  TTreeReaderValue< float >* PV_yPos;
  TTreeReaderValue< float >* PV_zPos;
  TTreeReaderValue< float >* PV_chi2;
  TTreeReaderValue< float >* PV_score;
  TTreeReaderValue< int >* nPV;
  TTreeReaderValue< int >* nGoodPV;

  //Muon
  // https://cms-nanoaod-integration.web.cern.ch/integration/master-cmsswmaster/mc102X_doc.html#Muon
  // https://github.com/cms-nanoAOD/cmssw/blob/master-cmsswmaster/PhysicsTools/NanoAOD/python/muons_cff.py
  // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2 :: for RunII only. UL is under processing
  TTreeReaderValue< unsigned int >* nMuon;
  TTreeReaderArray< float >* Muon_pt;
  TTreeReaderArray< float >* Muon_corrpt;
  TTreeReaderArray< float >* Muon_eta;
  TTreeReaderArray< float >* Muon_phi;
  TTreeReaderArray< int >* Muon_charge;
  TTreeReaderArray< float >* Muon_mass;
  TTreeReaderArray< float >* Muon_sip3d;
  TTreeReaderArray< int >* Muon_jetIdx;
  TTreeReaderArray< bool >* Muon_LooseId;
  TTreeReaderArray< bool >* Muon_MediumId;
  TTreeReaderArray< bool >* Muon_TightId;
  TTreeReaderArray< unsigned char >* Muon_mvaId;
  TTreeReaderArray< unsigned char >* Muon_highPtId; 
  TTreeReaderArray< float >* Muon_pfRelIso03_all;
  TTreeReaderArray< float >* Muon_pfRelIso04_all;
  TTreeReaderArray< float >* Muon_pfRelIso03_chg;
  TTreeReaderArray< int >* Muon_genPartIdx;
  TTreeReaderArray< unsigned char >* Muon_genPartFlv; 
  TTreeReaderArray< bool >* Muon_isGlobal;
  TTreeReaderArray< bool >* Muon_isPFcand;
  TTreeReaderArray< bool >* Muon_isTracker;
  TTreeReaderArray< float >* Muon_dxy;
  TTreeReaderArray< float >* Muon_dz;
  TTreeReaderArray< int >* Muon_tightCharge; 
  TTreeReaderArray< float >* Muon_miniPFRelIso_all;
  
  // Electron
  // https://cms-nanoaod-integration.web.cern.ch/integration/master-cmsswmaster/mc102X_doc.html#Electron
  // https://twiki.cern.ch/twiki/bin/view/CMS/EgammaUL2016To2018
  // https://github.com/cms-nanoAOD/cmssw/blob/master-cmsswmaster/PhysicsTools/NanoAOD/python/electrons_cff.py
  // https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Working_points_for_2016_data_for
  TTreeReaderValue< unsigned int >* nElectron;
  TTreeReaderArray< float >* Electron_pt;
  TTreeReaderArray< float >* Electron_eta;
  TTreeReaderArray< float >* Electron_phi;
  TTreeReaderArray< int >* Electron_charge;
  TTreeReaderArray< float >* Electron_mass;
  TTreeReaderArray< float >* Electron_sip3d;
  TTreeReaderArray< int >* Electron_jetIdx;
  TTreeReaderArray< int >* Electron_phoIdx;
  TTreeReaderArray< float >* Electron_pfRelIso03_all;
  TTreeReaderArray< float >* Electron_pfRelIso03_chg;
  TTreeReaderArray< bool >* Electron_mvaFall17V2Iso_WP80;
  TTreeReaderArray< bool >* Electron_mvaFall17V2Iso_WP90;
  TTreeReaderArray< bool >* Electron_mvaFall17V1Iso_WP90;
  TTreeReaderArray< bool >* Electron_mvaFall17V2noIso_WP80;
  TTreeReaderArray< bool >* Electron_mvaFall17V2noIso_WP90;
  TTreeReaderArray< int >* Electron_genPartIdx;
  TTreeReaderArray< unsigned char >* Electron_genPartFlv; 
  TTreeReaderArray< float >* Electron_dxy;
  TTreeReaderArray< float >* Electron_dz;
  TTreeReaderArray< bool >* Electron_mvaFall17V2noIso_WPL;
  TTreeReaderArray< bool >* Electron_mvaFall17V1noIso_WPL;
  TTreeReaderArray< unsigned char >*Electron_lostHits;
  TTreeReaderArray< bool >*Electron_convVeto;
  TTreeReaderArray< float >*Electron_deltaEtaSC;

  // Jet
  // https://github.com/cms-nanoAOD/cmssw/blob/master-cmsswmaster/PhysicsTools/NanoAOD/python/jets_cff.py
  TTreeReaderValue< unsigned int >* nJet;
  TTreeReaderArray< float >* Jet_pt;
  TTreeReaderArray< float >* Jet_nomPt;
  TTreeReaderArray< float >* Jet_eta;
  TTreeReaderArray< float >* Jet_phi;
  TTreeReaderArray< int >* Jet_nConstituents;
  TTreeReaderArray< int >* Jet_nMuons;
  TTreeReaderArray< int >* Jet_nElectrons;
  TTreeReaderArray< int >* Jet_puId;
  TTreeReaderArray< int >* Jet_muIdx1;
  TTreeReaderArray< int >* Jet_muIdx2;
  TTreeReaderArray< int >* Jet_elIdx1;
  TTreeReaderArray< int >* Jet_elIdx2;
  TTreeReaderArray< int >* Jet_jetId;
  TTreeReaderArray< float >* Jet_qgl; //Quark vs Gluon likelihood discriminator
  TTreeReaderArray< float >* Jet_neHEF; //neutral Hadron Energy Fraction
  TTreeReaderArray< float >* Jet_neEmEF; //neutral Electromagnetic Energy Fraction
  TTreeReaderArray< float >* Jet_chHEF; //charged Hadron Energy Fraction
  TTreeReaderArray< float >* Jet_chEmEF; //charged Electromagnetic Energy Fraction
  TTreeReaderArray< float >* Jet_area;
  TTreeReaderArray< float >* Jet_mass;
  TTreeReaderArray< float >* Jet_nomMass;
  TTreeReaderArray< float >* Jet_btagDeepFlavB;
  TTreeReaderArray< float >* Jet_btagCSVV2;

  //MET
  TTreeReaderValue< float >* Met_pt;
  TTreeReaderValue< float >* Met_nomPt;
  TTreeReaderValue< float >* Met_phi;
  TTreeReaderValue< float >* Met_nomPhi;
  TTreeReaderValue< float >* Met_significance;
  TTreeReaderValue< float >* Met_sumEt;



  //FatJet
  TTreeReaderValue< unsigned int >* nFatJet;
  TTreeReaderArray< float >* FatJet_pt;
  TTreeReaderArray< float >* FatJet_eta;
  TTreeReaderArray< float >* FatJet_phi;
  TTreeReaderArray< float >* FatJet_mass;
  TTreeReaderArray< int >* FatJet_jetId;
  TTreeReaderArray< float >* FatJet_btagDeepB;
  TTreeReaderArray< float >* FatJet_btagCSVV2;//pfCombinedInclusiveSecondaryVertexV2 b-tag discriminator (aka CSVV2)
  TTreeReaderArray< float >* FatJet_msoftdrop;
  TTreeReaderArray< float >* FatJet_n2b1;
  TTreeReaderArray< float >* FatJet_n3b1;
  TTreeReaderArray< unsigned char >* FatJet_nBHadrons;
  TTreeReaderArray< unsigned char >* FatJet_nCHadrons;
  TTreeReaderArray< int >* FatJet_hadronFlavour;
  TTreeReaderArray< float >* FatJet_rawFactor;
  TTreeReaderArray< int >* FatJet_subJetIdx1;
  TTreeReaderArray< int >* FatJet_subJetIdx2;
  TTreeReaderArray< float >* FatJet_tau1;
  TTreeReaderArray< float >* FatJet_tau2;
  TTreeReaderArray< float >* FatJet_tau3;
  TTreeReaderArray< float >* FatJet_tau4;
  TTreeReaderArray< float >* FatJet_deepTag_WvsQCD;
  TTreeReaderArray< float >* FatJet_deepTag_ZvsQCD;
  TTreeReaderArray< float >* FatJet_deepTag_TvsQCD;
  TTreeReaderArray< float >* FatJet_deepTagMD_WvsQCD;
  TTreeReaderArray< float >* FatJet_deepTagMD_ZvsQCD;
  TTreeReaderArray< float >* FatJet_deepTagMD_TvsQCD;
  TTreeReaderArray< int >* FatJet_electronIdx3SJ;
  TTreeReaderArray< int >* FatJet_muonIdx3SJ;

  //SubJet
  TTreeReaderValue< unsigned int >* nSubJet;
  TTreeReaderArray< float >* SubJet_pt;
  TTreeReaderArray< float >* SubJet_eta;
  TTreeReaderArray< float >* SubJet_phi;
  TTreeReaderArray< float >* SubJet_mass;
  TTreeReaderArray< float >* SubJet_btagDeepB;
  TTreeReaderArray< float >* SubJet_rawFactor;

  //Tau
  TTreeReaderValue< unsigned int >* nTau;
  TTreeReaderArray< float >* Tau_pt;
  TTreeReaderArray< float >* Tau_eta;
  TTreeReaderArray< float >* Tau_phi;
  TTreeReaderArray< float >* Tau_mass;
  TTreeReaderArray< float >* Tau_dxy;
  TTreeReaderArray< float >* Tau_dz;
  TTreeReaderArray< int >* Tau_charge;
  TTreeReaderArray< int >* Tau_decayMode;
  TTreeReaderArray< bool >* Tau_idDecayMode;
  TTreeReaderArray< bool >* Tau_idDecayModeNewDMs;
  TTreeReaderArray< int >* Tau_jetIdx; //index of the associated jet (-1 if none)
  TTreeReaderArray< unsigned char >* Tau_idAntiEle;//Anti-electron MVA discriminator V6: 
                                                  //bitmask 1 = VLoose, 2 = Loose, 4 = Medium, 8 = Tight, 16 = VTight
  TTreeReaderArray< unsigned char >* Tau_idAntiMu;//Anti-muon discriminator V3: : bitmask 1 = Loose, 2 = Tight
  TTreeReaderArray< unsigned char >* Tau_idMVAoldDM;//IsolationMVArun2v1DBoldDMwLT ID working point (2015): 
                                                    //bitmask 1 = VLoose, 2 = Loose, 4 = Medium, 8 = Tight, 16 = VTight, 32 = VVTight
  TTreeReaderArray< unsigned char >* Tau_idDeepTau2017v2VSjet;
  TTreeReaderArray< unsigned char >* Tau_idDeepTau2017v2VSe;
  TTreeReaderArray< unsigned char >* Tau_idDeepTau2017v2VSmu;
  TTreeReaderArray< unsigned char >* Tau_idDeepTau2017v2p1VSjet;
  TTreeReaderArray< unsigned char >* Tau_idDeepTau2017v2p1VSe;
  TTreeReaderArray< unsigned char >* Tau_idDeepTau2017v2p1VSmu;



  //GenParticle
  TTreeReaderValue< unsigned int >* nGenPart;
  TTreeReaderArray< float >* GenPart_pt;
  TTreeReaderArray< float >* GenPart_eta;
  TTreeReaderArray< float >* GenPart_phi;
  TTreeReaderArray< float >* GenPart_mass;
  TTreeReaderArray< int >* GenPart_motherIdx;
  TTreeReaderArray< int >* GenPart_pdgId;
  TTreeReaderArray< int >* GenPart_status;  
  TTreeReaderArray< int >* GenPart_statusFlags;  

  //LHE, LHEPart
  TTreeReaderValue< unsigned int >* nLHEPart;
  TTreeReaderArray< float >* LHEPart_pt;
  TTreeReaderArray< float >* LHEPart_eta;
  TTreeReaderArray< float >* LHEPart_phi;
  TTreeReaderArray< float >* LHEPart_mass;
  TTreeReaderArray< int >* LHEPart_pdgId;
  TTreeReaderValue< unsigned char >* LHEnJets;

  //userDefined 
  TTreeReaderValue< char >* HLT_SingleMu;
  TTreeReaderValue< char >* HLT_SingleEle;
  TTreeReaderValue< char >* HLT_DoubleMu;
  TTreeReaderValue< char >* HLT_DoubleEle;

  TTreeReaderValue< int >* evType; //:0,1,2=Same Charged Leptons :: -1=Opposite charged Leptons  
};
#endif
