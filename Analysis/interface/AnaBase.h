#ifndef __AnaBase__hh
#define __AnaBase__hh

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <map>
#include <fstream>

#include "TLorentzVector.h"
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include "PhysicsObjects.h"
#include "AnaUtil.h"
#include "ScaleFactorHandler.h"

using uint = unsigned int;
using ulong = unsigned long;

// REQUIRED, most probably to solve circular dependence problem!!!
class TChain;
class TFile;

template <class T>
class PtComparator {
public:
  bool operator()(const T &a, const T &b) const {
    return a.pt > b.pt;
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

// Specially for ZCandidate sorting
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
  bool branchFound(const std::string& b, const std::string& tname="Events");
  int getEntries() const;

  const std::unique_ptr<TChain>& chain() const {return chain_;}
  std::unique_ptr<TChain>& chain() {return chain_;}

  const std::unique_ptr<TChain>& chainRun() const {return chainRun_;}
  std::unique_ptr<TChain>& chainRun() {return chainRun_;}

  const std::unique_ptr<TFile>& histf() const {return histf_;}
  std::unique_ptr<TFile>& histf() {return histf_;}

  const std::unique_ptr<TTreeReader>& treeReader() const {return treeReader_;}
  std::unique_ptr<TTreeReader>& treeReader() {return treeReader_;}

  const std::unique_ptr<TTreeReader>& treeReaderRun() const {return treeReaderRun_;}
  std::unique_ptr<TTreeReader>& treeReaderRun() {return treeReaderRun_;}

  int nEvents() const {return nEvents_;}
  int firstEvent() const {return firstEvt_;}
  int lastEvent() const {return lastEvt_;}

  std::ofstream& fLog() {return fLog_;}
  std::ofstream& evLog() {return evLog_;}
  std::ofstream& selEvLog() {return selEvLog_;}

  const std::ofstream& fLog() const {return fLog_;}
  const std::ofstream& evLog() const {return evLog_;}
  const std::ofstream& selEvLog() const {return selEvLog_;}

  bool readGenInfo() const {return readGenInfo_;}
  bool isMC() const {return isMC_;}
  bool isSignal() const {return isSignal_;}
  int logOption() const {return logOption_;}
  double lumiWt(double evtWeightSum=-1, bool verbose=false) const;
  int getEra() const {return era_;}
  std::string getDatasetName() {return dataset_;}

  const std::map<std::string, double>& lumiWtMap() const {return AnaUtil::cutMap(hmap_, "lumiWtList");}
  const std::map<std::string, double>& vtxCutMap() const {return AnaUtil::cutMap(hmap_, "vtxCutList");}
  const std::map<std::string, double>& muonCutMap() const {return AnaUtil::cutMap(hmap_, "muonCutList");}
  const std::map<std::string, double>& photonCutMap() const {return AnaUtil::cutMap(hmap_, "photonCutList");}
  const std::map<std::string, double>& electronCutMap() const {return AnaUtil::cutMap(hmap_, "electronCutList");}
  const std::map<std::string, double>& tauCutMap() const {return AnaUtil::cutMap(hmap_, "tauCutList");}
  const std::map<std::string, double>& jetCutMap() const {return AnaUtil::cutMap(hmap_, "jetCutList");}
  const std::map<std::string, double>& fatJetCutMap() const {return AnaUtil::cutMap(hmap_, "fatJetCutList");}
  const std::map<std::string, double>& evselCutMap() const {return AnaUtil::cutMap(hmap_, "evselCutList");}

  const std::unordered_map<std::string, int>& eventIdMap() const {return eventIdMap_;}

  const std::vector<std::string> getSingleMuonHLTpaths() const {return singleMuonHltPathList_;}
  const std::vector<std::string> getDoubleMuonHLTpaths() const {return doubleMuonHltPathList_;}
  const std::vector<std::string> getSingleElectronHLTpaths() const {return singleElectronHltPathList_;}
  const std::vector<std::string> getDoubleEgHLTpaths() const {return doubleEgHltPathList_;}
  const std::vector<std::string> getMuonEgHLTpaths() const {return muonEgHltPathList_;}
  const std::vector<std::string> getSingleMuonHLTForFakepaths() const {return singleMuonHltForFakePathList_;}
  const std::vector<std::string> getSingleElectronHLTForFakepaths() const {return singleElectronHltForFakePathList_;}

  const std::vector<std::unique_ptr<TTreeReaderValue<bool>>>& getDoubleMuonHLTptrs() const {return doubleMuonHltPtrList_;}
  const std::vector<std::unique_ptr<TTreeReaderValue<bool>>>& getSingleMuonHLTptrs() const {return singleMuonHltPtrList_;}
  const std::vector<std::unique_ptr<TTreeReaderValue<bool>>>& getDoubleEgHLTptrs() const {return doubleEgHltPtrList_;}
  const std::vector<std::unique_ptr<TTreeReaderValue<bool>>>& getSingleElectronHLTptrs() const {return singleElectronHltPtrList_;}
  const std::vector<std::unique_ptr<TTreeReaderValue<bool>>>& getMuonEgHLTptrs() const {return muonEgHltPtrList_;}
  const std::vector<std::unique_ptr<TTreeReaderValue<bool>>>& getSingleElectronHLTForFakeptrs() const {return singleElectronHltForFakePtrList_;}
  const std::vector<std::unique_ptr<TTreeReaderValue<bool>>>& getSingleMuonHLTForFakeptrs() const {return singleMuonHltForFakePtrList_;}

  void setHltPtrList(const std::vector<std::string>& hltPathList, std::vector<std::unique_ptr<TTreeReaderValue<bool>>>& hltPtrList);

  static std::vector<bool> getHLTscores(const std::vector<std::unique_ptr<TTreeReaderValue<bool>>>& ptrs);
  static bool isDuplicate(bool passDoubleMuonHLT, 
		   bool passDoubleEgHLT, 
		   bool passMuonEgHLT, 
		   bool passSingleMuonHLT, 
		   bool passSingleEleHLT,
		   const std::string& dataset);
  static bool isTriggered(const std::vector<std::string>& paths, const std::vector<bool>& scores);

  const ScaleFactorHandler& SFHandler() const {return SFHandler_;}
  ScaleFactorHandler& SFHandler() {return SFHandler_;}

private:
  std::unique_ptr<TChain> chain_;                // chain contains a list of root files containing the Events tree
  std::unique_ptr<TChain> chainRun_;             // chain contains a list of root files containing the Runs tree
  std::unique_ptr<TTreeReader> treeReader_;      //!pointer to the analyzed TTree::Events
  std::unique_ptr<TTreeReader> treeReaderRun_;   //!pointer to the analyzed TTree::Runs
  std::unique_ptr<TFile> histf_;                 // The output file with histograms

  int nEvents_;
  std::ofstream fLog_;   
  std::ofstream evLog_;   
  std::ofstream selEvLog_;   

  bool isMC_ {false};
  bool isSignal_ {false};
  bool readGenInfo_ {false};
  std::vector<std::string> fileList_;
  std::vector<std::string> singleMuonHltPathList_;
  std::vector<std::string> singleElectronHltPathList_;
  std::vector<std::string> doubleMuonHltPathList_;
  std::vector<std::string> doubleEgHltPathList_;
  std::vector<std::string> muonEgHltPathList_;
  std::vector<std::string> singleMuonHltForFakePathList_;
  std::vector<std::string> singleElectronHltForFakePathList_;

  std::vector<std::unique_ptr<TTreeReaderValue<bool>>> doubleMuonHltPtrList_;
  std::vector<std::unique_ptr<TTreeReaderValue<bool>>> singleMuonHltPtrList_;
  std::vector<std::unique_ptr<TTreeReaderValue<bool>>> doubleEgHltPtrList_;
  std::vector<std::unique_ptr<TTreeReaderValue<bool>>> singleElectronHltPtrList_;
  std::vector<std::unique_ptr<TTreeReaderValue<bool>>> muonEgHltPtrList_;
  std::vector<std::unique_ptr<TTreeReaderValue<bool>>> singleMuonHltForFakePtrList_;
  std::vector<std::unique_ptr<TTreeReaderValue<bool>>> singleElectronHltForFakePtrList_;

  int logOption_ {0};

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

public:
  ScaleFactorHandler SFHandler_;

  // Required Branches

  //----------------------------------TTree::Runs-----------------------------------//
  std::unique_ptr<TTreeReaderValue<double>> genEventSumw_;
  std::unique_ptr<TTreeReaderValue<double>> genEventSumw2_;

  //----------------------------------TTree::Events---------------------------------//
  std::unique_ptr<TTreeReaderValue<uint>> run_;
  std::unique_ptr<TTreeReaderValue<unsigned long long>> event_;
  std::unique_ptr<TTreeReaderValue<uint>> lumis_;

  // weights
  std::unique_ptr<TTreeReaderValue<float>> PU_Weight;
  std::unique_ptr<TTreeReaderValue<float>> PU_WeightUp;
  std::unique_ptr<TTreeReaderValue<float>> PU_WeightDown;
  std::unique_ptr<TTreeReaderValue<float>> genEvWt;
  std::unique_ptr<TTreeReaderValue<float>> btagWeight_CSVV2;
  std::unique_ptr<TTreeReaderValue<float>> btagWeight_CMVA;

  // primary vertex
  std::unique_ptr<TTreeReaderValue<float>> PV_ndf;
  std::unique_ptr<TTreeReaderValue<float>> PV_xPos;
  std::unique_ptr<TTreeReaderValue<float>> PV_yPos;
  std::unique_ptr<TTreeReaderValue<float>> PV_zPos;
  std::unique_ptr<TTreeReaderValue<float>> PV_chi2;
  std::unique_ptr<TTreeReaderValue<float>> PV_score;
  std::unique_ptr<TTreeReaderValue<int>> nPV;
  std::unique_ptr<TTreeReaderValue<int>> nGoodPV;

  // Muon
  // https://cms-nanoaod-integration.web.cern.ch/integration/master-cmsswmaster/mc102X_doc.html#Muon
  // https://github.com/cms-nanoAOD/cmssw/blob/master-cmsswmaster/PhysicsTools/NanoAOD/python/muons_cff.py
  // https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2 :: for RunII only. UL is under processing
  std::unique_ptr<TTreeReaderValue<uint>> nMuon;
  std::unique_ptr<TTreeReaderArray<float>> Muon_pt;
  std::unique_ptr<TTreeReaderArray<float>> Muon_corrpt;
  std::unique_ptr<TTreeReaderArray<float>> Muon_eta;
  std::unique_ptr<TTreeReaderArray<float>> Muon_phi;
  std::unique_ptr<TTreeReaderArray<int>> Muon_charge;
  std::unique_ptr<TTreeReaderArray<float>> Muon_mass;
  std::unique_ptr<TTreeReaderArray<float>> Muon_sip3d;
  std::unique_ptr<TTreeReaderArray<int>> Muon_jetIdx;
  std::unique_ptr<TTreeReaderArray<bool>> Muon_LooseId;
  std::unique_ptr<TTreeReaderArray<bool>> Muon_MediumId;
  std::unique_ptr<TTreeReaderArray<bool>> Muon_TightId;
  std::unique_ptr<TTreeReaderArray<unsigned char>> Muon_mvaId;
  std::unique_ptr<TTreeReaderArray<unsigned char>> Muon_highPtId; 
  std::unique_ptr<TTreeReaderArray<float>> Muon_pfRelIso03_all;
  std::unique_ptr<TTreeReaderArray<float>> Muon_pfRelIso04_all;
  std::unique_ptr<TTreeReaderArray<float>> Muon_pfRelIso03_chg;
  std::unique_ptr<TTreeReaderArray<int>> Muon_genPartIdx;
  std::unique_ptr<TTreeReaderArray<unsigned char>> Muon_genPartFlv; 
  std::unique_ptr<TTreeReaderArray<bool>> Muon_isGlobal;
  std::unique_ptr<TTreeReaderArray<bool>> Muon_isPFcand;
  std::unique_ptr<TTreeReaderArray<bool>> Muon_isTracker;
  std::unique_ptr<TTreeReaderArray<float>> Muon_dxy;
  std::unique_ptr<TTreeReaderArray<float>> Muon_dz;
  std::unique_ptr<TTreeReaderArray<int>> Muon_tightCharge; //Tight charge criterion using pterr/pt of muonBestTrack (0:fail, 2:pass)
  std::unique_ptr<TTreeReaderArray<float>> Muon_miniPFRelIso_all;
  
  // Electron
  // https://cms-nanoaod-integration.web.cern.ch/integration/master-cmsswmaster/mc102X_doc.html#Electron
  // https://twiki.cern.ch/twiki/bin/view/CMS/EgammaUL2016To2018
  // https://github.com/cms-nanoAOD/cmssw/blob/master-cmsswmaster/PhysicsTools/NanoAOD/python/electrons_cff.py
  // https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Working_points_for_2016_data_for
  std::unique_ptr<TTreeReaderValue<uint>> nElectron;
  std::unique_ptr<TTreeReaderArray<float>> Electron_pt;
  std::unique_ptr<TTreeReaderArray<float>> Electron_eta;
  std::unique_ptr<TTreeReaderArray<float>> Electron_phi;
  std::unique_ptr<TTreeReaderArray<int>> Electron_charge;
  std::unique_ptr<TTreeReaderArray<float>> Electron_mass;
  std::unique_ptr<TTreeReaderArray<float>> Electron_sip3d;
  std::unique_ptr<TTreeReaderArray<int>> Electron_jetIdx;
  std::unique_ptr<TTreeReaderArray<int>> Electron_phoIdx;
  std::unique_ptr<TTreeReaderArray<float>> Electron_pfRelIso03_all;
  std::unique_ptr<TTreeReaderArray<float>> Electron_pfRelIso03_chg;
  std::unique_ptr<TTreeReaderArray<bool>> Electron_mvaFall17V2Iso_WP80;
  std::unique_ptr<TTreeReaderArray<bool>> Electron_mvaFall17V2Iso_WP90;
  std::unique_ptr<TTreeReaderArray<bool>> Electron_mvaFall17V1Iso_WP90;
  std::unique_ptr<TTreeReaderArray<bool>> Electron_mvaFall17V2noIso_WP80;
  std::unique_ptr<TTreeReaderArray<bool>> Electron_mvaFall17V2noIso_WP90;
  std::unique_ptr<TTreeReaderArray<int>> Electron_genPartIdx;
  std::unique_ptr<TTreeReaderArray<unsigned char>> Electron_genPartFlv; 
  std::unique_ptr<TTreeReaderArray<float>> Electron_dxy;
  std::unique_ptr<TTreeReaderArray<float>> Electron_dz;
  std::unique_ptr<TTreeReaderArray<bool>> Electron_mvaFall17V2noIso_WPL;
  std::unique_ptr<TTreeReaderArray<bool>> Electron_mvaFall17V1noIso_WPL;
  std::unique_ptr<TTreeReaderArray<unsigned char>> Electron_lostHits;
  std::unique_ptr<TTreeReaderArray<bool>> Electron_convVeto;
  std::unique_ptr<TTreeReaderArray<float>> Electron_deltaEtaSC;

  // https://twiki.cern.ch/twiki/bin/view/CMS/MultivariatePhotonIdentificationRun2#MVA_recipes_for_2016_data_and_Sp
  // Photon

  // Jet
  // https://github.com/cms-nanoAOD/cmssw/blob/master-cmsswmaster/PhysicsTools/NanoAOD/python/jets_cff.py
  std::unique_ptr<TTreeReaderValue<uint>> nJet;
  std::unique_ptr<TTreeReaderArray<float>> Jet_pt;
  std::unique_ptr<TTreeReaderArray<float>> Jet_nomPt;
  std::unique_ptr<TTreeReaderArray<float>> Jet_eta;
  std::unique_ptr<TTreeReaderArray<float>> Jet_phi;
  std::unique_ptr<TTreeReaderArray<int>> Jet_nConstituents;
  std::unique_ptr<TTreeReaderArray<int>> Jet_nMuons;
  std::unique_ptr<TTreeReaderArray<int>> Jet_nElectrons;
  std::unique_ptr<TTreeReaderArray<int>> Jet_puId;
  std::unique_ptr<TTreeReaderArray<int>> Jet_muIdx1;
  std::unique_ptr<TTreeReaderArray<int>> Jet_muIdx2;
  std::unique_ptr<TTreeReaderArray<int>> Jet_elIdx1;
  std::unique_ptr<TTreeReaderArray<int>> Jet_elIdx2;
  std::unique_ptr<TTreeReaderArray<int>> Jet_jetId; //https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD#Jets
  std::unique_ptr<TTreeReaderArray<float>> Jet_qgl; //Quark vs Gluon likelihood discriminator
  std::unique_ptr<TTreeReaderArray<float>> Jet_neHEF; //neutral Hadron Energy Fraction
  std::unique_ptr<TTreeReaderArray<float>> Jet_neEmEF; //neutral Electromagnetic Energy Fraction
  std::unique_ptr<TTreeReaderArray<float>> Jet_chHEF; //charged Hadron Energy Fraction
  std::unique_ptr<TTreeReaderArray<float>> Jet_chEmEF; //charged Electromagnetic Energy Fraction
  std::unique_ptr<TTreeReaderArray<float>> Jet_area;
  std::unique_ptr<TTreeReaderArray<float>> Jet_mass;
  std::unique_ptr<TTreeReaderArray<float>> Jet_nomMass;
  std::unique_ptr<TTreeReaderArray<float>> Jet_btagDeepFlavB;
  std::unique_ptr<TTreeReaderArray<float>> Jet_btagCSVV2;

  // MET
  std::unique_ptr<TTreeReaderValue<float>> Met_pt;
  std::unique_ptr<TTreeReaderValue<float>> Met_nomPt;
  std::unique_ptr<TTreeReaderValue<float>> Met_phi;
  std::unique_ptr<TTreeReaderValue<float>> Met_nomPhi;
  std::unique_ptr<TTreeReaderValue<float>> Met_significance;
  std::unique_ptr<TTreeReaderValue<float>> Met_sumEt;

  // FatJet
  std::unique_ptr<TTreeReaderValue<uint>> nFatJet;
  std::unique_ptr<TTreeReaderArray<float>> FatJet_pt;
  std::unique_ptr<TTreeReaderArray<float>> FatJet_eta;
  std::unique_ptr<TTreeReaderArray<float>> FatJet_phi;
  std::unique_ptr<TTreeReaderArray<float>> FatJet_mass;
  std::unique_ptr<TTreeReaderArray<int>> FatJet_jetId;
  std::unique_ptr<TTreeReaderArray<float>> FatJet_btagDeepB;
  std::unique_ptr<TTreeReaderArray<float>> FatJet_btagCSVV2; // pfCombinedInclusiveSecondaryVertexV2 b-tag discriminator (aka CSVV2)
  std::unique_ptr<TTreeReaderArray<float>> FatJet_msoftdrop;
  std::unique_ptr<TTreeReaderArray<float>> FatJet_n2b1;
  std::unique_ptr<TTreeReaderArray<float>> FatJet_n3b1;
  std::unique_ptr<TTreeReaderArray<unsigned char>> FatJet_nBHadrons;
  std::unique_ptr<TTreeReaderArray<unsigned char>> FatJet_nCHadrons;
  std::unique_ptr<TTreeReaderArray<int>> FatJet_hadronFlavour;
  std::unique_ptr<TTreeReaderArray<float>> FatJet_rawFactor;
  std::unique_ptr<TTreeReaderArray<int>> FatJet_subJetIdx1;
  std::unique_ptr<TTreeReaderArray<int>> FatJet_subJetIdx2;
  std::unique_ptr<TTreeReaderArray<float>> FatJet_tau1;
  std::unique_ptr<TTreeReaderArray<float>> FatJet_tau2;
  std::unique_ptr<TTreeReaderArray<float>> FatJet_tau3;
  std::unique_ptr<TTreeReaderArray<float>> FatJet_tau4;
  std::unique_ptr<TTreeReaderArray<float>> FatJet_deepTag_WvsQCD;
  std::unique_ptr<TTreeReaderArray<float>> FatJet_deepTag_ZvsQCD;
  std::unique_ptr<TTreeReaderArray<float>> FatJet_deepTag_TvsQCD;
  std::unique_ptr<TTreeReaderArray<float>> FatJet_deepTagMD_WvsQCD;
  std::unique_ptr<TTreeReaderArray<float>> FatJet_deepTagMD_ZvsQCD;
  std::unique_ptr<TTreeReaderArray<float>> FatJet_deepTagMD_TvsQCD;
  std::unique_ptr<TTreeReaderArray<int>> FatJet_electronIdx3SJ;
  std::unique_ptr<TTreeReaderArray<int>> FatJet_muonIdx3SJ;

  // SubJet
  std::unique_ptr<TTreeReaderValue<uint>> nSubJet;
  std::unique_ptr<TTreeReaderArray<float>> SubJet_pt;
  std::unique_ptr<TTreeReaderArray<float>> SubJet_eta;
  std::unique_ptr<TTreeReaderArray<float>> SubJet_phi;
  std::unique_ptr<TTreeReaderArray<float>> SubJet_mass;
  std::unique_ptr<TTreeReaderArray<float>> SubJet_btagDeepB;
  std::unique_ptr<TTreeReaderArray<float>> SubJet_rawFactor;

  // Tau
  std::unique_ptr<TTreeReaderValue<uint>> nTau;
  std::unique_ptr<TTreeReaderArray<float>> Tau_pt;
  std::unique_ptr<TTreeReaderArray<float>> Tau_eta;
  std::unique_ptr<TTreeReaderArray<float>> Tau_phi;
  std::unique_ptr<TTreeReaderArray<float>> Tau_mass;
  std::unique_ptr<TTreeReaderArray<float>> Tau_dxy;
  std::unique_ptr<TTreeReaderArray<float>> Tau_dz;
  std::unique_ptr<TTreeReaderArray<int>> Tau_charge;
  std::unique_ptr<TTreeReaderArray<int>> Tau_decayMode;
  std::unique_ptr<TTreeReaderArray<bool>> Tau_idDecayMode;
  std::unique_ptr<TTreeReaderArray<bool>> Tau_idDecayModeNewDMs;
  std::unique_ptr<TTreeReaderArray<int>> Tau_jetIdx;                 // index of the associated jet (-1 if none)
  std::unique_ptr<TTreeReaderArray<unsigned char>> Tau_idAntiEle;    // Anti-electron MVA discriminator V6: 
                                         // bitmask 1 = VLoose, 2 = Loose, 4 = Medium, 8 = Tight, 16 = VTight
  std::unique_ptr<TTreeReaderArray<unsigned char>> Tau_idAntiMu;     // Anti-muon discriminator V3: : bitmask 1 = Loose, 2 = Tight
  std::unique_ptr<TTreeReaderArray<unsigned char>> Tau_idMVAoldDM;   // IsolationMVArun2v1DBoldDMwLT ID working point (2015): 
                                         // bitmask 1 = VLoose, 2 = Loose, 4 = Medium, 8 = Tight, 16 = VTight, 32 = VVTight
  std::unique_ptr<TTreeReaderArray<unsigned char>> Tau_idDeepTau2017v2VSjet;
  std::unique_ptr<TTreeReaderArray<unsigned char>> Tau_idDeepTau2017v2VSmu;
  std::unique_ptr<TTreeReaderArray<unsigned char>> Tau_idDeepTau2017v2VSe;
  std::unique_ptr<TTreeReaderArray<unsigned char>> Tau_idDeepTau2017v2p1VSjet;
  std::unique_ptr<TTreeReaderArray<unsigned char>> Tau_idDeepTau2017v2p1VSmu;
  std::unique_ptr<TTreeReaderArray<unsigned char>> Tau_idDeepTau2017v2p1VSe;

  // GenParticle
  std::unique_ptr<TTreeReaderValue<uint>> nGenPart;
  std::unique_ptr<TTreeReaderArray<float>> GenPart_pt;
  std::unique_ptr<TTreeReaderArray<float>> GenPart_eta;
  std::unique_ptr<TTreeReaderArray<float>> GenPart_phi;
  std::unique_ptr<TTreeReaderArray<float>> GenPart_mass;
  std::unique_ptr<TTreeReaderArray<int>> GenPart_motherIdx;
  std::unique_ptr<TTreeReaderArray<int>> GenPart_pdgId;
  std::unique_ptr<TTreeReaderArray<int>> GenPart_status;  
  std::unique_ptr<TTreeReaderArray<int>> GenPart_statusFlags;  

  // LHE, LHEPart
  std::unique_ptr<TTreeReaderValue<uint>> nLHEPart;
  std::unique_ptr<TTreeReaderArray<float>> LHEPart_pt;
  std::unique_ptr<TTreeReaderArray<float>> LHEPart_eta;
  std::unique_ptr<TTreeReaderArray<float>> LHEPart_phi;
  std::unique_ptr<TTreeReaderArray<float>> LHEPart_mass;
  std::unique_ptr<TTreeReaderArray<int>> LHEPart_pdgId;
  std::unique_ptr<TTreeReaderValue<unsigned char>> LHEnJets;
};
#endif
