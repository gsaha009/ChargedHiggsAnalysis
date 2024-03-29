#ifndef __PhysicsObjSelector__hh
#define __PhysicsObjSelector__hh

#define NEL(x) (sizeof((x))/sizeof((x)[0]))

#include <string>
#include <vector>

#include "PhysicsObjects.h"
#include "AnaBase.h"
#include "LeptonCand.h"
#include "ZCandidate.h"

template <typename T>
void packLeptons(const std::vector<T>& lepList, std::vector<LeptonCand>& candList, bool ismc=false) {
  for (const auto& v: lepList) {
    LeptonCand lc;
    lc.index   = v.index;
    lc.pt      = v.pt;
    lc.eta     = v.eta;
    lc.SCeta   = v.SCeta;
    lc.phi     = v.phi;
    lc.mass    = v.mass;
    lc.charge  = v.charge;
    lc.flavour = (typeid(v) == typeid(vhtm::Muon)) ? 1 : 2;
    if (ismc) lc.genFlv  = v.genFlv;

    candList.push_back(lc);
  }
}
// create unique lepton pair combination giving a Z and add the candidate into a vector
template <typename T>
void ZSelector(const std::vector<T>& lepList, std::vector<ZSpace::ZCandidate>& candList) {
  for (unsigned int i = 0; i < lepList.size(); ++i) {
    const auto& ip = lepList[i];

    TLorentzVector lep1P4(AnaUtil::getP4(ip));
    for (unsigned int j = i+1; j < lepList.size(); ++j) {
      const auto& jp = lepList[j];
	
      // require opposite charges
      if ((ip.charge * jp.charge) > 0) continue; 
	
      TLorentzVector lep2P4(AnaUtil::getP4(jp));
      if (lep1P4.DeltaR(lep2P4) < 0.02) continue;

      ZSpace::ZCandidate ztmp;
      if (typeid(jp) == typeid(vhtm::Muon)) {
	ztmp.flavour = static_cast<int>(ZSpace::ZType::mumu);
      }
      else if (typeid(jp) == typeid(vhtm::Electron)) {
	ztmp.flavour = static_cast<int>(ZSpace::ZType::ee);
      }
      else 
	ztmp.flavour = static_cast<int>(ZSpace::ZType::unkwn);
	
      ztmp.l1Index = i;
      ztmp.l1P4 = lep1P4;
      ztmp.l1Charge = ip.charge;
	
      ztmp.l2Index = j;
      ztmp.l2P4 = lep2P4;
      ztmp.l2Charge = jp.charge;
	
      TLorentzVector p4 = lep1P4 + lep2P4;
      double Zmass = p4.M();
      ztmp.p4 = p4;
      ztmp.mass = Zmass;
      ztmp.massDiff = std::fabs(Zmass - ZSpace::MZnominal);
      ztmp.dEtall = ztmp.l1P4.Eta() - ztmp.l2P4.Eta();
      ztmp.dPhill = TVector2::Phi_mpi_pi(ztmp.l1P4.Phi() - ztmp.l2P4.Phi());
      ztmp.dRll   = ztmp.l1P4.DeltaR(ztmp.l2P4);

      candList.push_back(ztmp);
    }
  }
}

class PhysicsObjSelector: public AnaBase {
 public:
  PhysicsObjSelector();
  virtual ~PhysicsObjSelector() {
  }
  virtual void eventLoop() = 0;  // the main analysis
  virtual bool beginJob() override;
  virtual void endJob() override;
  virtual bool readJob(const std::string& jobFile, int& nFiles) override;
  virtual void bookHistograms();
  void objectEfficiency();
  bool jetLeptonCleaning(const vhtm::Jet& jet, double minDR=0.4) const;
  bool fatJetLeptonCleaning(const vhtm::FatJet& jet, double minDR=0.8) const;
  bool tauLeptonCleaning(const vhtm::Tau& tau, double minDR=0.3) const;
  bool thisElectronIsMuon(const vhtm::Electron& ele, bool VsLooseMuons, bool VsTightMuons, double minDR=0.3) const;  

  //-----Inline functions to get the object collections-----
  const std::vector<vhtm::Event>& getEventList() const {return eventList_;}
  const std::vector<vhtm::MET>& getMETList() const {return metList_;}

  // Jet Collections
  const std::vector<vhtm::Jet>& getPreSelJetList() const {return preSelJetList_;}
  const std::vector<vhtm::Jet>& getCleanJetList() const {return leptonCleanJetList_;}  
  const std::vector<vhtm::Jet>& getAk8CleanJetList() const {return leptonCleanJetListOutsideAk8_;}  
  const std::vector<vhtm::Jet>& getLooseBJetList() const {return looseBJetList_;}  
  const std::vector<vhtm::Jet>& getBJetList() const {return bJetList_;}  

  // FatJet Collections
  const std::vector<vhtm::FatJet>& getFatJetList() const {return fatJetList_;}
  const std::vector<vhtm::FatJet>& getCleanFatJetList() const {return cleanFatJetList_;}
  const std::vector<vhtm::FatJet>& getBTaggedFatJetList() const {return bTaggedFatJetList_;}
  
  // SubJet Collections
  const std::vector<vhtm::SubJet>& getSubJetList() const {return subJetList_;}

  // Tau Collections
  const std::vector<vhtm::Tau>& getTauList() const {return tauList_;}
  const std::vector<vhtm::Tau>& getLepCleanTauList() const {return leptonCleanTauList_;}

  // Muon Collections
  const std::vector<vhtm::Muon>& getPreSelMuList() const {return preSelMuList_;}
  const std::vector<vhtm::Muon>& getFakeableMuList() const {return fakeableMuList_;}
  const std::vector<vhtm::Muon>& getTightMuList() const {return tightMuList_;}

  // Electron Collection
  const std::vector<vhtm::Electron>& getPreSelEleList() const {return preSelEleList_;}
  const std::vector<vhtm::Electron>& getFakeableEleList() const {return fakeableEleList_;}
  const std::vector<vhtm::Electron>& getTightEleList() const {return tightEleList_;}

  // GenParticle Collections
  const std::vector<vhtm::GenParticle>& getGenList() const {return genParticleList_;}
  const std::vector<vhtm::LHEParticle>& getLHEList() const {return lheParticleList_;}

  //---------------------------------------------------------  
  bool findEventInfo();
  bool findGenPartInfo();
  bool findLHEPartInfo();
  void findObjects(); 
  double getGenSumW();

  void muonSelector();
  void tauSelector();
  void electronSelector();
  void metSelector();
  void jetSelector();
  void fatJetSelector();
  void subJetSelector();

  void dumpEvent(int evNo, std::ostream& os=std::cout) const;

  void clear();

 private:
  std::vector<vhtm::Event> eventList_;
  std::vector<vhtm::MET> metList_;
  std::vector<vhtm::Jet> preSelJetList_, 
    leptonCleanJetList_, 
    leptonCleanJetListOutsideAk8_, 
    looseBJetList_, 
    bJetList_; 
  std::vector<vhtm::FatJet> fatJetList_, cleanFatJetList_, bTaggedFatJetList_;
  std::vector<vhtm::SubJet> subJetList_;
  std::vector<vhtm::Tau> tauList_, leptonCleanTauList_;
  std::vector<vhtm::Muon> preSelMuList_, fakeableMuList_, tightMuList_;
  std::vector<vhtm::Electron> preSelEleList_, fakeableEleList_, tightEleList_; 
  std::vector<vhtm::GenParticle> genParticleList_;
  std::vector<vhtm::LHEParticle> lheParticleList_;

  bool searchedEle_{false}, 
       searchedMu_{false}, 
       searchedJet_{false}, 
       searchedFatJet_{false};
};
#endif
