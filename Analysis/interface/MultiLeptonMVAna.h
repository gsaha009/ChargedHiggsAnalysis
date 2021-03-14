#ifndef __MultiLeptonMVAna__hh
#define __MultiLeptonMVAna__hh

#define NEL(x) (sizeof((x))/sizeof((x)[0]))

#include "configana.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <bitset>

#include "TLorentzVector.h"
#include "TVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

#include "PhysicsObjects.h"
#include "PhysicsObjSelector.h"
#include "AnaUtil.h"
#include "MVASkim.h"
#include "MVAnalysis.h"
#include "LeptonCand.h"

//see AN-13-178
template<typename LVector1, typename LVector2>
  double Calculate_MT(const LVector1& lepton_p4, const LVector2& met_p4)
{
  const double delta_phi = TVector2::Phi_mpi_pi(lepton_p4.Phi() -
						met_p4.Phi());
  return std::sqrt( 2.0 * lepton_p4.Pt() * met_p4.Pt() * ( 1.0 -
							   std::cos(delta_phi) ) );
}

template<typename LVector1, typename LVector2, typename LVector3>
  double Calculate_TotalMT_Doubt(const LVector1& lepton1_p4, const LVector2&
			   lepton2_p4, const LVector3& met_p4)
{
  const double mt_1 = Calculate_MT(lepton1_p4, met_p4);
  const double mt_2 = Calculate_MT(lepton2_p4, met_p4);
  const double mt_ll = Calculate_MT(lepton1_p4, lepton2_p4);
  return std::sqrt(std::pow(mt_1, 2) + std::pow(mt_2, 2) +
		   std::pow(mt_ll, 2));
}
template<typename LVector1, typename LVector2, typename LVector3>
  double Calculate_TotalMT(const LVector1& lepton1_p4, const LVector2&
			   lepton2_p4, const LVector3& met_p4)
{
  const double mt_1 = Calculate_MT(lepton1_p4, met_p4);
  const double mt_2 = Calculate_MT(lepton2_p4, met_p4);
  return std::sqrt(std::pow(mt_1, 2) + std::pow(mt_2, 2));
}


class MultiLeptonMVAna : public PhysicsObjSelector {
    
public:
  enum class EventType {
    unkwn = -1, mmem = 0, eeem
  };
  MultiLeptonMVAna();
  virtual ~MultiLeptonMVAna();
  

  virtual void eventLoop();  // the main analysis 
  virtual bool beginJob();
  virtual void endJob();
  virtual bool readJob(const std::string& jobFile, int& nFiles);
  virtual void printJob(std::ostream& os=std::cout) const;
  virtual void bookHistograms();
  virtual void closeFiles();
    
  void clearLists();
  bool hasZcandidate(const std::vector<LeptonCand>& LepColl);
  bool hasLowMassResonance(const std::vector<LeptonCand>& LepColl);
private:
  double lumiFac {0.0};
  double evtWeightSum_ {0.0};
  bool dumpGenInfo_ {false};
  bool useEventList_ {false};
  bool skipDuplicate_ {false};
  bool selectPM_ {false};
  int nMEPartons_ {-1};
  ofstream syncDumpf_;
  std::string dumpFilename_ {"syncDumpFile.txt"};
  int dumpEventCount_ {0};
  std::vector<std::string> eventFilelist_;
  std::unordered_map<std::string, int> eventIdStore_;

  std::unique_ptr<MVAnalysis> mvaObj_ {nullptr};
  std::unique_ptr<MVASkim> skimObj_ {nullptr};
  bool _createMVATree {false};
  bool _readMVA {false};
  std::string _mvaInputFile {""};  
  std::string _MVAnetwork {""};
  std::string _MVAxmlFile {""};
};
#endif
