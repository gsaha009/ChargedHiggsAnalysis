#ifndef __MultiLeptonMVAna_DY__hh
#define __MultiLeptonMVAna_DY__hh

#define NEL(x) (sizeof((x))/sizeof((x)[0]))

#include "configana.h"

#include <memory>
#include <string>
#include <vector>

#include "PhysicsObjSelector.h"
#include "AnaUtil.h"
#include "MVASkim.h"
#include "MVAnalysis.h"
#include "LeptonCand.h"
#include "HistBooker.h"

class MultiLeptonMVAna_DY : public PhysicsObjSelector {
    
public:
  enum class EventType {
    unkwn = -1, mmem = 0, eeem
  };
  MultiLeptonMVAna_DY();
  virtual ~MultiLeptonMVAna_DY();  

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
  bool isPrompt(const LeptonCand& lep);

private:
  double lumiFac {0.0};
  double evtWeightSum_ {0.0};
  bool dumpGenInfo_ {false};
  bool useEventList_ {false};
  int dumpEventCount_ {0};
  std::vector<std::string> eventFilelist_;

  std::unique_ptr<MVAnalysis> mvaObj_ {nullptr};
  std::unique_ptr<MVASkim> skimObj_ {nullptr};

  bool _createMVATree {false};
  bool _readMVA {false};
  std::string _mvaInputFile {""};  
  std::string _MVAnetwork {""};
  std::string _MVAxmlFile {""};
};
#endif
