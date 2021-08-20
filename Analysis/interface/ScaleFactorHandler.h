#ifndef __SCALEFACTORHANDLER__H
#define __SCALEFACTORHANDLER__H

#include <string>
#include "TH2F.h"
#include "TH2D.h"

class ScaleFactorHandler {
public: 
  ScaleFactorHandler() {}
  virtual ~ScaleFactorHandler() {}
  bool openRootFiles();
  double getIdSF(const std::string& IdType, float pt, float eta, const std::string& Flav) const;
  double getIsoSF(const std::string& IsoType, float pt, float eta, const std::string& Flav) const;

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
};
#endif
