#include <memory>
#include <iostream>

#include "TSystem.h"
#include "TFile.h"
#include "ScaleFactorHandler.h"

using namespace std;
bool ScaleFactorHandler::openRootFiles() {
  // MUON POG ID SCALE FACTOR
  const char* fname_MuonIdSF = gSystem->ExpandPathName(muonIdSFRootFile_.c_str());
  if (gSystem->AccessPathName(fname_MuonIdSF)) {
    cerr << ">>> Warning: File <<" << muonIdSFRootFile_ << ">> not found!!" << endl;
    return false;
  }
  std::unique_ptr<TFile> file_MuonIdSF = std::make_unique<TFile>(fname_MuonIdSF);

  // --- Loose Muon ID SF
  file_MuonIdSF->GetObject(looseMuonIdSFhistName_.c_str(), looseMuonIdSFhist_);
  if (!looseMuonIdSFhist_) {
    cerr << ">>> Warning: Histogram <<" << looseMuonIdSFhistName_ << ">> not found!!" << endl;
    return false;
  }
  looseMuonIdSFhist_->SetDirectory(0);

  // --- Medium Muon ID SF
  file_MuonIdSF->GetObject(medMuonIdSFhistName_.c_str(), medMuonIdSFhist_);
  if (!medMuonIdSFhist_) {
    cerr << ">>> Warning: Histogram <<" << medMuonIdSFhistName_ << ">> not found!!" << endl;
    return false;
  }
  medMuonIdSFhist_->SetDirectory(0);

  // --- Tight Muon ID SF
  file_MuonIdSF->GetObject(tightMuonIdSFhistName_.c_str(), tightMuonIdSFhist_);
  if (!tightMuonIdSFhist_) {
    cerr << ">>> Warning: Histogram <<" << tightMuonIdSFhist_ << ">> not found!!" << endl;
    return false;
  }
  tightMuonIdSFhist_->SetDirectory(0);

  file_MuonIdSF->Close();

  // Electron POG ID SCALE FACTOR
  const char* fname_EleLooseIdSF = gSystem->ExpandPathName(electronLooseIdSFRootFile_.c_str());
  if (gSystem->AccessPathName(fname_EleLooseIdSF)) {
    cerr << ">>> Warning: File <<" << electronLooseIdSFRootFile_ << ">> not found!!" << endl;
    return false;
  }
  std::unique_ptr<TFile> file_EleLooseIdSF = std::make_unique<TFile>(fname_EleLooseIdSF);

  file_EleLooseIdSF->GetObject(looseEleIdSFhistName_.c_str(), looseEleIdSFhist_);
  if (!looseEleIdSFhist_) {
    cerr << ">>> Warning: Histogram <<" << looseEleIdSFhist_ << ">> not found!!" << endl;
    return false;
  }
  looseEleIdSFhist_->SetDirectory(0);

  file_EleLooseIdSF ->Close();

  const char* fname_EleTightIdSF = gSystem->ExpandPathName(electronTightIdSFRootFile_.c_str());
  if (gSystem->AccessPathName(fname_EleTightIdSF)) {
    cerr << ">>> Warning: File <<" << electronTightIdSFRootFile_ << ">> not found!!" << endl;
    return false;
  }
  std::unique_ptr<TFile> file_EleTightIdSF = std::make_unique<TFile>(fname_EleTightIdSF);

  file_EleTightIdSF->GetObject(tightEleIdSFhistName_.c_str(), tightEleIdSFhist_);
  if (!tightEleIdSFhist_) {
    cerr << ">>> Warning: Histogram <<" << tightEleIdSFhist_ << ">> not found!!" << endl;
    return false;
  }
  tightEleIdSFhist_->SetDirectory(0);

  file_EleTightIdSF->Close();

  // Iso Scale factor
  const char* fname_MuonTightIsoSF = gSystem->ExpandPathName(muonTightIsoSFRootFile_.c_str());
  if (gSystem->AccessPathName(fname_MuonTightIsoSF)) {
    cerr << ">>> Warning: File <<" << muonTightIsoSFRootFile_ << ">> not found!!" << endl;
    return false;
  }
  std::unique_ptr<TFile> file_MuTightIsoSF = std::make_unique<TFile>(fname_MuonTightIsoSF);

  file_MuTightIsoSF->GetObject(tightMuIsoSFhistName_.c_str(), tightMuIsoSFhist_);
  if (!tightMuIsoSFhist_) {
    cerr << ">>> Warning: Histogram <<" << tightMuIsoSFhist_ << ">> not found!!" << endl;
    return false;
  }
  tightMuIsoSFhist_->SetDirectory(0);

  file_MuTightIsoSF->Close();

  return true;
}
double ScaleFactorHandler::getIdSF(const std::string& IdType, float pt, float eta, const std::string& Flav) const {
  double SF = 0.0;
  if (Flav == "Muon") {
    if (IdType == "Loose") {
      if (pt < looseMuonIdSFhist_->GetXaxis()->GetXmax()) {
	int binX = looseMuonIdSFhist_->GetXaxis()->FindBin(pt);
	int binY = looseMuonIdSFhist_->GetYaxis()->FindBin(std::abs(eta));
	SF = looseMuonIdSFhist_->GetBinContent(binX, binY);
      }
      else SF = 1.0;
    }
    else if (IdType == "Medium") {
      if (pt < medMuonIdSFhist_->GetXaxis()->GetXmax()) {
	int binX = medMuonIdSFhist_->GetXaxis()->FindBin(pt);
	int binY = medMuonIdSFhist_->GetYaxis()->FindBin(std::abs(eta));
	SF = medMuonIdSFhist_->GetBinContent(binX, binY);
      }
      else SF = 1.0;
    }
    else if (IdType == "Tight") {
      if (pt < tightMuonIdSFhist_->GetXaxis()->GetXmax()) {
	int binX = tightMuonIdSFhist_->GetXaxis()->FindBin(pt);
	int binY = tightMuonIdSFhist_->GetYaxis()->FindBin(std::abs(eta));
	SF = tightMuonIdSFhist_->GetBinContent(binX, binY);
      }
    }
  }
  else if (Flav == "Electron") {
    if (IdType == "Loose") {
      if (pt < looseEleIdSFhist_->GetYaxis()->GetXmax()) {
	int binX = looseEleIdSFhist_->GetXaxis()->FindBin(std::abs(eta));
	int binY = looseEleIdSFhist_->GetYaxis()->FindBin(pt);
	SF = looseEleIdSFhist_->GetBinContent(binX, binY);
      }
      else SF = 1.00;
    }
    else if (IdType == "Tight") {
      if (pt < tightEleIdSFhist_->GetYaxis()->GetXmax()) {
	int binX = tightEleIdSFhist_->GetXaxis()->FindBin(std::abs(eta));
	int binY = tightEleIdSFhist_->GetYaxis()->FindBin(pt);
	SF = tightEleIdSFhist_->GetBinContent(binX, binY);
      }
      else SF = 1.0;
    }
  }

  return SF;
}
double ScaleFactorHandler::getIsoSF(const std::string& IsoType, float pt, float eta, const std::string& Flav) const {
  double SF = 0.0;
  if (Flav == "Muon") {
    if (IsoType == "Tight") {
      if (pt < tightMuIsoSFhist_->GetXaxis()->GetXmax()) {
        int binX = tightMuIsoSFhist_->GetXaxis()->FindBin(pt);
        int binY = tightMuIsoSFhist_->GetYaxis()->FindBin(std::abs(eta));
        SF = tightMuIsoSFhist_->GetBinContent(binX, binY);
      }
      else SF = 1.0;
    }
  }
  return SF;
}
