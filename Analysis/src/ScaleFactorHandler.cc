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

  // FakeRates
  const char* fname_FR = gSystem->ExpandPathName(FRRootFile_.c_str());
  if (gSystem->AccessPathName(fname_FR)) {
    cerr << ">>> Warning: File <<" << FRRootFile_ << ">> not found!!" << endl;
    return false;
  }
  std::unique_ptr<TFile> file_FR = std::make_unique<TFile>(fname_FR);

  // Muon FakeRate hist
  file_FR -> GetObject(muonFRhistName_.c_str(), muonFRhist_);
  if (!muonFRhist_) {
    cerr << ">>> Warning: Histogram <<" << muonFRhistName_ << ">> not found!!" << endl;
    return false;
  }
  muonFRhist_->SetDirectory(0);

  // Electron FakeRate hist
  file_FR -> GetObject(electronFRhistName_.c_str(), electronFRhist_);
  if (!electronFRhist_) {
    cerr << ">>> Warning: Histogram <<" << electronFRhistName_ << ">> not found!!" << endl;
    return false;
  }
  electronFRhist_->SetDirectory(0);

  file_FR ->Close();

  return true;
}
double ScaleFactorHandler::getIdSF(const std::string& IdType, float pt, float eta, const std::string& Flav) const {
  double SF = 0.0;
  if (Flav == "Muon") {
    if (IdType == "Loose") {
      if (pt >= looseMuonIdSFhist_-> GetXaxis()->GetXmax()) {
	int binX = looseMuonIdSFhist_->GetXaxis()->GetLast();
	int binY = looseMuonIdSFhist_->GetYaxis()->FindBin(std::abs(eta));
	SF = looseMuonIdSFhist_->GetBinContent(binX, binY);
      }
      else {
	int binX = looseMuonIdSFhist_->GetXaxis()->FindBin(pt);
	int binY = looseMuonIdSFhist_->GetYaxis()->FindBin(std::abs(eta));
	SF = looseMuonIdSFhist_->GetBinContent(binX, binY);
      }
    }
    else if (IdType == "Medium") {
      if (pt >= medMuonIdSFhist_->GetXaxis()->GetXmax()) {
	int binX = medMuonIdSFhist_->GetXaxis()->GetLast();
	int binY = medMuonIdSFhist_->GetYaxis()->FindBin(std::abs(eta));
	SF = medMuonIdSFhist_->GetBinContent(binX, binY);
      }
      else {
	int binX = medMuonIdSFhist_->GetXaxis()->FindBin(pt);
	int binY = medMuonIdSFhist_->GetYaxis()->FindBin(std::abs(eta));
	SF = medMuonIdSFhist_->GetBinContent(binX, binY);
      }
    }
    else if (IdType == "Tight") {
      if (pt >= tightMuonIdSFhist_->GetXaxis()->GetXmax()) {
	int binX = tightMuonIdSFhist_->GetXaxis()->GetLast();
	int binY = tightMuonIdSFhist_->GetYaxis()->FindBin(std::abs(eta));
	SF = tightMuonIdSFhist_->GetBinContent(binX, binY);
      }
      else {
	int binX = tightMuonIdSFhist_->GetXaxis()->FindBin(pt);
	int binY = tightMuonIdSFhist_->GetYaxis()->FindBin(std::abs(eta));
	SF = tightMuonIdSFhist_->GetBinContent(binX, binY);
      } 
    }
  }
  else if (Flav == "Electron") {
    if (IdType == "Loose") {
      if (pt >= looseEleIdSFhist_->GetYaxis()->GetXmax()) {
	int binX = looseEleIdSFhist_->GetXaxis()->FindBin(eta);
	int binY = looseEleIdSFhist_->GetYaxis()->GetLast();
	SF = looseEleIdSFhist_->GetBinContent(binX, binY);
      }
      else {
	int binX = looseEleIdSFhist_->GetXaxis()->FindBin(eta);
	int binY = looseEleIdSFhist_->GetYaxis()->FindBin(pt);
	SF = looseEleIdSFhist_->GetBinContent(binX, binY);
      }
    }
    else if (IdType == "Tight") {
      if (pt >= tightEleIdSFhist_->GetYaxis()->GetXmax()) {
	int binX = tightEleIdSFhist_->GetXaxis()->FindBin(eta);
	int binY = tightEleIdSFhist_->GetYaxis()->GetLast();
	SF = tightEleIdSFhist_->GetBinContent(binX, binY);
      }
      else {
	int binX = tightEleIdSFhist_->GetXaxis()->FindBin(eta);
	int binY = tightEleIdSFhist_->GetYaxis()->FindBin(pt);
	SF = tightEleIdSFhist_->GetBinContent(binX, binY);
      }
    }
  }
  return SF;
}

double ScaleFactorHandler::getIsoSF(const std::string& IsoType, float pt, float eta, const std::string& Flav) const {
  double SF = 0.0;
  if (Flav == "Muon") {
    if (IsoType == "Tight") {
      if (pt >= tightMuIsoSFhist_->GetXaxis()->GetXmax()) {
        int binX = tightMuIsoSFhist_->GetXaxis()->GetLast();
        int binY = tightMuIsoSFhist_->GetYaxis()->FindBin(std::abs(eta));
        SF = tightMuIsoSFhist_->GetBinContent(binX, binY);
      }
      else {
        int binX = tightMuIsoSFhist_->GetXaxis()->FindBin(pt);
        int binY = tightMuIsoSFhist_->GetYaxis()->FindBin(std::abs(eta));
        SF = tightMuIsoSFhist_->GetBinContent(binX, binY);
      }
    }
  }
  return SF;
}

double ScaleFactorHandler::getFF(float pt, float eta, const std::string& Flav) const {
  double FR = 0.0;
  if (Flav == "Muon") {
    if (pt >= muonFRhist_->GetXaxis()->GetXmax()) {
      int binX = muonFRhist_->GetXaxis()->GetLast();
      int binY = muonFRhist_->GetYaxis()->FindBin(std::fabs(eta));
      FR = muonFRhist_->GetBinContent(binX, binY);
    }
    else {
      int binX = muonFRhist_->GetXaxis()->FindBin(pt);
      int binY = muonFRhist_->GetYaxis()->FindBin(std::fabs(eta));
      FR = muonFRhist_->GetBinContent(binX, binY);
    }
  }
  else if (Flav == "Electron") {
    if (pt >= electronFRhist_->GetXaxis()->GetXmax()) {
      int binX = electronFRhist_->GetXaxis()->GetLast();
      int binY = electronFRhist_->GetYaxis()->FindBin(std::fabs(eta));
      FR = electronFRhist_->GetBinContent(binX, binY);
    }
    else {
      int binX = electronFRhist_->GetXaxis()->FindBin(pt);
      int binY = electronFRhist_->GetYaxis()->FindBin(std::fabs(eta));
      FR = electronFRhist_->GetBinContent(binX, binY);
    }
  }
  double FF = FR/(1-FR);
  //std::cout<<"Flavour : "<<Flav<<"\t"<<"pt : "<<pt<<"\t"<<"eta : "<<std::abs(eta)<<"\t"<<"FakeRate : "<<FR<<"\t"<<"FakeFactor : "<<FF<<"\n";
  return FF;
}
