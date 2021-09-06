#include <iostream>
#include <memory>

#include "TFile.h"
#include "TTree.h"
#include "MVASkim.h"

using std::string;
using std::cout;
using std::endl;

MVASkim::MVASkim(const string& filename) {
  _mvaFile = std::make_unique<TFile>(filename.c_str(), "RECREATE", "Skimmed Tree");

  // Branches for Resolved tree
  _treeR = std::make_unique<TTree>("Resolved", "Resolved tree");
  _treeR->Branch("MCweight",            &_varListR.MCweight,          "MCweight/F");
  _treeR->Branch("Channel",             &_varListR.Channel,           "Channel/F");

  _treeR->Branch("px_lep1",             &_varListR.px_lep1,           "px_lep1/F");    
  _treeR->Branch("py_lep1",             &_varListR.py_lep1,           "py_lep1/F");    
  _treeR->Branch("pz_lep1",             &_varListR.pz_lep1,           "pz_lep1/F");    
  _treeR->Branch("E_lep1",              &_varListR.E_lep1,            "E_lep1/F");    
  _treeR->Branch("pt_lep1",             &_varListR.pt_lep1,           "pt_lep1/F");    
  _treeR->Branch("eta_lep1",            &_varListR.eta_lep1,          "eta_lep1/F");    
  _treeR->Branch("phi_lep1",            &_varListR.phi_lep1,          "phi_lep1/F");    
  _treeR->Branch("charge_lep1",         &_varListR.charge_lep1,       "charge_lep1/I");
  _treeR->Branch("pdgid_lep1",          &_varListR.pdgid_lep1,        "pdgid_lep1/I");

  _treeR->Branch("px_lep2",             &_varListR.px_lep2,           "px_lep2/F");    
  _treeR->Branch("py_lep2",             &_varListR.py_lep2,           "py_lep2/F");    
  _treeR->Branch("pz_lep2",             &_varListR.pz_lep2,           "pz_lep2/F");    
  _treeR->Branch("E_lep2",              &_varListR.E_lep2,            "E_lep2/F");    
  _treeR->Branch("pt_lep2",             &_varListR.pt_lep2,           "pt_lep2/F");    
  _treeR->Branch("eta_lep2",            &_varListR.eta_lep2,          "eta_lep2/F");    
  _treeR->Branch("phi_lep2",            &_varListR.phi_lep2,          "phi_lep2/F");    
  _treeR->Branch("charge_lep2",         &_varListR.charge_lep2,       "charge_lep2/I");
  _treeR->Branch("pdgid_lep2",          &_varListR.pdgid_lep2,        "pdgid_lep2/I");

  _treeR->Branch("px_jet1",             &_varListR.px_jet1,           "px_jet1/F");    
  _treeR->Branch("py_jet1",             &_varListR.py_jet1,           "py_jet1/F");    
  _treeR->Branch("pz_jet1",             &_varListR.pz_jet1,           "pz_jet1/F");    
  _treeR->Branch("E_jet1",              &_varListR.E_jet1,            "E_jet1/F");    
  _treeR->Branch("pt_jet1",             &_varListR.pt_jet1,           "pt_jet1/F");    
  _treeR->Branch("eta_jet1",            &_varListR.eta_jet1,          "eta_jet1/F");    
  _treeR->Branch("phi_jet1",            &_varListR.phi_jet1,          "phi_jet1/F");    
  _treeR->Branch("btag_jet1",           &_varListR.btag_jet1,         "btag_jet1/I");

  _treeR->Branch("px_jet2",             &_varListR.px_jet2,           "px_jet2/F");    
  _treeR->Branch("py_jet2",             &_varListR.py_jet2,           "py_jet2/F");    
  _treeR->Branch("pz_jet2",             &_varListR.pz_jet2,           "pz_jet2/F");    
  _treeR->Branch("E_jet2",              &_varListR.E_jet2,            "E_jet2/F");    
  _treeR->Branch("pt_jet2",             &_varListR.pt_jet2,           "pt_jet2/F");    
  _treeR->Branch("eta_jet2",            &_varListR.eta_jet2,          "eta_jet2/F");    
  _treeR->Branch("phi_jet2",            &_varListR.phi_jet2,          "phi_jet2/F");    
  _treeR->Branch("btag_jet2",           &_varListR.btag_jet2,         "btag_jet2/I");

  _treeR->Branch("px_jet3",             &_varListR.px_jet3,           "px_jet3/F");    
  _treeR->Branch("py_jet3",             &_varListR.py_jet3,           "py_jet3/F");    
  _treeR->Branch("pz_jet3",             &_varListR.pz_jet3,           "pz_jet3/F");    
  _treeR->Branch("E_jet3",              &_varListR.E_jet3,            "E_jet3/F");    
  _treeR->Branch("pt_jet3",             &_varListR.pt_jet3,           "pt_jet3/F");    
  _treeR->Branch("eta_jet3",            &_varListR.eta_jet3,          "eta_jet3/F");    
  _treeR->Branch("phi_jet3",            &_varListR.phi_jet3,          "phi_jet3/F");    
  _treeR->Branch("btag_jet3",           &_varListR.btag_jet3,         "btag_jet3/I");

  _treeR->Branch("px_jet4",             &_varListR.px_jet4,           "px_jet4/F");    
  _treeR->Branch("py_jet4",             &_varListR.py_jet4,           "py_jet4/F");    
  _treeR->Branch("pz_jet4",             &_varListR.pz_jet4,           "pz_jet4/F");    
  _treeR->Branch("E_jet4",              &_varListR.E_jet4,            "E_jet4/F");    
  _treeR->Branch("pt_jet4",             &_varListR.pt_jet4,           "pt_jet4/F");    
  _treeR->Branch("eta_jet4",            &_varListR.eta_jet4,          "eta_jet4/F");    
  _treeR->Branch("phi_jet4",            &_varListR.phi_jet4,          "phi_jet4/F");    
  _treeR->Branch("btag_jet4",           &_varListR.btag_jet4,         "btag_jet4/I");

  _treeR->Branch("px_met",             &_varListR.px_met,           "px_met/F");    
  _treeR->Branch("py_met",             &_varListR.py_met,           "py_met/F");    
  _treeR->Branch("pz_met",             &_varListR.pz_met,           "pz_met/F");    
  _treeR->Branch("E_met",              &_varListR.E_met,            "E_met/F");    
  _treeR->Branch("pt_met",             &_varListR.pt_met,           "pt_met/F");    
  _treeR->Branch("eta_met",            &_varListR.eta_met,          "eta_met/F");    
  _treeR->Branch("phi_met",            &_varListR.phi_met,          "phi_met/F");    


  // Branches for Boosted tree
  _treeB = std::make_unique<TTree>("Boosted", "Boosted tree");
  _treeB->Branch("MCweight", &_varListB.MCweight, "MCweight/F");
  _treeB->Branch("Channel",  &_varListB.Channel,  "Channel/F");
  _treeB->Branch("pt_lep1",  &_varListB.pt_lep1,  "pt_lep1/F");    
  _treeB->Branch("pt_lep2",  &_varListB.pt_lep2,  "pt_lep2/F");    

  _mvaFile->ls();
}
void MVASkim::fill(const TreeVariablesResolved& varList) {
  memcpy(&_varListR, &varList, sizeof(varList));
  _treeR->Fill();
}
void MVASkim::fill(const TreeVariablesBoosted& varList) {
  memcpy(&_varListB, &varList, sizeof(varList));
  _treeB->Fill();
}
void MVASkim::close() {
  _mvaFile->cd();
  _treeR->Write();
  _treeB->Write();
}
