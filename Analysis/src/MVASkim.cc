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
  _treeR->Branch("MCweight",  &_varListR.MCweight,  "MCweight/F");
  _treeR->Branch("Channel",   &_varListR.Channel,   "Channel/F");
  _treeR->Branch("pt_lep1",   &_varListR.pt_lep1,   "pt_lep1/F");    
  _treeR->Branch("pt_lep2",   &_varListR.pt_lep2,   "pt_lep2/F");    
  _treeR->Branch("invM_jets", &_varListR.invM_jets, "invM_jets/F");    

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
  _mvaFile->Write();
  _mvaFile->Close();
}
