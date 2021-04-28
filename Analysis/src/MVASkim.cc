#include <iostream>
#include <memory>
#include "TFile.h"
#include "TTree.h"
#include "MVASkim.h"

using std::string;
using std::cout;
using std::endl;

MVASkim::MVASkim(const string& filename) {
  _mvaFile = TFile::Open(filename.c_str(), "RECREATE", "Skimmed Tree");
  _mvaFile->cd();

  // Branches for Resolved tree
  _treeR = new TTree("Resolved", "Resolved tree");
  _treeR -> Branch("MCweight",           &_varListR.MCweight,       "MCweight/F");
  _treeR -> Branch("Channel",            &_varListR.Channel,         "Channel/F");
  _treeR -> Branch("pt_lep1",            &_varListR.pt_lep1,         "pt_lep1/F");    
  _treeR -> Branch("pt_lep2",            &_varListR.pt_lep2,         "pt_lep2/F");    
  _treeR -> Branch("invM_jets",          &_varListR.invM_jets,       "invM_jets/F");    

  // Branches for Boosted tree
  _treeB = new TTree("Boosted",  "Boosted tree");
  _treeB -> Branch("MCweight",           &_varListB.MCweight,       "MCweight/F");
  _treeB -> Branch("Channel",            &_varListB.Channel,         "Channel/F");
  _treeB -> Branch("pt_lep1",            &_varListB.pt_lep1,         "pt_lep1/F");    
  _treeB -> Branch("pt_lep2",            &_varListB.pt_lep2,         "pt_lep2/F");    

  _mvaFile->ls();
}

MVASkim::~MVASkim() {
  if (_treeR)   delete _treeR;  
  if (_treeB)   delete _treeB;  
  if (_mvaFile) delete _mvaFile;
}
/*
void MVASkim::fill(const TreeVariablesResolved& varListR, const TreeVariablesBoosted& varListB,
		   bool isResolved = false, bool isBoosted = false) {
  if (isResolved)  memcpy(&_varListR, &varListR, sizeof(varListR));
  else if (isBoosted) memcpy(&_varListB, &varListB, sizeof(varListB));
  _mvaFile->cd();
  if (isResolved) _treeR->Fill();
  else if (isBoosted) _treeB->Fill();
}
*/

void MVASkim::close() {
  _mvaFile->cd();
  _treeR->Write();
  _treeB->Write();
  //_mvaFile->Save();
  //_mvaFile->Write();
  //_mvaFile->Close();
}
