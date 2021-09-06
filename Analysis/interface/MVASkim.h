#ifndef __MVASkim__h
#define __MVASkim__h

#include <memory>
#include <string>

class TTree;
class TFile;

typedef struct  
{
  float MCweight;
  float Channel;
  // ---------- Low level variables ---------- //
  // lep1 variables
  float px_lep1;
  float py_lep1;
  float pz_lep1;
  float E_lep1;
  float pt_lep1;
  float eta_lep1;
  float phi_lep1;
  int charge_lep1;
  int pdgid_lep1;
  // lep2 variables
  float px_lep2;
  float py_lep2;
  float pz_lep2;
  float E_lep2;
  float pt_lep2;
  float eta_lep2;
  float phi_lep2;
  int charge_lep2;
  int pdgid_lep2;
  // jets variables
  // jet1
  float px_jet1;
  float py_jet1;
  float pz_jet1;
  float E_jet1;
  float pt_jet1;
  float eta_jet1;
  float phi_jet1;
  float btag_jet1;
  // jet2
  float px_jet2;
  float py_jet2;
  float pz_jet2;
  float E_jet2;
  float pt_jet2;
  float eta_jet2;
  float phi_jet2;
  float btag_jet2;
  // jet3
  float px_jet3;
  float py_jet3;
  float pz_jet3;
  float E_jet3;
  float pt_jet3;
  float eta_jet3;
  float phi_jet3;
  float btag_jet3;
  // jet4
  float px_jet4;
  float py_jet4;
  float pz_jet4;
  float E_jet4;
  float pt_jet4;
  float eta_jet4;
  float phi_jet4;
  float btag_jet4;
  // MET
  float px_met;
  float py_met;
  float pz_met;
  float E_met;
  float pt_met;
  float eta_met;
  float phi_met;

} TreeVariablesResolved;

typedef struct  
{
  float MCweight;
  float Channel;
  float pt_lep1;    
  float pt_lep2;    
} TreeVariablesBoosted;

class MVASkim {
    
public:

  MVASkim(const std::string& filename);
  virtual ~MVASkim() {}

  void fill(const TreeVariablesResolved& varList);
  void fill(const TreeVariablesBoosted& varList);

  void close();

private:
  std::unique_ptr<TFile> _mvaFile; 
  std::unique_ptr<TTree> _treeR; 
  std::unique_ptr<TTree> _treeB; 

  TreeVariablesResolved _varListR;
  TreeVariablesBoosted  _varListB;
};
#endif
