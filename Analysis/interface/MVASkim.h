#ifndef __MVASkim__h
#define __MVASkim__h

#include <fstream>
#include <string>

class TTree;
class TFile;

typedef struct  
{
  float MCweight;
  float Channel;
  float pt_lep1;    
  float pt_lep2;
  float invM_jets;
    
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
  virtual ~MVASkim();

  //void fill(const TreeVariablesResolved& varListR, const TreeVariablesBoosted& varListB,
  //	    bool isResolved = false, bool isBoosted = false);

  TFile* _mvaFile;
  TTree* _treeR;
  TTree* _treeB;

  template <typename T>
    void fill(const T& varList, bool isResolved = false, bool isBoosted = false) {
    if (isResolved) memcpy(&_varListR, &varList, sizeof(varList));
    else if (isBoosted) memcpy(&_varListB, &varList, sizeof(varList));
    //_mvaFile->cd();
    if (isResolved) _treeR->Fill();
    else if (isBoosted) _treeB->Fill();
  }
  void close();

  TreeVariablesResolved _varListR;
  TreeVariablesBoosted  _varListB;
};
#endif
