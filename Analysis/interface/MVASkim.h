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
