#ifndef __MVASkim__h
#define __MVASkim__h

#include <fstream>
#include <string>

class TTree;
class TFile;

typedef struct  
{
  float puwt;
  float evwt;
  float btagwt;
  int nLep;     
  int nJet;     
  float met;       
  float jet1pt;   
  float jet2pt;   
  float HT;        
  float HTvec;     
  float HTMET;        
  float HTMETvec;     
  float HTfrac;        
  float HTMETfrac;     
  float jInvM;     
  float j1j2DR;    
  float j1j2DEta;  
  float j1j2DPhi;  
  float j1j3DR;    
  float j1j3DEta;  
  float j1j3DPhi;  
  float j2j3DR;    
  float j2j3DEta;  
  float j2j3DPhi;  

  float minJJDR;    
  float minJJDPhi;  
  float maxJJDR;    
  float maxJJDPhi;  


  float lep1pt;    
  float lep2pt;    
  float LT;        
  float LTMET;    
  float LTMETvec;  
  float LTMETfrac;  
  float LTMETfracInv;  
  float MT;        
  float MTovLT;    
  float l1l2InvM;  
  float l1l2DR;    
  float l1l2DEta;  
  float l1l2DPhi;  

  float l1MetDPhi;
  float l2MetDPhi;
  float j1MetDPhi;
  float j2MetDPhi;
  float j3MetDPhi;

  float j1l1DR;    
  float j1l1DEta;  
  float j1l1DPhi;  

  float j1l2DR;    
  float j1l2DEta;  
  float j1l2DPhi;  

  float j2l1DR;    
  float j2l1DEta;  
  float j2l1DPhi;  

  float j2l2DR;    
  float j2l2DEta;  
  float j2l2DPhi;  


  float j3l1DR;    
  float j3l1DEta;  
  float j3l1DPhi;  

  float j3l2DR;    
  float j3l2DEta;  
  float j3l2DPhi;  

  float maxJLDR;
  float minJLDR;
  float maxJLDPhi;
  float minJLDPhi;

  float minLepMetDPhi;
  float maxLepMetDPhi;
  float minJetMetDPhi;
  float maxJetMetDPhi;

  float METovSQHT;
  float METovSQLT;
  float ST;
  float LTovSqrtST;
  float HTovSqrtST;
  float STvec;
  float STMETvec;

} TreeVariables;

class MVASkim {
    
public:

  MVASkim(const std::string& filename);
  virtual ~MVASkim();

  void fill(const TreeVariables& varList);
  void close();

  TFile* _mvaFile;
  TTree* _tree;

  TreeVariables _varList;
};
#endif
