#ifndef __MVAnalysis__h
#define __MVAnalysis__h

#include <fstream>
#include <string>

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

typedef struct  
{
  //  float met;
  float HT;
  //float HTvec;
  //   float HTMETfrac;
  float j1j2DR;
  float j1j3DR;
  //float minJJDR;
  float maxJJDR;
  //float minJJDPhi;
  float LT;
  //float LTMETvec;
  //float LTMETfracInv;
  float l1l2InvM;
  float l1l2DR;
  //   float l1l2DPhi;
  float l1MetDPhi;
  //   float l2MetDPhi;
  float j1MetDPhi;
  float j2MetDPhi;
  //float j3MetDPhi;
  float j1l1DPhi;
  //   float j1l2DPhi;
  float minJLDPhi;
  //float maxJLDPhi;
  float maxLepMetDPhi;
  float minJetMetDPhi;

} InputVariables;

class MVAnalysis {
    
public:

  MVAnalysis(const std::string& mva_algo, const std::string& xmlfile);
  virtual ~MVAnalysis() {}

  double evaluate(const std::string& tag, const InputVariables& varList);

  InputVariables varList_;
  std::unique_ptr<TMVA::Reader> reader_;
};
#endif
