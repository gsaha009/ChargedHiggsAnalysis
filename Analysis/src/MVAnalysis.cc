#include <iostream>
#include <memory>
#include "MVAnalysis.h"

using std::string;
using std::cout;
using std::endl;

MVAnalysis::MVAnalysis(const string& mva_algo, const string& xmlfile) 
{
  reader_ = std::make_unique<TMVA::Reader>("!Silent");

  //  reader_->AddVariable("met", &varList_.met);
  reader_->AddVariable("HT", &varList_.HT);
  //  reader_->AddVariable("HTvec", &varList_.HTvec);
  //   reader_->AddVariable("HTMETfrac", &varList_.HTMETfrac);
  reader_->AddVariable("j1j2DR", &varList_.j1j2DR);
  reader_->AddVariable("j1j3DR", &varList_.j1j3DR);
  //reader_->AddVariable("minJJDR", &varList_.minJJDR);
  reader_->AddVariable("maxJJDR", &varList_.maxJJDR);
  // reader_->AddVariable("minJJDPhi", &varList_.minJJDPhi);
  reader_->AddVariable("LT", &varList_.LT);
  //reader_->AddVariable("LTMETvec", &varList_.LTMETvec);
  //reader_->AddVariable("LTMETfracInv", &varList_.LTMETfracInv);
  reader_->AddVariable("l1l2InvM", &varList_.l1l2InvM);
  reader_->AddVariable("l1l2DR", &varList_.l1l2DR);
  //   reader_->AddVariable("l1l2DPhi", &varList_.l1l2DPhi);
  reader_->AddVariable("l1MetDPhi", &varList_.l1MetDPhi);
  //   reader_->AddVariable("l2MetDPhi", &varList_.l2MetDPhi);
  reader_->AddVariable("j1MetDPhi", &varList_.j1MetDPhi);
  reader_->AddVariable("j2MetDPhi", &varList_.j2MetDPhi);
  //reader_->AddVariable("j3MetDPhi", &varList_.j3MetDPhi);
  reader_->AddVariable("j1l1DPhi", &varList_.j1l1DPhi);
  //   reader_->AddVariable("j1l2DPhi", &varList_.j1l2DPhi);
  reader_->AddVariable("minJLDPhi", &varList_.minJLDPhi);
  //reader_->AddVariable("maxJLDPhi", &varList_.maxJLDPhi);
  reader_->AddVariable("maxLepMetDPhi", &varList_.maxLepMetDPhi);
  reader_->AddVariable("minJetMetDPhi", &varList_.minJetMetDPhi);

  reader_->BookMVA(mva_algo.c_str(), xmlfile);
}
double MVAnalysis::evaluate(const string& mva_algo, const InputVariables& varList) {
  memcpy(&varList_, &varList, sizeof(varList)); // use move syntax here
  return reader_->EvaluateMVA(mva_algo.c_str());
}
