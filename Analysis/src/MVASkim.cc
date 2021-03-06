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
  _tree = new TTree("RTree", "RTree");
  
  _tree->Branch("puwt", &_varList.puwt,  "puwt/F");
  _tree->Branch("evwt", &_varList.evwt,  "evwt/F");
  _tree->Branch("btagwt", &_varList.btagwt,  "btagwt/F");
  _tree->Branch("nLep", &_varList.nLep  ,  "nLep/I");     
  _tree->Branch("nJet", &_varList.nJet  ,  "nJet/I");     
  _tree->Branch("met", &_varList.met,  "met/F");       
  _tree->Branch("jet1pt", &_varList.jet1pt,  "jet1pt/F");   
  _tree->Branch("jet2pt", &_varList.jet2pt  ,  "jet2pt/F");   
  _tree->Branch("HT", &_varList.HT,  "HT/F");        
  _tree->Branch("HTvec", &_varList.HTvec  ,  "HTvec/F");     
  _tree->Branch("HTMET", &_varList.HTMET  ,  "HTMET/F");     
  _tree->Branch("HTMETvec", &_varList.HTMETvec  ,  "HTMETvec/F");     
  _tree->Branch("HTfrac", &_varList.HTfrac,  "HTfrac/F");        
  _tree->Branch("HTMETfrac", &_varList.HTMETfrac  ,  "HTMETfrac/F");     
  _tree->Branch("jInvM", &_varList.jInvM  ,  "jInvM/F");     
  _tree->Branch("j1j2DR", &_varList.j1j2DR  ,  "j1j2DR/F");     
  _tree->Branch("j1j2DEta", &_varList.j1j2DEta  ,  "j1j2DEta/F");  
  _tree->Branch("j1j2DPhi", &_varList.j1j2DPhi  ,  "j1j2DPhi/F");  
  _tree->Branch("j1j3DR", &_varList.j1j3DR  ,  "j1j3DR/F");     
  _tree->Branch("j1j3DEta", &_varList.j1j3DEta  ,  "j1j3DEta/F");  
  _tree->Branch("j1j3DPhi", &_varList.j1j3DPhi  ,  "j1j3DPhi/F");  
  _tree->Branch("j2j3DR", &_varList.j2j3DR  ,  "j2j3DR/F");     
  _tree->Branch("j2j3DEta", &_varList.j2j3DEta  ,  "j2j3DEta/F");  
  _tree->Branch("j2j3DPhi", &_varList.j2j3DPhi  ,  "j2j3DPhi/F");  

  _tree->Branch("minJJDR", &_varList.minJJDR  ,  "minJJDR/F");     
  _tree->Branch("minJJDPhi", &_varList.minJJDPhi  ,  "minJJDPhi/F");  
  _tree->Branch("maxJJDR", &_varList.maxJJDR  ,  "maxJJDR/F");     
  _tree->Branch("maxJJDPhi", &_varList.maxJJDPhi  ,  "maxJJDPhi/F");  


  _tree->Branch("lep1pt", &_varList.lep1pt  ,  "lep1pt/F");    
  _tree->Branch("lep2pt", &_varList.lep2pt  ,  "lep2pt/F");    
  _tree->Branch("LT", &_varList.LT  ,  "LT/F");        
  _tree->Branch("LTMET", &_varList.LTMET  ,  "LTMET/F");     
  _tree->Branch("LTMETvec", &_varList.LTMETvec  ,  "LTMETvec/F");
  _tree->Branch("LTMETfrac", &_varList.LTMETfrac  ,  "LTMETfrac/F");     
  _tree->Branch("LTMETfracInv", &_varList.LTMETfracInv  ,  "LTMETfracInv/F");       
  _tree->Branch("MT", &_varList.MT  ,  "MT/F");        
  _tree->Branch("MTovLT", &_varList.MTovLT  ,  "MTovLT/F");    
  _tree->Branch("l1l2InvM", &_varList.l1l2InvM  ,  "l1l2InvM/F");  
  _tree->Branch("l1l2DR", &_varList.l1l2DR  ,  "l1l2DR/F");    
  _tree->Branch("l1l2DEta", &_varList.l1l2DEta  ,  "l1l2DEta/F");  
  _tree->Branch("l1l2DPhi", &_varList.l1l2DPhi  ,  "l1l2DPhi/F");  

  _tree->Branch("l1MetDPhi", &_varList.l1MetDPhi  ,  "l1MetDPhi/F");
  _tree->Branch("l2MetDPhi", &_varList.l2MetDPhi  ,  "l2MetDPhi/F");
  _tree->Branch("j1MetDPhi", &_varList.j1MetDPhi  ,  "j1MetDPhi/F");
  _tree->Branch("j2MetDPhi", &_varList.j2MetDPhi  ,  "j2MetDPhi/F");
  _tree->Branch("j3MetDPhi", &_varList.j3MetDPhi  ,  "j3MetDPhi/F");

  _tree->Branch("j1l1DR", &_varList.j1l1DR  ,  "j1l1DR/F");    
  _tree->Branch("j1l1DEta", &_varList.j1l1DEta  ,  "j1l1DEta/F");  
  _tree->Branch("j1l1DPhi", &_varList.j1l1DPhi  ,  "j1l1DPhi/F");  

  _tree->Branch("j1l2DR", &_varList.j1l2DR  ,  "j1l2DR/F");    
  _tree->Branch("j1l2DEta", &_varList.j1l2DEta ,  "j1l2DEta/F");  
  _tree->Branch("j1l2DPhi", &_varList.j1l2DPhi  ,  "j1l2DPhi/F");  

  _tree->Branch("j2l1DR", &_varList.j2l1DR  ,  "j2l1DR/F");    
  _tree->Branch("j2l1DEta", &_varList.j2l1DEta  ,  "j2l1DEta/F");  
  _tree->Branch("j2l1DPhi", &_varList.j2l1DPhi  ,  "j2l1DPhi/F");  

  _tree->Branch("j2l2DR", &_varList.j2l2DR  ,  "j2l2DR/F");    
  _tree->Branch("j2l2DEta", &_varList.j2l2DEta  ,  "j2l2DEta/F");  
  _tree->Branch("j2l2DPhi", &_varList.j2l2DPhi  ,  "j2l2DPhi/F");  

  _tree->Branch("j3l1DR", &_varList.j3l1DR  ,  "j3l1DR/F");    
  _tree->Branch("j3l1DEta", &_varList.j3l1DEta  ,  "j3l1DEta/F");  
  _tree->Branch("j3l1DPhi", &_varList.j3l1DPhi  ,  "j3l1DPhi/F");  

  _tree->Branch("j3l2DR", &_varList.j3l2DR  ,  "j3l2DR/F");    
  _tree->Branch("j3l2DEta", &_varList.j3l2DEta  ,  "j3l2DEta/F");  
  _tree->Branch("j3l2DPhi", &_varList.j3l2DPhi  ,  "j3l2DPhi/F");  

  _tree->Branch("maxJLDR", &_varList.maxJLDR  ,  "maxJLDR/F");  
  _tree->Branch("minJLDR", &_varList.minJLDR  ,  "minJLDR/F");  
  _tree->Branch("maxJLDPhi", &_varList.maxJLDPhi  ,  "maxJLDPhi/F");  
  _tree->Branch("minJLDPhi", &_varList.minJLDPhi  ,  "minJLDPhi/F");  

  _tree->Branch("minLepMetDPhi", &_varList.minLepMetDPhi  ,  "minLepMetDPhi/F");
  _tree->Branch("maxLepMetDPhi", &_varList.maxLepMetDPhi  ,  "maxLepMetDPhi/F");
  _tree->Branch("minJetMetDPhi", &_varList.minJetMetDPhi  ,  "minJetMetDPhi/F");
  _tree->Branch("maxJetMetDPhi", &_varList.maxJetMetDPhi  ,  "maxJetMetDPhi/F");

  _tree->Branch("METovSQHT", &_varList.METovSQHT, "METovSQHT/F");
  _tree->Branch("METovSQLT", &_varList.METovSQLT, "METovSQLT/F");
  _tree->Branch("ST", &_varList.ST, "ST/F");
  _tree->Branch("LTovSqrtST", &_varList.LTovSqrtST, "LTovSqrtST/F");
  _tree->Branch("HTovSqrtST", &_varList.HTovSqrtST, "HTovSqrtST/F");
  _tree->Branch("STvec", &_varList.STvec,   "STvec/F");
  _tree->Branch("STMETvec", &_varList.STMETvec,  "STMETvec/F");


  _mvaFile->ls();
}
MVASkim::~MVASkim() {
  if (_tree) delete _tree;  
  if (_mvaFile) delete _mvaFile;
}
void MVASkim::fill(const TreeVariables& varList) {
  memcpy(&_varList, &varList, sizeof(varList));
  _mvaFile->cd();
  _tree->Fill();
}
void MVASkim::close() {
  //_mvaFile = TFile::Open(filename.c_str(), "RECREATE", "Skimmed Tree");
  _mvaFile->cd();
  //  _tree->Print();
  _tree->Write();
  _mvaFile->Save();
  _mvaFile->Write();
  //_mvaFile->Close();
}
