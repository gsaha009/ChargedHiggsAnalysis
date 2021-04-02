#ifndef __PhysicsObjects_h
#define __PhysicsObjects_h

#include <vector>
#include <map>
#include <string>
#include "TLorentzVector.h"
#include "TObject.h"

namespace vhtm {
  class Event: public TObject {
  public:
    Event();
    ~Event() {}

    //info
    unsigned int run;
    unsigned long long event;
    unsigned int lumis;
    //Weight
    float genEvWt;
    float PUWeight;
    float btagWeight_CSVV2;
    float btagWeight_CMVA;
    //Primary Vetrex
    int nPV;
    int nGoodPV;
    float PVscore;
    float PVchi2;
    float PVndf;
    float PVx;
    float PVy;
    float PVz;
    //LHE Info
    int nLHEParticles;
    int nLHEJets;
    //Gen Info
    int nGenParticles;
    //UserDefined
    char HLT_SingleMu;
    char HLT_SingleEle;
    char HLT_DoubleMu;
    char HLT_DoubleEle;
    int evType;

    ClassDef(Event, 1)
  };

  class GenParticle: public TObject {
  public:
    GenParticle();
    ~GenParticle() {}
    
    unsigned int index;
    float pt;
    float eta;
    float phi;
    float mass;
    int motherIdx;
    int pdgId;
    int status;
    int statusFlags;
  
    ClassDef(GenParticle, 1)
  };

  class LHEParticle: public TObject {
  public:
    LHEParticle();
    ~LHEParticle() {}

    unsigned int index;
    float pt;
    float eta;
    float phi;
    float mass;
    int pdgId;

    ClassDef(LHEParticle, 1)
  };

  class Muon: public TObject {
  public:
    Muon();
    ~Muon() {}

    unsigned int index;
    float pt;
    float eta;
    float SCeta;
    float phi;
    float mass;
    int charge;
    int jetIdx;
    int genIdx;
    unsigned char genFlv;
    bool mediumId;
    bool tightId;
    unsigned char highPtId;
    float pfRelIso03_all;
    float pfRelIso04_all;

    ClassDef(Muon, 1)
  };


  class Electron: public TObject {
  public:
    Electron();
    ~Electron() {}

    unsigned int index;
    float pt;
    float eta;
    float SCdEta;
    float SCeta;
    float phi;
    float mass;
    int charge;
    int jetIdx;
    int phoIdx;
    int genIdx;
    unsigned char genFlv;

    ClassDef(Electron, 1)
  };

  class MET: public TObject {
  public:
    MET();
    ~MET() {}
  
    float pt;
    float phi;
    float signf;
    float sumEt;
  
    ClassDef(MET, 1)
  };
  class Tau: public TObject {
  public:
    Tau();
    ~Tau() {}
  
    unsigned int index;
    float pt;
    float eta;
    float phi;
    float mass;
    int charge;
    int jetIdx;
  
    ClassDef(Tau, 1)
  };

  class Jet: public TObject {
  public:
    Jet();
    ~Jet() {}
    
    unsigned int index;
    float pt;
    float eta;
    float phi;
    float mass;
    int muIdx1;
    int muIdx2;
    int elIdx1;
    int elIdx2;
    float btagDeepFlavB;
    float btagCSVV2;
  
    ClassDef(Jet, 1)
  };

  class FatJet: public TObject {
  public:
    FatJet();
    ~FatJet() {}
    
    unsigned int index;
    float pt;
    float eta;
    float phi;
    float mass;
    float softDropMass;
    float n2b1;
    float n3b1;
    unsigned char nBHadrons;
    unsigned char nCHadrons;
    int hadronFlavour;
    int subJetIdx1;
    int subJetIdx2;
    float btagDeepB;
    float btagCSVV2;
    float rawFactor;
    float tau1;  
    float tau2;  
    float tau3;  
    float tau4;  
    float deepTag_WvsQCD;
    float deepTag_ZvsQCD;
    float deepTag_TvsQCD;
    float deepTagMD_WvsQCD;
    float deepTagMD_ZvsQCD;
    float deepTagMD_TvsQCD;
    int electronIdx3SJ;
    int muonIdx3SJ;

    ClassDef(FatJet, 1)
  };

  class SubJet: public TObject {
  public:
    SubJet();
    ~SubJet() {}
    
    unsigned int index;
    float pt;
    float eta;
    float phi;
    float mass;
    float btagDeepB;
    float rawFactor;

    ClassDef(SubJet, 1)
  };

}
#endif
