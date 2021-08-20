#include "interface/PhysicsObjects.h"

#define NEL(x) (sizeof((x))/sizeof((x)[0]))

ClassImp(vhtm::Event)
ClassImp(vhtm::GenParticle)
ClassImp(vhtm::LHEParticle)
ClassImp(vhtm::Muon)
ClassImp(vhtm::Electron)
ClassImp(vhtm::MET)
ClassImp(vhtm::Tau)
ClassImp(vhtm::Jet)
ClassImp(vhtm::FatJet)
ClassImp(vhtm::SubJet)

vhtm::Event::Event():
  run(0),
  event(0),
  lumis(0),
  genEvWt(-1),
  PUWeight(-999),
  btagWeight_CSVV2(-1),
  btagWeight_CMVA(-1),
  nPV(-1),
  nGoodPV(-1),
  PVscore(-999),
  PVchi2(-999),
  PVndf(-999),
  PVx(-999),
  PVy(-999),
  PVz(-999),
  nLHEParticles(-999),
  nLHEJets(-999),
  nGenParticles(-999)
{}

vhtm::GenParticle::GenParticle():
  index(-1),
  pt(-999),
  eta(-999),
  phi(-999),
  mass(-999),
  motherIdx(-1),
  pdgId(-1),
  status(-1),
  statusFlags(-1)
{}

vhtm::LHEParticle::LHEParticle():
  index(-1),
  pt(-999),
  eta(-999),
  phi(-999),
  mass(-999),
  pdgId(-1)
{}

vhtm::Muon::Muon():
  index(-1),
  pt(-999),
  eta(-999),
  SCeta(-999),
  phi(-999),
  mass(-999),
  charge(-999),
  jetIdx(-111),
  genIdx(-111),
  genFlv(-111),
  mediumId(false),
  tightId(false),
  highPtId(-111),
  pfRelIso03_all(-999),
  pfRelIso04_all(-999)
{}

vhtm::Electron::Electron():
  index(-1),
  pt(-999),
  eta(-999),
  SCdEta(-999),
  SCeta(-999),
  phi(-999),
  mass(-999),
  charge(-999),
  jetIdx(-111),
  genIdx(-111),
  genFlv(-111)
{}


vhtm::MET::MET():
  pt(-999),
  phi(-999),
  signf(-999),
  sumEt(-999)
{}

vhtm::Tau::Tau():
  index(-1),
  pt(-999),
  eta(-999),
  phi(-999),
  mass(-999),
  charge(-999),
  jetIdx(-111)
{}

vhtm::Jet::Jet():
  index(-1),
  pt(-999),
  eta(-999),
  phi(-999),
  mass(-999),
  muIdx1(-1),
  muIdx2(-1),
  elIdx1(-1),
  elIdx2(-1),
  btagDeepFlavB(-999),
  btagCSVV2(-999)
{}

vhtm::FatJet::FatJet():
  index(-1),
  pt(-999),
  eta(-999),
  phi(-999),
  mass(-999),
  softDropMass(-999),
  n2b1(-999),
  n3b1(-999),
  nBHadrons(-1),
  nCHadrons(-1),
  subJetIdx1(-999),
  subJetIdx2(-999),
  btagCSVV2(-999),
  tau1(-1),
  tau2(-1),
  tau3(-1),
  tau4(-1),
  deepTag_WvsQCD(-999),
  deepTag_ZvsQCD(-999),
  deepTag_TvsQCD(-999),
  deepTagMD_WvsQCD(-999),
  deepTagMD_ZvsQCD(-999),
  deepTagMD_TvsQCD(-999),
  electronIdx3SJ(-111),
  muonIdx3SJ(-111)
{}

vhtm::SubJet::SubJet():
  index(-1),
  pt(-999),
  eta(-999),
  phi(-999),
  mass(-999),
  btagDeepB(-999),
  rawFactor(-999)
{}
