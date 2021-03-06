#ifndef __ZCandidate__hh
#define __ZCandidate__hh

#include "TLorentzVector.h"

namespace ZSpace {
  constexpr double MZnominal = 91.1876;
  enum class ZType: int {
    unkwn = -1, mumu = 0, ee
  };
  enum class llType: int {
    unkwn = -1, mue = 0, emu
  };
  enum class lTauType: int {
    unkwn = -1, mutau = 0, etau
  };
  struct ZCandidate {
    int l1Index;
    TLorentzVector l1P4;
    int l1Charge;
    double l1Isolation;
  
    int l2Index;
    TLorentzVector l2P4;
    int l2Charge;
    double l2Isolation;
  
    int flavour; // -1: unknown, 0: mumu, 1: ee etc.
    TLorentzVector p4;
    double mass;
    double massDiff;
    double dEtall;
    double dPhill;
    double dRll;
  };
}
#endif
