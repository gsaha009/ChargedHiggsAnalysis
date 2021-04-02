#ifndef __LeptonCand__hh
#define __LeptonCand__hh

#include "TLorentzVector.h"

struct LeptonCand {

  unsigned int index {100};
  float pt {-999};
  float eta {-999};
  float SCeta {-999};
  float phi {-999};
  float mass {-999};
  int charge {-999};
  int flavour {-111}; // 1:muon, 2:electron, -1:others

};
#endif
