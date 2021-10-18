#ifndef __HISTBOOKER__H
#define __HISTBOOKER__H

#include <string>
#include "AnaUtil.h"

class HistBooker {
public: 
  HistBooker() {}
  virtual ~HistBooker() {}
  void bookHist1D(const char* hname,
                  const char* htitle,
                  int nbins, float xlow, float xhigh);
  void bookHist1D(const char* hname,
                  const char* htitle,
                  int nbins, float xlow, float xhigh,
                  const std::vector<std::string> &Flags1,
                  const std::vector<std::string> &Flags2);

  void bookHistograms(bool isMC);
  void bookHistograms_DY(bool isMC);
};
#endif
