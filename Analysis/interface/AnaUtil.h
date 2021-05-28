#ifndef __ANAUTIL__HH
#define __ANAUTIL__HH

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>

#include "TMath.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TH1K.h"

using std::ostream;

namespace AnaUtil {
  // Templated functioned must be defined in the header itself
  void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters=" ");
  void bit_print(int value, int pos=32, ostream& os=std::cout);
  template <typename T>
  T deltaPhiT(T phi1, T phi2) {
    T result = phi1 - phi2;
    while (result > TMath::Pi()) result -= 2*TMath::Pi();
    while (result <= -TMath::Pi()) result += 2*TMath::Pi();
    return result;
  }  
  // get p4
  template <class T>
    TLorentzVector getP4(const T& obj) {
    TLorentzVector lv;
    lv.SetPtEtaPhiM(obj.pt, obj.eta, obj.phi, obj.mass);
    return lv;
  }
  //  
  template <typename T1, typename T2>
    bool within(T1 pt_, T2 pt1, T2 pt2) {
    if (pt_ > pt1 && pt_ <= pt2) return true;
    else return false;
  }
  // Transverse mass
  template<typename LVector1, typename LVector2>
    double Calculate_MT(const LVector1& lepton_p4, const LVector2& met_p4){
    const double delta_phi = TVector2::Phi_mpi_pi(lepton_p4.Phi() - met_p4.Phi());
    return std::sqrt( 2.0 * lepton_p4.Pt() * met_p4.Pt() * ( 1.0 - std::cos(delta_phi) ) );
  }
  template<typename LVector1, typename LVector2>
    double Calculate_MTfix(const LVector1& lepton_p4, const LVector2& met_p4) {
    const double delta_phi = TVector2::Phi_mpi_pi(lepton_p4.Phi() - met_p4.Phi());
    return std::sqrt( 2.0 * 35.0 * met_p4.Pt() * ( 1.0 - std::cos(delta_phi)));
  }
  template<typename LVector1, typename LVector2, typename LVector3>
    double Calculate_TotalMT(const LVector1& lepton1_p4, const LVector2&
			     lepton2_p4, const LVector3& met_p4)
  {
    const double mt_1 = Calculate_MT(lepton1_p4, met_p4);
    const double mt_2 = Calculate_MT(lepton2_p4, met_p4);
    const double mt_ll = Calculate_MT(lepton1_p4, lepton2_p4);
    return std::sqrt(std::pow(mt_1, 2) + std::pow(mt_2, 2) +
		     std::pow(mt_ll, 2));
  }
  // if an event has paased HLT condition or not
  template<typename T1, typename T2>
    bool isTriggered(const std::vector<T1>& paths, const std::vector<T2>& scores){
    for (size_t i = 0; i < paths.size(); ++i)
      if (scores[i]) return true;
    return false;
  }

  double deltaPhi(double phia, double phib);
  double deltaPhi(const TLorentzVector& a, const TLorentzVector& b);
  double deltaR(const TLorentzVector& a, const TLorentzVector& b);
  bool sameObject(const TLorentzVector& lv1, const TLorentzVector& lv2);
  double cutValue(const std::map<std::string, double>& m, const std::string& cname);
  const std::map<std::string, double>& cutMap(const std::map<std::string, std::map<std::string, double>>& hmap, const std::string& mkey);
  void buildList(const std::vector<std::string>& tokens, std::vector<std::string>& list);
  void buildMap(const std::vector<std::string>& tokens, std::map<std::string, int>& hmap);
  void buildMap(const std::vector<std::string>& tokens, std::unordered_map<std::string, int>& hmap);
  void storeCuts(const std::vector<std::string>& tokens, std::map<std::string, std::map<std::string, double>>& hmap);
  void showCuts(const std::map<std::string, std::map<std::string, double> >& hmap, ostream& os=std::cout);
  // ------------------------------------------------------------------------
  // Convenience routine for filling 1D histograms. We rely on root to keep 
  // track of all the histograms that are booked all over so that we do not 
  // have to use any global variables to save the histogram pointers. Instead, 
  // we use the name of the histograms and gROOT to retrieve them from the 
  // Root object pool whenever necessary. This is the closest one can go to 
  // hbook and ID based histogramming
  // -------------------------------------------------------------------------
  // dynamic approach for 1d histogram booking and filling
  TH1* getHist1D(const char* hname, int nbins, float xlow, float xhigh, const char* channel);
  TH1* getHist1D(const std::string& hname, int nbins, float xlow, float xhigh, const std::string& channel);
  template <class T>
    bool fillHist1D(const char* hname, T value, int nbins, float xlow, float xhigh,
		    const std::map<std::string, bool> &channelFlags, double w=1.0, bool selection=true) {
    if (!selection) return false;
    for (auto const& x : channelFlags) {
      std::string channel = x.first;
      bool check = x.second;
      if (check) {
	TH1* h = getHist1D(hname, nbins, xlow, xhigh, channel.c_str());
	if (h == nullptr) continue;
	h->Fill(value, w);
      }
    }
    return true;
  }
  template <class T>
    bool fillHist1D(const std::string& hname, T value, int nbins, float xlow, float xhigh,
		    const std::map<std::string, bool> channelFlags, double w=1.0, bool selection=true) {
    return fillHist1D(hname.c_str(), value, nbins, xlow, xhigh, channelFlags, w, selection);
  }

  // static approach for 1d histogram booking and filling
  TH1* getHist1DBasic(const char* hname);
  TH1* getHist1DBasic(const std::string& hname);
  template <class T>
    bool fillHist1DBasic(const char* hname, T value, double w=1.0, bool flag=true) {
    if (flag) {
      TH1* h = getHist1DBasic(hname);
      if (h == nullptr) return false;
      h->Fill(value, w);
    }
    return true;
  }
  template <class T>
    bool fillHist1DBasic(const std::string& hname, T value, double w=1.0, bool flag=true) {
    return fillHist1DBasic(hname.c_str(), value, w, flag);
  }
  // ---------------------------------------------
  // Convenience routine for filling 2D histograms
  // ---------------------------------------------
  // dynamic approach for 2d histogram booking and filling
  TH2* getHist2D(const char* hname, int nbinsX, float xlow, float xhigh, int nbinsY, float ylow, float yhigh, /*const char* region,*/ const char* channel);
  TH2* getHist2D(const std::string& hname, int nbinsX, float xlow, float xhigh, int nbinsY, float ylow, float yhigh, /*const std::string& region,*/ const std::string& channel);
  template <class T1, class T2>
    bool fillHist2D(const char* hname, T1 xvalue, T2 yvalue, int nbinsX, float xlow, float xhigh, int nbinsY, float ylow, float yhigh,
		    /*const char* region,*/ const std::map<std::string, bool> &channelFlags, double w=1.0, bool selection=true) {
    if (!selection) return false;
    for (auto const& x : channelFlags) {
      std::string channel = x.first;
      bool check = x.second;
      if (check) {
	TH2* h = getHist2D(hname, nbinsX, xlow, xhigh, nbinsY, ylow, yhigh, /*region, */channel.c_str());
	if (h == nullptr) continue;
	h->Fill(xvalue, yvalue, w);
      }
    }
    return true;
  }
  template <class T1, class T2>
    bool fillHist2D(const std::string& hname, T1 xvalue, T2 yvalue, int nbinsX, float xlow, float xhigh, int nbinsY, float ylow, float yhigh, 
		    /*const std::string& region,*/ const std::map<std::string, bool> channelFlags, double w=1.0, bool selection=true) {
    return fillHist2D(hname.c_str(), xvalue, yvalue, nbinsX, xlow, xhigh, nbinsY, ylow, yhigh, /*region.c_str(), */channelFlags, w, selection);
  }

  // static approach for 2d histogram booking and filling
  TH2* getHist2DBasic(const char* hname);
  TH2* getHist2DBasic(const std::string& hname);
  template <class T1, class T2>
    bool fillHist2DBasic(const char* hname, T1 xvalue, T2 yvalue, double w=1.0, bool flag=true) {
    if (flag) {
      TH2* h = getHist2DBasic(hname);
      if (h == nullptr) return false;
      h->Fill(xvalue, yvalue, w);
    }
    return true;
  }
  template <class T1, class T2>
    bool fillHist2DBasic(const std::string& hname, T1 xvalue, T2 yvalue, double w=1.0, bool flag=true) {
    return fillHist2DBasic(hname.c_str(), xvalue, yvalue, w, flag);
  }
  // ---------------------------------------------
  // Convenience routine for filling 3D histograms // TODO :: dynamic approach
  // ---------------------------------------------
  TH3* getHist3D(const char* hname);
  TH3* getHist3D(const std::string& hname);
  template <class T1, class T2, class T3>
  bool fillHist3D(const char* hname, T1 xvalue, T2 yvalue, T3 zvalue, double w=1.0) {
    TH3* h = getHist3D(hname);
    if (h == nullptr) return false;
    h->Fill(xvalue, yvalue, zvalue, w);
    return true;
  }
  template <class T1, class T2, class T3>
  bool fillHist3D(const std::string& hname, T1 xvalue, T2 yvalue, T3 zvalue, double w=1.0) {
    return fillHist3D(hname.c_str(), xvalue, yvalue, zvalue, w);
  }

  // --------------------------------------------------
  // Convenience routine for filling profile histograms
  // --------------------------------------------------
  TProfile* getProfile(const char* hname);
  TProfile* getProfile(const std::string& hname);
  bool fillProfile(const char *hname, float xvalue, float yvalue, double w=1.0);
  bool fillProfile(const std::string& hname, float xvalue, float yvalue, double w=1.0);

  /* print_list_elements()
   * - prints optional C-string optcstr followed by
   * - all elements of the collection coll
   */
  template <class T>
  void showList(const T& coll, const char* optcstr="", std::ostream& os=std::cout) {
    os << optcstr << ", Total # = " << coll.size() << ":" << std::endl;
    for (auto const& v: coll)
      os << v << std::endl;
  }  
  template <class T1, class T2>
  void showMap(const std::map<T1,T2>& m, const char* optcstr="", std::ostream& os=std::cout) {
    os << optcstr << std::endl;
    for (auto const& k: m)
      os << k.first << std::endl;
  }
  template <class T1, class T2>
  void showMap(const std::unordered_map<T1,T2>& m, const char* optcstr="", std::ostream& os=std::cout) {
    os << optcstr << std::endl;
    for (auto const& k: m)
      os << k.first << std::endl;
  }
  /* copyList */
  template <class T>
  void copyList (const T& sourceColl, T& destColl) {
    destColl.clear();
    for (auto const& v: sourceColl)   
      destColl.push_back(v); 
  }
  void scaleHistogram(const std::string& hname, double fac);
  void showEfficiency(const std::string& hname,
		      const std::vector<std::string>& slist,
		      const std::string& header,
		      const std::string& tag="Events", std::ostream& os=std::cout);
  void showYield(const std::string& hname,
		 const std::vector<std::string>& slist,
		 const std::string& header,
		 const std::string& tag="Events", std::ostream& os=std::cout);
}
#endif
