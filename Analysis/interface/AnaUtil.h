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
  void bit_print(int value, int pos=32, std::ostream& os=std::cout);
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
  template <typename T>
  bool within(T pt_, T pt1, T pt2) {
    if (pt_ > pt1 && pt_ <= pt2) return true;
    else return false;
  }
  // see AN-13-178
  // Transverse mass
  template<typename T>
  double Calculate_MT(const T& lepton_p4, const T& met_p4){
    double delta_phi = TVector2::Phi_mpi_pi(lepton_p4.Phi() - met_p4.Phi());
    return std::sqrt( 2.0 * lepton_p4.Pt() * met_p4.Pt() * ( 1.0 - std::cos(delta_phi) ) );
  }
  template<typename T>
  double Calculate_MTfix(const T& lepton_p4, const T& met_p4) {
    double delta_phi = TVector2::Phi_mpi_pi(lepton_p4.Phi() - met_p4.Phi());
    return std::sqrt( 2.0 * 35.0 * met_p4.Pt() * ( 1.0 - std::cos(delta_phi)));
  }
  template<typename T>
  double Calculate_TotalMT(const T& lepton1_p4, const T& lepton2_p4, const T& met_p4)
  {
    double mt_1 = Calculate_MT(lepton1_p4, met_p4);
    double mt_2 = Calculate_MT(lepton2_p4, met_p4);
    double mt_ll = Calculate_MT(lepton1_p4, lepton2_p4);
    return std::sqrt(std::pow(mt_1, 2) + std::pow(mt_2, 2) + std::pow(mt_ll, 2));
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
  void showCuts(const std::map<std::string, std::map<std::string, double> >& hmap, std::ostream& os=std::cout);
  // ------------------------------------------------------------------------
  // Convenience routine for filling 1D histograms. We rely on root to keep 
  // track of all the histograms that are booked all over so that we do not 
  // have to use any global variables to save the histogram pointers. Instead, 
  // we use the name of the histograms and gROOT to retrieve them from the 
  // Root object pool whenever necessary. This is the closest one can go to 
  // hbook and ID based histogramming
  // -------------------------------------------------------------------------
  // dynamic approach for 1d histogram booking and filling
  TH1* getHist1D(const char* hname, const char* htitle, int nbins, float xlow, float xhigh, const char* region, const char* channel);
  TH1* getHist1D(const std::string& hname, const std::string& htitle, int nbins, float xlow, float xhigh, const std::string& region, const std::string& channel);
  template <class T>
    bool fillHist1D(const char* hname, const char* htitle, 
		    T value, int nbins, float xlow, float xhigh, 
		    const std::map<std::string, bool> &regionFlags, 
		    const std::map<std::string, bool> &channelFlags, 
		    double w=1.0) {
    for (auto const& x : channelFlags) {
      std::string channel = x.first;
      bool check1 = x.second;
      for (auto const& y : regionFlags) {
	std::string region = y.first;
	bool check2 = y.second;
	TH1* h = getHist1D(hname, htitle, nbins, xlow, xhigh, region.c_str(), channel.c_str());
	if (h == nullptr) continue;
	if (check1 && check2) h->Fill(value, w);
      }
    }
    return true;
  }

  template <class T>
    bool fillHist1D(const std::string& hname, T value, int nbins, float xlow, float xhigh, const std::map<std::string, bool> regionFlags, 
		    const std::map<std::string, bool> channelFlags, double w=1.0) {
    return fillHist1D(hname.c_str(), value, nbins, xlow, xhigh, regionFlags, channelFlags, w);
  }

  // static approach for 1d histogram booking and filling
  TH1* getHist1D(const char* hname);
  TH1* getHist1D(const std::string& hname);
  template <class T>
    bool fillHist1D(const char* hname, T value, double w=1.0) {
    TH1* h = getHist1D(hname);
    if (h == nullptr) return false;
    h->Fill(value, w);
    return true;
  }
  template <class T>
    bool fillHist1D(const std::string& hname, T value, double w=1.0) {
    return fillHist1D(hname.c_str(), value, w);
  }
  // ---------------------------------------------
  // Convenience routine for filling 2D histograms
  // ---------------------------------------------
  // dynamic approach for 2d histogram booking and filling
  TH2* getHist2D(const char* hname, const char* htitle, 
		 int nbinsX, float xlow, float xhigh, int nbinsY, float ylow, float yhigh, 
		 const char* region, const char* channel);
  TH2* getHist2D(const std::string& hname, const std::string& htitle, 
		 int nbinsX, float xlow, float xhigh, int nbinsY, float ylow, float yhigh, 
		 const std::string& region, const std::string& channel);
  template <class T>
    bool fillHist2D(const char* hname, const char* htitle, 
		    T xvalue, T yvalue, int nbinsX, float xlow, float xhigh, int nbinsY, float ylow, float yhigh,
		    const std::map<std::string, bool> &regionFlags, 
		    const std::map<std::string, bool> &channelFlags, 
		    double w=1.0) {
    for (auto const& x : channelFlags) {
      std::string channel = x.first;
      bool check1 = x.second;
      for (auto const& y : regionFlags) {
	std::string region = y.first;
	bool check2 = y.second;
	TH1* h = getHist2D(hname, htitle, nbinsX, xlow, xhigh, nbinsY, ylow, yhigh, region.c_str(), channel.c_str());
	if (h == nullptr) continue;
	if (check1 && check2) h->Fill(xvalue, yvalue, w);
      }
    }
    return true;
  }
  template <class T>
    bool fillHist2D(const std::string& hname, const std::string& htitle, 
		    T xvalue, T yvalue, int nbinsX, float xlow, float xhigh, int nbinsY, float ylow, float yhigh, 
		    const std::map<std::string, bool> regionFlags, const std::map<std::string, bool> channelFlags, double w=1.0) {
    return fillHist2D(hname.c_str(), htitle.c_str(), xvalue, yvalue, nbinsX, xlow, xhigh, nbinsY, ylow, yhigh, regionFlags, channelFlags, w);
  }

  // static approach for 2d histogram booking and filling
  TH2* getHist2D(const char* hname);
  TH2* getHist2D(const std::string& hname);
  template <class T>
    bool fillHist2D(const char* hname, T xvalue, T yvalue, double w=1.0) {
    TH2* h = getHist2D(hname);
    if (h == nullptr) return false;
    h->Fill(xvalue, yvalue, w);
    return true;
  }
  template <class T>
    bool fillHist2D(const std::string& hname, T xvalue, T yvalue, double w=1.0) {
    return fillHist2D(hname.c_str(), xvalue, yvalue, w);
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
  void SetEvtCutFlowBinLabels(const std::string& hname, 
			      const std::vector<std::string>& slist);
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
