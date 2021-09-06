#include <memory>
#include <iostream>

#include "HistBooker.h"

void HistBooker::bookHist1D(const char* hname,
			    const char* htitle, 
			    int nbins, float xlow, float xhigh) {
  new TH1D(hname, htitle, nbins, xlow, xhigh);
}

void HistBooker::bookHist1D(const char* hname,
			    const char* htitle,
			    int nbins, float xlow, float xhigh,
			    const std::vector<std::string> &Flags1,
			    const std::vector<std::string> &Flags2) {
  for (auto &flag1 : Flags1) {
    for(auto &flag2 : Flags2) {
      std::string hname_ = std::string(flag1)+std::string(flag2)+std::string(hname);
      new TH1D(hname_.c_str(), htitle, nbins, xlow, xhigh);
    }
  }
}

void HistBooker::bookHistograms (bool isMC) {
  // book basic histograms to be filled at different stages
  bookHist1D("evtCutFlow",          "Event CutFlow",                                    15, -0.5, 14.5);
  bookHist1D("SR_yield",            "Yield in Signal Region",                             9, -0.5, 8.5);
  bookHist1D("SB_yield",            "Yield in Fake Extrapolated Region",                  9, -0.5, 8.5);
  if (isMC) {
    bookHist1D("EventWtSum",        "Event Weight Sum",                                   1, -0.5, 0.5);
    bookHist1D("evtCutFlowWt",      "Event CutFlow (Weighted)",                         15, -0.5, 14.5);
    bookHist1D("SR_yieldWt",        "Yield in Signal Region Weighted",                    9, -0.5, 8.5);
    bookHist1D("SB_yieldWt",        "Yield in Fake Extrapolated Region Weighted",         9, -0.5, 8.5);
  }

  std::vector<std::string>RegionFlags_ = {"_SR_","_SB_"};
  std::vector<std::string>ChannelFlags_ = {"EleEle","EleMu", "MuMu"};
  // Book any hists with basic channel and region flags if needed.

  std::vector<std::string>ResolvedFlags_ = AnaUtil::mergeToList(RegionFlags_, "IsResolved_");
  // MET plots
  bookHist1D("MetPt",              "Missing E_{T} (GeV)",                                  50, 0, 300,      ChannelFlags_, ResolvedFlags_);
  bookHist1D("MetPhi",             "Missing E_{T} #Phi",                                   32, -3.2, 3.2,   ChannelFlags_, ResolvedFlags_);
  // Di-lepton plots
  bookHist1D("Lep1pt",             "Leading lepton p_{T} (GeV)",                           50, 0, 300,      ChannelFlags_, ResolvedFlags_);
  bookHist1D("Lep2pt",             "Sub-leading lepton p_{T} (GeV)",                       50, 0, 300,      ChannelFlags_, ResolvedFlags_);
  bookHist1D("DiLepPt",            "Di-lepton p_{T} (GeV)",                                50, 0, 300,      ChannelFlags_, ResolvedFlags_);
  bookHist1D("Lep1eta",            "Leading lepton #eta",                                  25, -2.5, 2.5,   ChannelFlags_, ResolvedFlags_);
  bookHist1D("Lep2eta",            "Sub-leading lepton #eta",                              25, -2.5, 2.5,   ChannelFlags_, ResolvedFlags_);
  bookHist1D("Lep1phi",            "Leading lepton #Phi",                                  32, -3.2, 3.2,   ChannelFlags_, ResolvedFlags_);
  bookHist1D("Lep2phi",            "Sub-leading lepton #Phi",                              50, -3.2, 3.2,   ChannelFlags_, ResolvedFlags_);
  bookHist1D("InvM_l1l2",          "Invariant mass [lep1, lep2] (GeV)",                    50, 0, 200,      ChannelFlags_, ResolvedFlags_);
  bookHist1D("DR_lep1lep2",        "#Delta R [lep1, lep2]",                                25, 0, 5,        ChannelFlags_, ResolvedFlags_);
  bookHist1D("DPhi_lep1lep2",      "#Delta#Phi [lep1, lep2]",                              32, -3.2, 3.2,   ChannelFlags_, ResolvedFlags_);
  bookHist1D("DEta_lep1lep2",      "#Delta#eta [lep1, lep2]",                              25, -2.5, 2.5,   ChannelFlags_, ResolvedFlags_);
  // Ak4-Jets plots
  bookHist1D("NoAk4Jets",          "Number of Ak4 jets",                                   10, 0, 10,       ChannelFlags_, ResolvedFlags_);
  bookHist1D("Ak4Jet1Pt",          "Leading Ak4 jet p_{T} (GeV)",                          20, 0, 300,      ChannelFlags_, ResolvedFlags_);
  bookHist1D("Ak4Jet2Pt",          "2ndLeading Ak4 jet p_{T} (GeV)",                       20, 0, 300,      ChannelFlags_, ResolvedFlags_);
  bookHist1D("Ak4Jet3Pt",          "3rdLeading Ak4 jet p_{T} (GeV)",                       20, 0, 300,      ChannelFlags_, ResolvedFlags_);
  bookHist1D("Ak4Jet4Pt",          "4thLeading Ak4 jet p_{T} (GeV)",                       20, 0, 300,      ChannelFlags_, ResolvedFlags_);
  bookHist1D("Ak4Jet1Eta",         "Leading Ak4 jet #eta",                                 25, -5, 5,       ChannelFlags_, ResolvedFlags_);
  bookHist1D("Ak4Jet2Eta",         "2ndLeading Ak4 jet #eta",                              25, -5, 5,       ChannelFlags_, ResolvedFlags_);
  bookHist1D("Ak4Jet3Eta",         "3rdLeading Ak4 jet #eta",                              25, -5, 5,       ChannelFlags_, ResolvedFlags_);
  bookHist1D("Ak4Jet4Eta",         "4thLeading Ak4 jet #eta",                              25, -5, 5,       ChannelFlags_, ResolvedFlags_);
  bookHist1D("Ak4Jet1Phi",         "Leading Ak4 jet #Phi",                                 28, -3.5, 3.5,   ChannelFlags_, ResolvedFlags_);
  bookHist1D("Ak4Jet2Phi",         "2ndLeading Ak4 jet #Phi",                              28, -3.5, 3.5,   ChannelFlags_, ResolvedFlags_);
  bookHist1D("Ak4Jet3Phi",         "3rdLeading Ak4 jet #Phi",                              28, -3.5, 3.5,   ChannelFlags_, ResolvedFlags_);
  bookHist1D("Ak4Jet4Phi",         "4thLeading Ak4 jet #Phi",                              28, -3.5, 3.5,   ChannelFlags_, ResolvedFlags_);
  
  bookHist1D("Ak4Jet1btagDeepFlvB","bTagDeepFlvScore of 1st Ak4 jet",                      20, 0, 1,        ChannelFlags_, ResolvedFlags_);
  bookHist1D("Ak4Jet2btagDeepFlvB","bTagDeepFlvScore of 2nd Ak4 jet",                      20, 0, 1,        ChannelFlags_, ResolvedFlags_);
  bookHist1D("Ak4Jet3btagDeepFlvB","bTagDeepFlvScore of 3rd Ak4 jet",                      20, 0, 1,        ChannelFlags_, ResolvedFlags_);
  bookHist1D("Ak4Jet4btagDeepFlvB","bTagDeepFlvScore of 4th Ak4 jet",                      20, 0, 1,        ChannelFlags_, ResolvedFlags_);

  bookHist1D("InvM_jet1_jet2",     "Invariant mass [jet1, jet2] (GeV)",                    100, 0, 500,     ChannelFlags_, ResolvedFlags_);
  bookHist1D("InvM_jet1_jet3",     "Invariant mass [jet1, jet3] (GeV)",                    100, 0, 500,     ChannelFlags_, ResolvedFlags_);
  bookHist1D("InvM_jet1_jet4",     "Invariant mass [jet1, jet4] (GeV)",                    100, 0, 500,     ChannelFlags_, ResolvedFlags_);
  bookHist1D("InvM_jet2_jet3",     "Invariant mass [jet2, jet3] (GeV)",                    100, 0, 500,     ChannelFlags_, ResolvedFlags_);
  bookHist1D("InvM_jet2_jet4",     "Invariant mass [jet2, jet4] (GeV)",                    100, 0, 500,     ChannelFlags_, ResolvedFlags_);
  bookHist1D("InvM_jet3_jet4",     "Invariant mass [jet3, jet4] (GeV)",                    100, 0, 500,     ChannelFlags_, ResolvedFlags_);
  bookHist1D("DR_jet1_jet2",       "#Delta R [jet1, jet2]",                                25, 0, 5,        ChannelFlags_, ResolvedFlags_);
  bookHist1D("DR_jet1_jet3",       "#Delta R [jet1, jet3]",                                25, 0, 5,        ChannelFlags_, ResolvedFlags_);
  bookHist1D("DR_jet1_jet4",       "#Delta R [jet1, jet4]",                                25, 0, 5,        ChannelFlags_, ResolvedFlags_);
  bookHist1D("DR_jet2_jet3",       "#Delta R [jet2, jet3]",                                25, 0, 5,        ChannelFlags_, ResolvedFlags_);
  bookHist1D("DR_jet2_jet4",       "#Delta R [jet2, jet4]",                                25, 0, 5,        ChannelFlags_, ResolvedFlags_);
  bookHist1D("DR_jet3_jet4",       "#Delta R [jet3, jet4]",                                25, 0, 5,        ChannelFlags_, ResolvedFlags_);
  bookHist1D("DPhi_jet1_jet2",     "#Delta#Phi [jet1, jet2]",                              32, -3.2, 3.2,   ChannelFlags_, ResolvedFlags_);
  bookHist1D("DPhi_jet1_jet3",     "#Delta#Phi [jet1, jet3]",                              32, -3.2, 3.2,   ChannelFlags_, ResolvedFlags_);
  bookHist1D("DPhi_jet1_jet4",     "#Delta#Phi [jet1, jet4]",                              32, -3.2, 3.2,   ChannelFlags_, ResolvedFlags_);
  bookHist1D("DPhi_jet2_jet3",     "#Delta#Phi [jet2, jet3]",                              32, -3.2, 3.2,   ChannelFlags_, ResolvedFlags_);
  bookHist1D("DPhi_jet2_jet4",     "#Delta#Phi [jet2, jet4]",                              32, -3.2, 3.2,   ChannelFlags_, ResolvedFlags_);
  bookHist1D("DPhi_jet3_jet4",     "#Delta#Phi [jet3, jet4]",                              32, -3.2, 3.2,   ChannelFlags_, ResolvedFlags_);
  bookHist1D("DEta_jet1_jet2",     "#Delta#eta [jet1, jet2]",                              50, -5, 5,       ChannelFlags_, ResolvedFlags_);
  bookHist1D("DEta_jet1_jet3",     "#Delta#eta [jet1, jet3]",                              50, -5, 5,       ChannelFlags_, ResolvedFlags_);
  bookHist1D("DEta_jet1_jet4",     "#Delta#eta [jet1, jet4]",                              50, -5, 5,       ChannelFlags_, ResolvedFlags_);
  bookHist1D("DEta_jet2_jet3",     "#Delta#eta [jet2, jet3]",                              50, -5, 5,       ChannelFlags_, ResolvedFlags_);
  bookHist1D("DEta_jet2_jet4",     "#Delta#eta [jet2, jet4]",                              50, -5, 5,       ChannelFlags_, ResolvedFlags_);
  bookHist1D("DEta_jet3_jet4",     "#Delta#eta [jet3, jet4]",                              50, -5, 5,       ChannelFlags_, ResolvedFlags_);

  bookHist1D("JetsInvMass",        "Invariant mass of max 3 jets (GeV)",                   100, 0, 500,     ChannelFlags_, ResolvedFlags_);
  bookHist1D("HT",                 "Jets scalar sum p_{T} (GeV)",                          50, 0, 700,      ChannelFlags_, ResolvedFlags_);
  bookHist1D("HT_vectSum",         "Jets vector sum p_{T} (GeV)",                          50, 0, 700,      ChannelFlags_, ResolvedFlags_);

  // HL plots
  bookHist1D("DR_lep1jet1",        "#Delta R [lep1, jet1]",                                25, 0, 5,        ChannelFlags_, ResolvedFlags_);
  bookHist1D("DR_lep1jet2",        "#Delta R [lep1, jet2]",                                25, 0, 5,        ChannelFlags_, ResolvedFlags_);
  bookHist1D("DR_lep1jet3",        "#Delta R [lep1, jet3]",                                25, 0, 5,        ChannelFlags_, ResolvedFlags_);
  bookHist1D("DR_lep1jet4",        "#Delta R [lep1, jet4]",                                25, 0, 5,        ChannelFlags_, ResolvedFlags_);

  bookHist1D("DEta_lep1jet1",      "#Delta#eta [lep1, jet1]",                              50, -5, 5,       ChannelFlags_, ResolvedFlags_);
  bookHist1D("DEta_lep1jet2",      "#Delta#eta [lep1, jet2]",                              50, -5, 5,       ChannelFlags_, ResolvedFlags_);
  bookHist1D("DEta_lep1jet3",      "#Delta#eta [lep1, jet3]",                              50, -5, 5,       ChannelFlags_, ResolvedFlags_);
  bookHist1D("DEta_lep1jet4",      "#Delta#eta [lep1, jet4]",                              50, -5, 5,       ChannelFlags_, ResolvedFlags_);

  bookHist1D("DPhi_lep1jet1",      "#Delta#Phi [lep1, jet1]",                              32, -3.2, 3.2,   ChannelFlags_, ResolvedFlags_);
  bookHist1D("DPhi_lep1jet2",      "#Delta#Phi [lep1, jet2]",                              32, -3.2, 3.2,   ChannelFlags_, ResolvedFlags_);
  bookHist1D("DPhi_lep1jet3",      "#Delta#Phi [lep1, jet3]",                              32, -3.2, 3.2,   ChannelFlags_, ResolvedFlags_);
  bookHist1D("DPhi_lep1jet4",      "#Delta#Phi [lep1, jet4]",                              32, -3.2, 3.2,   ChannelFlags_, ResolvedFlags_);

  bookHist1D("DR_lep2jet1",        "#Delta R [lep2, jet1]",                                25, 0, 5,        ChannelFlags_, ResolvedFlags_);
  bookHist1D("DR_lep2jet2",        "#Delta R [lep2, jet2]",                                25, 0, 5,        ChannelFlags_, ResolvedFlags_);
  bookHist1D("DR_lep2jet3",        "#Delta R [lep2, jet3]",                                25, 0, 5,        ChannelFlags_, ResolvedFlags_);
  bookHist1D("DR_lep2jet4",        "#Delta R [lep2, jet4]",                                25, 0, 5,        ChannelFlags_, ResolvedFlags_);

  bookHist1D("DEta_lep2jet1",      "#Delta#eta [lep2, jet1]",                              50, -5, 5,       ChannelFlags_, ResolvedFlags_);
  bookHist1D("DEta_lep2jet2",      "#Delta#eta [lep2, jet2]",                              50, -5, 5,       ChannelFlags_, ResolvedFlags_);
  bookHist1D("DEta_lep2jet3",      "#Delta#eta [lep2, jet3]",                              50, -5, 5,       ChannelFlags_, ResolvedFlags_);
  bookHist1D("DEta_lep2jet4",      "#Delta#eta [lep2, jet4]",                              50, -5, 5,       ChannelFlags_, ResolvedFlags_);

  bookHist1D("DPhi_lep2jet1",      "#Delta#Phi [lep2, jet1]",                              32, -3.2, 3.2,   ChannelFlags_, ResolvedFlags_);
  bookHist1D("DPhi_lep2jet2",      "#Delta#Phi [lep2, jet2]",                              32, -3.2, 3.2,   ChannelFlags_, ResolvedFlags_);
  bookHist1D("DPhi_lep2jet3",      "#Delta#Phi [lep2, jet3]",                              32, -3.2, 3.2,   ChannelFlags_, ResolvedFlags_);
  bookHist1D("DPhi_lep2jet4",      "#Delta#Phi [lep2, jet4]",                              32, -3.2, 3.2,   ChannelFlags_, ResolvedFlags_);

  bookHist1D("ST",                 "Jets+leptons+Missing E_{T} scalar sum p_{T} (GeV)",    50, 0, 700,      ChannelFlags_, ResolvedFlags_);


  std::vector<std::string>BoostedFlags_ = AnaUtil::mergeToList(RegionFlags_, "IsBoosted_");
  // MET plots
  bookHist1D("MetPt",              "Missing E_{T} (GeV)",                                  50, 0, 300,      ChannelFlags_, BoostedFlags_);
  bookHist1D("MetPhi",             "Missing E_{T} #Phi",                                   32, -3.2, 3.2,   ChannelFlags_, BoostedFlags_);
  // Di-lepton plots
  bookHist1D("Lep1pt",             "Leading lepton p_{T} (GeV)",                           50, 0, 300,      ChannelFlags_, BoostedFlags_);
  bookHist1D("Lep2pt",             "Sub-leading lepton p_{T} (GeV)",                       50, 0, 300,      ChannelFlags_, BoostedFlags_);
  bookHist1D("DiLepPt",            "Di-lepton p_{T} (GeV)",                                50, 0, 300,      ChannelFlags_, BoostedFlags_);
  bookHist1D("Lep1eta",            "Leading lepton #eta",                                  25, -2.5, 2.5,   ChannelFlags_, BoostedFlags_);
  bookHist1D("Lep2eta",            "Sub-leading lepton #eta",                              25, -2.5, 2.5,   ChannelFlags_, BoostedFlags_);
  bookHist1D("Lep1phi",            "Leading lepton #Phi",                                  32, -3.2, 3.2,   ChannelFlags_, BoostedFlags_);
  bookHist1D("Lep2phi",            "Sub-leading lepton #Phi",                              50, -3.2, 3.2,   ChannelFlags_, BoostedFlags_);
  bookHist1D("InvM_l1l2",          "Invariant mass [lep1, lep2] (GeV)",                    50, 0, 200,      ChannelFlags_, BoostedFlags_);
  bookHist1D("DR_lep1lep2",        "#Delta R [lep1, lep2]",                                25, 0, 5,        ChannelFlags_, BoostedFlags_);
  bookHist1D("DPhi_lep1lep2",      "#Delta#Phi [lep1, lep2]",                              32, -3.2, 3.2,   ChannelFlags_, BoostedFlags_);
  bookHist1D("DEta_lep1lep2",      "#Delta#eta [lep1, lep2]",                              25, -2.5, 2.5,   ChannelFlags_, BoostedFlags_);
  // Ak8 jets plots
  bookHist1D("Ak8Jet1Pt",          "Leading Ak8 jet p_{T} (GeV)",                          25, 0, 400,      ChannelFlags_, BoostedFlags_);
  bookHist1D("Ak8Jet2Pt",          "Sub-leading Ak8 jet p_{T} (GeV)",                      25, 0, 400,      ChannelFlags_, BoostedFlags_);
  bookHist1D("Ak8Jet1Eta",         "Leading Ak8 jet #eta",                                 50, -5, 5,       ChannelFlags_, BoostedFlags_);
  bookHist1D("Ak8Jet2Eta",         "Sub-leading Ak8 jet #eta",                             50, -5, 5,       ChannelFlags_, BoostedFlags_);
  bookHist1D("Ak8Jet1Phi",         "Leading Ak8 jet #Phi",                                 32, -3.2, 3.2,   ChannelFlags_, BoostedFlags_);
  bookHist1D("Ak8Jet2Phi",         "Sub-leading Ak8 jet #Phi",                             32, -3.2, 3.2,   ChannelFlags_, BoostedFlags_);

  // Ak4 jets plots
  bookHist1D("NoAk4Jets", "No. of ak4 jets", 10, 0, 10, ChannelFlags_, BoostedFlags_);
  bookHist1D("Ak4Jet1Pt", "Leading ak4 jet p_{T} (GeV)", 50, 0, 300, ChannelFlags_, BoostedFlags_);
  bookHist1D("Ak4Jet2Pt", "Sub-leading ak4 jet p_{T} (GeV)", 50, 0, 300, ChannelFlags_, BoostedFlags_);
  bookHist1D("NoAk4JetsHas1FatJet", "No. of ak4 jets with one ak8", 10, 0, 10, ChannelFlags_, BoostedFlags_);
  bookHist1D("NoAk4JetsHas2orMoreFatJet", "No. of ak4 jets with ak8", 10, 0, 10, ChannelFlags_, BoostedFlags_);
}
