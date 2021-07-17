#####################################################
###           A script to produce plots           ###
###                 Gourab Saha                   ###
#####################################################
import os
import yaml
import json
import argparse
import shutil
import fileinput
import sys
import ROOT
import stat
import copy
from subprocess import Popen, PIPE
from collections import defaultdict
import logging
logging.basicConfig(level=logging.DEBUG,format='%(asctime)s - %(levelname)s - %(message)s',datefmt='%m/%d/%Y %H:%M:%S')


class plotter:
    def __init__(self, era, lumi, fileDict, xsecDict, plotYaml):
        self.era       = era
        self.lumi      = lumi
        self.fileDict  = fileDict
        self.xsecDict  = xsecDict
        self.plotYaml  = plotYaml

    def openPlotYaml(self):
        with open(self.plotYaml,'w') as file:
            return yaml.safe_load(file)

    def getListOfHistograms(self):
        index = 0
        histList = []
        for group, fileList in self.fileDict.items():
            if len(fileList) == 0:
                continue
            elif index == 1:
                break
            else:
                file = fileList[0]
                print('Getting list of Histograms from : {}'.format(file))
                outfile = ROOT.TFile(file,"READ")
                histList = [key.GetName() for key in outfile.GetListOfKeys()]
                #histList.append(outfile.ls())
                print(histList)
                index=index+1
    
    def doStackPlots(self, histList):
        for histName in histList:
            print(histName)
            for group, fileList in self.fileDict.items(): 
                print(group)
                for file in fileList:
                    file_ = ROOT.TFile.Open(file,'r')
                    hist_ = file_.Get(histName)
                    xsec = self.xsecDict.get(file)
                    print(file, hist_.GetEntries(), xsec)
        
    def doNormPlots(self, histList):
        for histName in histList:
            print(histName)
            for group, fileList in self.fileDict.items(): 
                print(group)
                for file in fileList:
                    file_ = ROOT.TFile.Open(file,'r')
                    hist_ = file_.Get(histName)
                    

                    xsec = self.xsecDict.get(file)
                    print(file, hist_.GetEntries(), xsec)

    def getCommonInfoDict(self):
        return {
            "configuration": {
                "blinded-range-fill-color": "#FDFBFB",
                "blinded-range-fill-style": 4050,
                "eras": [self.era],
                "experiment": "CMS",
                "extra-label": "Preliminary results",
                "height": 600,
                "luminosity": {
                    str(self.era): self.lumi
                },
                "luminosity-label": "%1$.2f fb^{-1} (13 TeV)",
                "margin-bottom": 0.13,
                "margin-left": 0.15,
                "margin-right": 0.03,
                "margin-top": 0.05,
                "root": "results",
                "show-overflow": "true",
                "width": 800,
                "yields-table-align": "v"
            }
        }
        
    def getFileInfoDict(self):
        fileDictToYaml = dict()
        fileInfoDict   = dict()

        for group, fileList in self.fileDict.items():
            testDict = dict()
            for file in fileList:
                idx = 1
                fileName = os.path.basename(file)
                xsec = self.xsecDict.get(file)[0]
                datatype = self.xsecDict.get(file)[1].lower()
                if not datatype == 'data' :
                    rootFile     = ROOT.TFile(file,"READ") 
                    EvtWtSumHist = copy.deepcopy(rootFile.Get('EventWtSum'))
                    print(type(EvtWtSumHist))
                    rootFile.Close()
                    generatedEvents = EvtWtSumHist.GetBinContent(EvtWtSumHist.FindBin(0))
                    
                    fileInfo = {
                        "cross-section":xsec, 
                        "era": self.era,
                        "generated-events":generatedEvents,
                        "group": group,
                        "type": datatype
                    }
                    if 'FakeExtrapolation' in fileName:
                        fileInfo = {
                            "branching-ratio": -1.0,
                            "cross-section":xsec, 
                            "era": self.era,
                            "generated-events":generatedEvents,
                            "group": "Fake",
                            "type": datatype
                        }
                    elif datatype == 'signal':
                        fileInfo = {
                            "branching-ratio": 10.0,
                            "cross-section":xsec, 
                            "era": self.era,
                            "generated-events":generatedEvents,
                            "legend": fileName.split('.')[0],
                            "line-color":'#8f0a1e',
                            "line-type": idx,
                            "order": idx,
                            "type": datatype
                        }
                        idx = idx + 1
                else :
                    fileInfo = {
                        "era": self.era,
                        "group": datatype,
                        "type": datatype
                    }
                    if 'FakeExtrapolation' in fileName:
                        fileInfo = {
                            "cross-section": 1.0,
                            "era": self.era,
                            "group": "Fake",
                            "type": datatype
                        }
                    
                testDict[fileName] = fileInfo   
                fileInfoDict.update(testDict)
        fileDictToYaml['files'] = fileInfoDict
        return fileDictToYaml

    def getHistogramDict(self):
        pass
        # TODO :: modify and finish the input to plot-It
