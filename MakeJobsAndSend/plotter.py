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
    def __init__(self, era, lumi, fileDict, xsecDict, plotYaml, histDir):
        self.era       = era
        self.lumi      = lumi
        self.fileDict  = fileDict
        self.xsecDict  = xsecDict
        self.plotYaml  = plotYaml
        self.histDir   = histDir

    def openPlotYaml(self):
        with open(self.plotYaml,'w') as file:
            return yaml.safe_load(file)

    def getListOfHistograms(self):
        killhistos = ['FakeExtrapolation','ObjectSelection','evtCutFlow','evtCutFlowWt','EventWtSum']
        lambda_matched = lambda sList, name : [True for s in sList if s in name]
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
                listOfKeys = [key for key in outfile.GetListOfKeys() if not any(lambda_matched(killhistos, key.GetName()))]
                histNameList  = [key.GetName() for key in listOfKeys]
                histTitleList = [key.GetTitle() for key in listOfKeys]
                index=index+1

        return histNameList, histTitleList

    def getCommonInfoDict(self):
        configDict = dict()
        commonDict = yaml.safe_load(
            '''
            blinded-range-fill-color: "#FDFBFB"
            blinded-range-fill-style: 4050
            #put list of eras
            eras: 
            experiment: "CMS"
            extra-label: "Preliminary results"
            height: 600
            #lumi is a dict {str(self.era): self.lumi}
            luminosity: 
            luminosity-label: "%1$.2f fb^{-1} (13 TeV)"
            margin-bottom: 0.13
            margin-left: 0.15
            margin-right: 0.03
            margin-top: 0.05
            root:
            show-overflow: "true"
            width: 800
            yields-table-align: "v"
            '''
        )
        commonDict['eras'] = [self.era]
        commonDict['luminosity'] = {self.era : self.lumi}
        commonDict['root'] = str(self.histDir)

        configDict['configuration'] = commonDict
        print(yaml.dump(configDict,default_flow_style=False))
        return configDict


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
        print(yaml.dump(fileDictToYaml,default_flow_style=False))
        return fileDictToYaml
    
    @staticmethod
    def plotInfo_template():
        return yaml.safe_load(
            '''
            #blinded-range:
            #labels:
            log-y: both
            ratio-y-axis: '#frac{Data}{MC}'
            ratio-y-axis-range:
            show-ratio: true
            sort-by-yields: true
            x-axis:
            x-axis-range:
            y-axis: Events
            y-axis-show-zero: true
            '''
        )

    def getHistogramDict(self, histList, titleList):
        plotYamlObj = self.plotInfo_template()
        print(f'\n===>> templated_plot_yaml_object : \n{yaml.dump(plotYamlObj,default_flow_style=False)}')
        histInfoDictToYaml= dict()
        histInfoDict = dict()
        histDict = dict()
        for i,hist in enumerate(histList):
            histDict[hist] = dict(plotYamlObj)
            histDict[hist]['x-axis'] = titleList[i]
            for group, fileList in self.fileDict.items():
                idx = 0
                #print(fileList)
                for file in fileList:
                    idx = idx + 1
                    rootFile     = ROOT.TFile(file,"READ") 
                    histogram    = copy.deepcopy(rootFile.Get(hist))
                    xLow         = histogram.GetXaxis().GetXmin()
                    xHigh        = histogram.GetXaxis().GetXmax()
                    #print(idx, file, xLow, xHigh)
                    rootFile.Close()
                    if (idx > 0):
                        break
                if (idx > 0):
                    break
            histDict[hist]['x-axis-range'] = [xLow, xHigh]
            histDict[hist]['ratio-y-axis-range'] = [0.5, 1.5]
            histInfoDict.update(histDict)
        #print(histInfoDict)
        histInfoDictToYaml['plots'] = histInfoDict
        print(yaml.dump(histInfoDictToYaml,default_flow_style=False))
        return histInfoDictToYaml
