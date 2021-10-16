#####################################################
###      A script to produce yaml for plots       ###
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
#logging.basicConfig(level=logging.DEBUG,format='%(asctime)s - %(levelname)s - %(message)s',datefmt='%m/%d/%Y %H:%M:%S')
logger = logging.getLogger("postlog")

class yamlMaker:
    def __init__(self, era, lumi, fileDict, xsecDict, legendDict, histDir, isNorm=False):
        self.era        = era
        self.lumi       = lumi
        self.fileDict   = fileDict
        self.xsecDict   = xsecDict
        self.legendDict = legendDict
        self.histDir    = histDir
        self.isNorm     = isNorm

    def getListOfHistograms(self):
        killhistos = ['FakeExtrapolation','ObjectSelection','evtCutFlow','EventWtSum','yield']
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
                if len(listOfKeys) == 0:
                    print(f'WARNING! {file} doesnt have any histograms !!!')
                    continue
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
            extra-label: "Work in progress"
            height: 600
            #lumi is a dict {str(self.era): self.lumi}
            luminosity: 
            luminosity-label: "%1$.2f fb^{-1} (13 TeV)"
            margin-bottom: 0.13
            margin-left: 0.15
            margin-right: 0.03
            margin-top: 0.05
            root: 
            show-overflow: true
            width: 800
            yields-table-align: v
            '''
        )
        commonDict['eras'] = [self.era]
        commonDict['luminosity'] = {self.era : self.lumi}
        commonDict['root'] = str(os.path.join(self.histDir, 'results'))

        configDict['configuration'] = commonDict
        #print(yaml.dump(configDict,default_flow_style=False))
        return configDict

    '''
    def getLegendInfoDict(self, nofake=False):
        import random
        groupLegendDict = dict()
        groupLegendInfo = dict()
        for group, fileList in self.fileDict.items():
            if nofake : 
                if group == 'Fake':
                    continue
            if 'Signal' in group:
                continue
            if group == 'data':
                groupLegendInfo[group] = {
                    "legend":'data'
                }
            else:
                random_number = random.randint(1100000,16777215)
                hex_number = str(hex(random_number))
                hex_number ='#'+ hex_number[2:]
                stack_order = random.randint(1,15)
                groupLegendInfo[group] = {
                    "fill-color": hex_number,
                    "legend": self.legendDict.get(group),
                    "order": stack_order
                }
        groupLegendDict['groups'] = groupLegendInfo
        return groupLegendDict
    '''
    def getLegendInfoDict(self):
        groupLegendDict = dict()
        groupLegendInfo = dict()
        for group, leginfoList in self.legendDict.items():
            if not group == 'data':
                groupLegendInfo[group] = {
                    "fill-color": leginfoList[0],
                    "legend": leginfoList[1],
                    "order": leginfoList[2]
                }
            else:
                groupLegendInfo[group] = {
                    "legend": leginfoList[0]
                }
        groupLegendDict['groups'] = groupLegendInfo
        return groupLegendDict


    def getFileInfoDict(self, lumi, nofake=False):
        fileDictToYaml = dict()
        fileInfoDict   = dict()

        for group, fileList in self.fileDict.items():
            testDict = dict()
            for file in fileList:
                idx = 1
                fileName = os.path.basename(file)
                fileHandle = file.replace('_FakeExtrapolation','') if '_FakeExtrapolation' in fileName else file
                xsec = self.xsecDict.get(fileHandle)[0]
                datatype = self.xsecDict.get(fileHandle)[1].lower()
                if nofake:
                    if 'FakeExtrapolation' in fileName:
                        continue
                #if not datatype == 'data' :
                if datatype == 'mc' :
                    rootFile     = ROOT.TFile(fileHandle,"READ") 
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
                            "group": group,
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
                            "generated-events": lumi,
                            "group": group,
                            "type": 'mc'
                        }
                        
                testDict[fileName] = fileInfo   
                fileInfoDict.update(testDict)
        fileDictToYaml['files'] = fileInfoDict
        #print(yaml.dump(fileDictToYaml,default_flow_style=False))
        return fileDictToYaml
    
    @staticmethod
    def plotInfo_template():
        return yaml.safe_load(
            '''
            #blinded-range:
            labels:
            log-y: both
            ratio-y-axis: '#frac{Data}{MC}'
            ratio-y-axis-range:
            show-ratio: true
            sort-by-yields: true
            x-axis:
            x-axis-range:
            y-axis: Events
            #y-axis-range: 
            y-axis-show-zero: true
            legend-columns: 3
            normalized:
            no-data: false
            save-extensions: [pdf]
            '''
        )

    def getHistogramDict(self, histList, titleList):
        plotYamlObj = self.plotInfo_template()
        #print(f'\n===>> templated_plot_yaml_object : \n{yaml.dump(plotYamlObj,default_flow_style=False)}')
        histInfoDictToYaml= dict()
        histInfoDict = dict()
        histDict = dict()
        for i,hist in enumerate(histList):
            histDict[hist] = dict(plotYamlObj)
            if 'EleEle' in hist:
                label = [{'text': 'e^{#pm}e^{#mp}', 'size': 36, 'position': [0.20,0.87]}]
            elif 'EleMu' in hist:
                label = [{'text': 'e^{#pm}#mu^{#mp}', 'size': 36, 'position': [0.20,0.87]}]
            elif 'MuMu' in hist:
                label = [{'text': '#mu^{#pm}#mu^{#mp}', 'size': 36, 'position': [0.20,0.87]}]
            elif 'Ele' in hist:
                label = [{'text': 'e^{#pm}', 'size': 36, 'position': [0.20,0.87]}]
            elif 'Mu' in hist:
                label = [{'text': '#mu^{#pm}', 'size': 36, 'position': [0.20,0.87]}]
            else:
                label = [{'text': 'All channel', 'size': 36, 'position': [0.20,0.87]}]

            histDict[hist]['labels'] = label
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
                    yHigh        = histogram.GetMaximum()
                    #print(idx, file, xLow, xHigh)
                    rootFile.Close()
                    if (idx > 0):
                        break
                if (idx > 0):
                    break
            histDict[hist]['x-axis-range'] = [xLow, xHigh]
            #histDict[hist]['y-axis-range'] = [0., 1.25*yHigh]
            histDict[hist]['ratio-y-axis-range'] = [0.5, 1.5]
            histDict[hist]['normalized'] = self.isNorm
            histInfoDict.update(histDict)
        #print(histInfoDict)
        histInfoDictToYaml['plots'] = histInfoDict
        #print(yaml.dump(histInfoDictToYaml,default_flow_style=False))
        return histInfoDictToYaml
