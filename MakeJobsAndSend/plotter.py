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
from subprocess import Popen, PIPE
from collections import defaultdict
import logging
logging.basicConfig(level=logging.DEBUG,format='%(asctime)s - %(levelname)s - %(message)s',datefmt='%m/%d/%Y %H:%M:%S')


class plotter:
    def __init__(self, fileDict, xsecDict, histYaml):
        self.fileDict  = fileDict
        self.xsecDict  = xsecDict
        self.histYaml  = histYaml

    def loadHistYaml(self):
        with open(self.histYaml,'r') as file:
            return yaml.safe_load(file)

    def getListOfHistograms(self):
        with open(self.histYaml,'r') as file:
            config=yaml.safe_load(file)
        histList = config.get('Histograms')
        return histList
    '''
    def getListOfHistograms(self):
        index = 0
        histList = []
        for group, fileList in self.refDict.items():
            if len(fileList) == 0:
                continue
            elif index == 1:
                break
            else:
                file = fileList[0]
                print('Getting list of Histograms from : {}'.format(file))
                outfile = ROOT.TFile.Open(file,'r')
                histList.append(outfile.ls())
                print(histList)
                index=index+1
    '''
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

