import os
import shutil
import fileinput
import sys
import ROOT
import stat
import copy
import math
from collections import defaultdict
from prettytable import PrettyTable

def getEfficiencies(ibin, evtCutFlow_hist):
    if ibin == 0:
        abs_eff = 1.0
        rel_eff = 1.0
    else:
        den_abs = evtCutFlow_hist.GetBinContent(evtCutFlow_hist.FindBin(0))
        den_rel = evtCutFlow_hist.GetBinContent(evtCutFlow_hist.FindBin(ibin-1))
        num     = evtCutFlow_hist.GetBinContent(evtCutFlow_hist.FindBin(ibin))

        abs_eff = num/den_abs if den_abs > 0 else math.nan
        rel_eff = num/den_rel if den_rel > 0 else math.nan

    return [abs_eff, rel_eff]

def getYields(cutFlowHistDict, cutFlowXLabels, xsec, era):
    tables = []
    totalObservedResolved = 0
    totalObservedBoosted  = 0
    finalYieldTable = PrettyTable(['Process', 'Resolved', 'Boosted'])
    finalYieldTable.title = f'Yield [era : {era}]'
    labels = list(cutFlowHistDict.keys())
    for label, histList in cutFlowHistDict.items():
        t0 = PrettyTable(['Selections','nEvents','nEvents:Weighted','Relative Efficiency','Absolute Efficiency','Yield'])
        t0.title = f'Process : {label}, X-sec : {xsec} pb'
        for i in range(histList[0].GetNbinsX()):
            effList = getEfficiencies(i, histList[0])
            if not 'DATA' in label:
                t0.add_row([cutFlowXLabels[i], 
                            round(histList[0].GetBinContent(histList[0].FindBin(i)),2), 
                            round(histList[1].GetBinContent(histList[1].FindBin(i)),2), 
                            '%.3f'%effList[1],
                            '%.3f'%effList[0],
                            round(histList[2].GetBinContent(histList[2].FindBin(i)),2)])
            else:
                t0.add_row([cutFlowXLabels[i], 
                            round(histList[0].GetBinContent(histList[0].FindBin(i)),2), 
                            math.nan, 
                            '%.3f'%effList[1],
                            '%.3f'%effList[0],
                            round(histList[0].GetBinContent(histList[0].FindBin(i)),2)])

        tables.append(t0)

        if not 'DATA' in label:
            finalYieldTable.add_row([label, 
                                     round(histList[2].GetBinContent(len(cutFlowXLabels)-1),3), 
                                     round(histList[2].GetBinContent(len(cutFlowXLabels)),3)])
        else:
            totalObservedResolved += histList[0].GetBinContent(len(cutFlowXLabels)-1)
            totalObservedBoosted  += histList[0].GetBinContent(len(cutFlowXLabels))

    finalYieldTable.add_row(['Observed', totalObservedResolved, totalObservedBoosted])
            
    return labels, tables, finalYieldTable

