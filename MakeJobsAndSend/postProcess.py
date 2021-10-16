#####################################################
### A script to postprocess the output root files ###
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
from alive_progress import alive_bar
from time import sleep
import logging
from yamlMaker import yamlMaker
from prettytable import PrettyTable
import datetime as cal
import getYields
import makeFakeFileSeparate

def runShellCmd(cmdList):
    process = Popen(cmdList, stdout=PIPE, stderr=PIPE)
    while True:
        output = process.stdout.readline()
        if process.poll() is not None:
            break
        if output:
            print(output.strip().decode("utf-8"))
    rc = process.poll()


def main():
    #logging.basicConfig(level   = logging.DEBUG,
    #                    format  = '%(asctime)s - %(levelname)s - %(message)s',
    #                    datefmt = '%m/%d/%Y %H:%M:%S')
    #logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(description='Make Jobs and Send')    
    parser.add_argument('--yaml', 
                        action='store', 
                        required=True, type=str, 
                        help='Name of the config')
    parser.add_argument('--histdir', 
                        action='store', 
                        required=True, type=str, 
                        help='just name of the output directory')
    parser.add_argument('--runplotit',
                        action='store_true',
                        required=False,default=False,
                        help="Run PlotIt to produce plots")
    parser.add_argument('--hasskim', 
                        action='store_true',
                        required=False,default=False,
                        help="Use this if skimmed ntuples need to be h-added")
    parser.add_argument('--norm',
                        action='store_true',
                        required=False,default=False,
                        help="To produce stack plots but normalised to data")
    parser.add_argument('--nofake',
                        action='store_true',
                        required=False,default=False,
                        help="To produce normalised plots")
    parser.add_argument('--verbose',
                        action='store_true',
                        required=False,default=False,
                        help="verbosity")
    parser.add_argument('--checkinputs',
                        action='store_true',
                        required=False,default=False,
                        help="Check inputfiles if any of the files is corrupted or not")
    parser.add_argument('--forcehadd',
                        action='store_true',
                        required=False,default=False,
                        help="do hadd even if some files are missing")


    args = parser.parse_args()
    
    with open(args.yaml, 'r') as config:
        configDict = yaml.safe_load(config)

    outdir  = configDict.get('outDir')
    histDir = os.path.join(outdir, args.histdir)
    if not os.path.isdir(histDir):
        raise RuntimeError(f'Couldnt find {histDir}')
    dirKey_condor = '_'.join(args.histdir.split('_')[-2:])

    batchOutDir = os.path.join(histDir, 'batch')
    if not os.path.isdir(batchOutDir):
        raise RuntimeError(f'batch dir not found : {batchOutDir}')

    # Logger
    logFormatter = logging.Formatter("%(asctime)s -- [%(levelname)s] -- %(message)s",
                                     "%Y-%m-%d %H:%M:%S")
    logger       = logging.getLogger("postlog")
    logger.setLevel(logging.DEBUG)
    # logger for console
    consoleHandler = logging.StreamHandler(sys.stdout)
    consoleHandler.setFormatter(logFormatter)
    logger.addHandler(consoleHandler)
    # Adding file-handler for logging
    now    = cal.datetime.now()
    logfiletag = now.strftime("%d%m%Y_%H%M%S")
    fileHandler = logging.FileHandler(os.path.join(histDir, 'postprocess_'+logfiletag+'.log'))
    fileHandler.setFormatter(logFormatter)
    logger.addHandler(fileHandler)

    resultDir = os.path.join(histDir, 'results')
    if os.path.isdir(resultDir):
        logger.info(f'{resultDir} exits !!!')
    else:
        os.mkdir(resultDir)

    skimDir = os.path.join(histDir, 'skims')
    if os.path.isdir(skimDir):
        logger.info(f'{skimDir} exits !!!')
    else:
        os.mkdir(skimDir)

    yieldDir = os.path.join(histDir, 'yields')
    if os.path.isdir(yieldDir):
        logger.info(f'{yieldDir} exits !!!')
    else:
        os.mkdir(yieldDir)

    pwd = os.getcwd()
    logger.info('present working dir : {}'.format(pwd))
    
    keyList = [str(key) for key in configDict.keys()]
    nofake  = args.nofake 

    era            = configDict.get('era')
    lumi           = configDict.get('lumi')

    logger.info('era  : {}'.format(era))
    logger.info('lumi : {} pb-1'.format(lumi))

    jobdir = os.path.join(os.getcwd(), 'JobCards_'+str(era)+'_'+dirKey_condor)

    samplesDict       = configDict.get('samplesDict')
    dataTypes         = [str(sample) for sample in samplesDict.keys()]

    if not os.path.isdir(histDir):
        logger.info('{} : output directory doesnt exist !!! '.format(histDir))

    groupLegendInfoDict = configDict.get('groupLegendInfo')
    legendPosDict = configDict.get('legendPos')
    resultDict = defaultdict(list)
    xsecDict   = dict()
    legendDict = dict()
    logger.info('Doing hadd ===>')
    evtCutFlowDict   = defaultdict(list)
    evtCutFlowLabels = []
    combinedDictForYaml = dict()
    resubmitList = []
    cutflowTableFile   = open(os.path.join(yieldDir, f'CutFlow_{era}.txt'),'w')
    yieldTableFile     = open(os.path.join(yieldDir, f'Yields_{era}.txt'),'w')

    for dataType, valDict in samplesDict.items():
        '''
        e.g. 
        dataType : MC
        valDict:
        TTToSemiLeptonic:
          filedirs: ['/eos/user/g/gsaha3/Exotic/MC_UL2017/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/']
          genEvtWtSum: 'genEventSumw'
          xsec: 365.52
          filesPerJob: 2
          group: 'TT'
        
        '''
        ismc     = False
        isdata   = False
        issignal = False
        if dataType == 'MC' : 
            ismc = True
            hasskim = args.hasskim
            forcehadd = args.forcehadd
        elif dataType == 'DATA':
            isdata = True
            hasskim = False
            forcehadd = False
        elif dataType == 'SIGNAL':
            issignal = True
            hasskim = args.hasskim
            forcehadd = args.forcehadd
        else:
            raise RuntimeError('Unknown datatype. Please mention MC, DATA or SIGNAL in the main yaml')
          
        # if verbose
        if args.verbose:
            logger.info(f'dataType : {dataType}')
            logger.info(f'Looking for skimmed ntuples? {hasskim}')

  
        if valDict == None or len(valDict) == 0:
            logger.info(f'Sorry! No {dataType} samples are present in yaml')
        else:
            logger.info(f'Looking for {dataType} - output root files >>------>')
            logger.debug('{}_Samples : {}'.format(str(dataType), [sample for sample in valDict.keys()]))
            # Loop over valDict
            for key, val in valDict.items():
                '''
                KEY = sample name e.g. TTToSemiLeptonic
                VAL = a dictionary with keys and values 
                filedirs: ['/eos/user/g/gsaha3/Exotic/MC_UL2017/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/']
                genEvtWtSum: 'genEventSumw'
                xsec: 365.52
                filesPerJob: 2
                group: 'TT'
                legend: 'tt+jets'
                '''
                logger.info(' Process : {}'.format(key))
                xsec = val.get('xsec')
                filePathList = val.get('filedirs')
                filesPerJob  = int(val.get('filesPerJob'))
                files        = []
                condorDir = os.path.join(jobdir, key+'_condorJobs')
                for item in filePathList:
                    logger.info('\t {}'.format(item))
                    files += [os.path.join(item,rfile) for rfile in os.listdir(item) if '.root' in rfile]

                logger.info('\t nFiles : {}'.format(len(files)))
                logger.info('\t filesPerJob : {}'.format(filesPerJob))
                infileListPerJob = [files[i:i+filesPerJob] for i in range(0, len(files), filesPerJob)]
                nJobs = len(infileListPerJob)
                logger.info('\t nJobs : {}'.format(nJobs))

                # if verbose or checkinputs
                if args.verbose or args.checkinputs:
                    logger.info('\t   Input INFO')
                    for i,jobfileList in enumerate(infileListPerJob):
                        logger.info(f'\t\t Job : {i}')
                        for ele in jobfileList:
                            logger.info('\t\t\t '+ele)
                            elefile = ROOT.Open(ele)
                            if elefile.IsZombie():
                                logger.error('Z O M B I E. Please check this inputfile, especially for data samples !')
                            elefiletrees = [key.GetName() for key in elefile.GetListOfKeys() if key.GetClassName() == 'TTree']
                            if elefiletrees.count('Events') <= 0:
                                logger.warning('Events tree found missing. Please check this inputfile !')
                            if elefiletrees.count('Runs') <= 0:
                                logger.warning('Runs tree found missing. Please check this inputfile !')

                batchHistDir = os.path.join(batchOutDir, key)
                if not os.path.isdir(batchHistDir):
                    raise RuntimeError(f'{batchHistDir} not found !')

                tobehadd         = []
                posthaddfile     = os.path.join(batchHistDir, str(key)+'_hist.root')
                haddcmd_         = ['hadd', posthaddfile]
                if hasskim:
                    tobehaddskim         = []
                    posthaddfileskim     = os.path.join(skimDir, str(key)+'_skim.root')
                    haddcmdskim_         = ['hadd', posthaddfileskim]

                fileabsent = False
                filesfoundmissing = False
                haddCond = os.path.exists(posthaddfile) and os.path.exists(posthaddfileskim) if hasskim else os.path.exists(posthaddfile)

                if haddCond :
                    logger.info(f'hadded files ---> {posthaddfile} & {posthaddfileskim} already exist!') if hasskim else logger.info(f'hadded file ---> {posthaddfile} already exist!')
                else:
                    #logger.info('\t   Output INFO')
                    #logger.info('\t   h-adding >>--------->')
                    #logger.info(f'\t   hadded output to be produced : {os.path.basename(posthaddfile)} using {nJobs} output root files')
                    with alive_bar(nJobs, title='collecting files to hadd ', length=50, enrich_print=True) as bar:
                        for i in range(nJobs) :
                            fileabsent = False
                            rootfile  = os.path.join(batchHistDir, str(key)+'_'+str(i)+'_hist.root')
                            #logger.info(f'\t\t Job {i} output : {rootfile}')
                            tfile    = ROOT.TFile(rootfile,"READ")
                            if not os.path.exists(rootfile): 
                                logger.warning(f'{rootfile} >>-x-x-x-> Missing')
                                fileabsent = True if not forcehadd else False
                                filesfoundmissing = True if not forcehadd else False
                            elif tfile.IsZombie():
                                logger.warning(f'{rootfile} is a Zombie! Please produce this file again.')
                                fileabsent = True if not forcehadd else False
                                filesfoundmissing = True if not forcehadd else False
                            else:
                                listOfHistos = [key for key in tfile.GetListOfKeys()]
                                if len(listOfHistos) <= 1:
                                    logger.error(f'check the logfile of {rootfile} for any possible segfault')
                                    fileabsent = True if not forcehadd else False
                                    filesfoundmissing = True if not forcehadd else False
                                else:
                                    sleep(0.03)
                                    tobehadd.append(rootfile)
                            tfile.Close()
                            if hasskim:
                                rootfileS = os.path.join(batchHistDir, str(key)+'_'+str(i)+'_skim.root')
                                tfileS    = ROOT.TFile(rootfileS,"READ")
                                if not os.path.exists(rootfileS): 
                                    logger.warning(f'{rootfileS} >>-x-x-x-> Missing')
                                    fileabsent = True if not forcehadd else False
                                    filesfoundmissing = True if not forcehadd else False
                                elif tfileS.IsZombie():
                                    logger.warning(f'{rootfileS} is a Zombie! Please produce this file again.')
                                    fileabsent = True if not forcehadd else False
                                    filesfoundmissing = True if not forcehadd else False
                                else:
                                    sleep(0.03)
                                    tobehaddskim.append(rootfileS)
                                tfileS.Close()
                            if fileabsent:
                                if not forcehadd:
                                    resubmitList.append(['condor_submit', os.path.join(condorDir,str(key)+'_'+str(i)+'.sub')])
                            bar()

                    if len(tobehadd) == 0:
                        logger.error(f'{key} job output root files are not present')
                        sys.exit(f'TERMINATED!!! no post hadd files for {key}. Check the input files and then check the scripts :( ')  
                    elif nJobs == len(tobehadd):
                        haddcmd = haddcmd_ + tobehadd
                        runShellCmd(haddcmd)
                    elif nJobs > len(tobehadd):
                        logger.warning('some output root files are missing ...')
                        if forcehadd:
                            logger.warning('[--forcehadd] doing force-hadd !!!')
                            haddcmd = haddcmd_ + tobehadd
                            runShellCmd(haddcmd)
                    else:
                        logger.warning('Skipping hadd hist files because files are found missing')
                        
                    if hasskim:
                        if len(tobehaddskim) == 0:
                            logger.error(f'{key} job output skim root files are not present')
                            sys.exit(f'TERMINATED!!! no post hadd files for {key}. Check the input files and then check the scripts :( ')  
                        elif nJobs == len(tobehaddskim):
                            haddcmdskim = haddcmdskim_ + tobehaddskim
                            runShellCmd(haddcmdskim)
                        elif nJobs > len(tobehaddskim):
                            logger.warning('some output root files are missing ...')
                            if forcehadd:
                                logger.warning('[--forcehadd] doing force-hadd !!!')
                                haddcmdskim = haddcmdskim_ + tobehaddskim
                                runShellCmd(haddcmdskim)
                        else:
                            logger.warning('Skipping hadd skim files because some files are missing')

                logger.info('hadding done ... |')# hadding done .................
                
                # If all files are present, start collecting event-cut-flow histogram
                if not filesfoundmissing:
                    outrootfile  = ROOT.TFile(posthaddfile,'read')
                    if outrootfile.IsZombie():
                        raiseRuntimeError (f'{posthaddfile} is a Zombie !!!')

                    cutFlowHist  = outrootfile.Get('evtCutFlow')
                    if len(evtCutFlowLabels) == 0:
                        for ibin in range(cutFlowHist.GetNbinsX()):
                            evtCutFlowLabels.append(cutFlowHist.GetXaxis().GetBinLabel(cutFlowHist.FindBin(ibin)))

                    if dataType == 'MC' or dataType == 'SIGNAL':
                        evtCutFlowDict[key].append(cutFlowHist)
                        cutFlowWtHist = outrootfile.Get('evtCutFlowWt')
                        cutFlowWtHist_ls = copy.deepcopy(cutFlowWtHist)
                        evtWtSumHist = outrootfile.Get('EventWtSum')
                        if cutFlowWtHist_ls.Integral() > 0:
                            cutFlowWtHist_ls.Scale(lumi*xsec/evtWtSumHist.GetBinContent(1))                        
                        evtCutFlowDict[key].append(cutFlowWtHist)
                        evtCutFlowDict[key].append(cutFlowWtHist_ls)
                        cutFlowHist.SetDirectory(0)
                        cutFlowWtHist.SetDirectory(0)
                        cutFlowWtHist_ls.SetDirectory(0)

                    elif dataType == 'DATA':
                        evtCutFlowDict[dataType+'_'+key].append(cutFlowHist)
                        cutFlowHist.SetDirectory(0)
                        
                    else:
                        raise RuntimeError('please mention correct datatype : MC or DATA or SIGNAL ...')

                    # making SR and FR files separate
                    SR_File, Fake_File = makeFakeFileSeparate.getSRandFakeRootFiles(posthaddfile, resultDir, noFake=nofake, isSignal=issignal)
                    
                    group = val.get('group')
                    if os.path.exists(SR_File):
                        resultDict[str(group)].append(SR_File)
                        #legendDict[str(group)] = val.get('legend')
                        xsecDict[str(SR_File)]=[xsec, dataType]
                    else:
                        logging.warning('hadded file |{}| is absent'.format(SR_File))
                    if not issignal and os.path.exists(Fake_File):
                        resultDict['Fake'].append(Fake_File)
                        #legendDict['Fake'] = 'Fake'
                        
    # Re-submission
    if len(resubmitList) > 0:
        logger.info('Resubmit the following jobs or use the sh script in jobdir')
        with open(os.path.join(jobdir,'resend.sh'), 'w') as outf:
            outf.write('# !/bin/sh \n\n')
            for cmd in resubmitList:
                print(' '.join(cmd))
                outf.write(' '.join(cmd)+'\n')

        logger.info('After re-running these failed jobs, do proceed again with the post-processing.')
        sys.exit()
    else:
        logger.info('All h-added files are present. So moving further --->')
     
    #Preparing legendDict
    for key, items in groupLegendInfoDict.items():
        group     = key
        legend    = items.get('legend') 
        if not group == 'data':
            fillcolor = items.get('fill-color')
            order     = items.get('order')
            legendDict[str(group)] = [fillcolor, legend, order]
        else:
            legendDict[str(group)] = [legend]
    if not args.nofake:
        legendDict['Fake'] = ['#cc6666', 'Fake', 20]


    #logger.debug(resultDict)
    yamlObj = yamlMaker(era,lumi,resultDict,xsecDict,legendDict,histDir,args.norm)
    histograms, histTitles = yamlObj.getListOfHistograms()
    #logger.info(f'List of histograms : \n {histograms} \n')
    commonInfoDict = yamlObj.getCommonInfoDict()
    logger.info('saving commonInfo to plots yaml object')
    combinedDictForYaml.update(commonInfoDict)

    fileInfoDict = yamlObj.getFileInfoDict(lumi,nofake)
    logger.info('saving fileInfo to plots yaml object')
    combinedDictForYaml.update(fileInfoDict)

    logger.info('saving legendInfo to plots yaml object')
    #groupLegendDict = yamlObj.getLegendInfoDict(nofake)
    groupLegendDict = yamlObj.getLegendInfoDict()
    combinedDictForYaml.update(groupLegendDict)

    logger.info('saving legendPosInfo to plots yaml object')
    combinedDictForYaml.update(legendPosDict)

    #logger.info('saving plotGroupInfo to plots yaml object')
    #combinedDictForYaml.update(groupLegendInfoDict)

    plotInfoDict = yamlObj.getHistogramDict(histograms, histTitles)
    logger.info('saving plotsInfo to plots yaml object')
    combinedDictForYaml.update(plotInfoDict)

    logger.info('preparing plot yml file')
    plotyaml = os.path.join(histDir, 'plots.yml') if nofake else os.path.join(histDir, 'plots_FakeExtrapolation.yml')
    with open(plotyaml, 'w') as ymlfile:
        documents = yaml.dump(combinedDictForYaml, ymlfile, default_flow_style=False)

    if args.runplotit:
        logger.info('running PlotIt ... \n')
        plotDir = os.path.join(histDir, 'Plots') if nofake else os.path.join(histDir, 'Plots_FakeExtrapolation')
        if not os.path.exists(plotDir):
            os.mkdir(plotDir)
        
        if not os.path.exists(plotyaml):
            raise RuntimeError(f'{plotyaml} doesnt exist !!!')
        # run plotIt
        plotItCmds = ['plotIt','-y','-o',plotDir, plotyaml]
        '''
        proc       = Popen(plotItCmds, stdout=PIPE)
        # Poll proc.stdout to show stdout live
        while True:
            output = proc.stdout.readline()
            if proc.poll() is not None:
                break
            if output:
                logging.info(output.strip().decode("utf-8"))
        rc = proc.poll()
        '''
        runShellCmd(plotItCmds)

    # Producing yield tables
    '''
    evtCutFlowLabels.insert(0, 'Process')
    yieldTable = PrettyTable(evtCutFlowLabels)
    yieldTableWt = PrettyTable(evtCutFlowLabels)
    yieldTableWtLs = PrettyTable(evtCutFlowLabels)
    for label, histList in evtCutFlowDict.items():
        entries = []
        hist = histList[0]
        for i in range(hist.GetNbinsX()):
            entries.append(round(hist.GetBinContent(hist.FindBin(i)),3))
        yieldTable.add_row([label]+entries)

        if not label == 'DATA':
            histWt = histList[1]
            histWtLs = histList[2]
            entriesWt = []
            entriesWtLs = []
            for i in range(hist.GetNbinsX()):
                entriesWt.append(round(histWt.GetBinContent(histWt.FindBin(i)),3))
                entriesWtLs.append(round(histWtLs.GetBinContent(histWtLs.FindBin(i)),3))
            yieldTableWt.add_row([label]+entriesWt)
            yieldTableWtLs.add_row([label]+entriesWtLs)

    with open (os.path.join(histDir, 'YieldTable.txt'), 'w') as yfile:
        yfile.write(' -------------- | For unweighted events | ---------------- \n')
        yfile.write(str(yieldTable)+'\n')
        yfile.write(' -------------- | For weighted events | ---------------- \n')
        yfile.write(str(yieldTableWt)+'\n')
        yfile.write(f' -------------- | For weighted events (at {lumi} pb-1) | ---------------- \n')
        yfile.write(str(yieldTableWtLs)+'\n')

    yieldTable = PrettyTable()
    yieldTableWt = PrettyTable()
    yieldTableWtLs = PrettyTable()
    yieldTable.add_column('Selections',evtCutFlowLabels)
    yieldTableWt.add_column('Selections',evtCutFlowLabels)
    yieldTableWtLs.add_column('Selections',evtCutFlowLabels)
    for label, histList in evtCutFlowDict.items():
        entries = []
        hist = histList[0]
        for i in range(hist.GetNbinsX()):
            if not label == 'Observed':
                entries.append(round(hist.GetBinContent(hist.FindBin(i)),3))
            else:
                content = 0.0
                for ih in histList:
                    content += ih.GetBinContent(ih.FindBin(i))
                entries.append(content)

        #if not label == 'Observed' :
        #    yieldTable.add_column(label, entries)
        if label == 'Observed' :
            yieldTableWtLs.add_column('Observed',entries)

        #if not label == 'Observed':
        else:
            yieldTable.add_column(label, entries)
            histWt = histList[1]
            histWtLs = histList[2]
            entriesWt = []
            entriesWtLs = []
            for i in range(hist.GetNbinsX()):
                entriesWt.append(round(histWt.GetBinContent(histWt.FindBin(i)),3))
                entriesWtLs.append(round(histWtLs.GetBinContent(histWtLs.FindBin(i)),3))
            yieldTableWt.add_column(label,entriesWt)
            yieldTableWtLs.add_column(label,entriesWtLs)

    with open (os.path.join(histDir, 'YieldTable.txt'), 'w') as yfile:
        yfile.write(' -------------- | For unweighted events | ---------------- \n')
        yfile.write(str(yieldTable)+'\n')
        yfile.write(' -------------- | For weighted events | ---------------- \n')
        yfile.write(str(yieldTableWt)+'\n')
        yfile.write(f' -------------- | For weighted events (at {lumi} pb-1) | ---------------- \n')
        yfile.write(str(yieldTableWtLs)+'\n')
'''
    labels, cutFlowTables, yieldTable = getYields.getYields(evtCutFlowDict, evtCutFlowLabels, xsec, era)
    for i,tab in enumerate(cutFlowTables):
        cutflowTableFile.write(f' =========== Process : {labels[i]} ========== \n')
        cutflowTableFile.write(str(tab)+'\n')
    yieldTableFile.write(str(yieldTable))

if __name__ == "__main__":
    main()


