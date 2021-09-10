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
import getYields
import makeFakeFileSeparate

def realTimeLogging(proc):
    while True:
        output = proc.stdout.readline()
        if proc.poll() is not None:
            break
        if output:
            logging.info(output.strip().decode("utf-8"))
    rc = proc.poll()


def main():
    logging.basicConfig(level   = logging.DEBUG,
                        format  = '%(asctime)s - %(levelname)s - %(message)s',
                        datefmt = '%m/%d/%Y %H:%M:%S')
    logger = logging.getLogger(__name__)

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
    parser.add_argument('--hasSkim', 
                        action='store_true',
                        required=False,default=False,
                        help="Use this if skimmed ntuples need to be h-added")
    parser.add_argument('--norm',
                        action='store_true',
                        required=False,default=False,
                        help="To produce stack plots but normalised to data")
    parser.add_argument('--noFake',
                        action='store_true',
                        required=False,default=False,
                        help="To produce normalised plots")

    args = parser.parse_args()
    
    with open(args.yaml, 'r') as config:
        configDict = yaml.safe_load(config)

    outdir  = configDict.get('outDir')
    histDir = os.path.join(outdir, args.histdir)
    if not os.path.isdir(histDir):
        raise RuntimeError(f'Couldnt find {histDir}')
    dirKey_condor = '_'.join(args.histdir.split('_')[-2:])

    batchOutDir = os.path.join(histDir, 'batch') # NEW
    if not os.path.isdir(batchOutDir):
        raise RuntimeError(f'batch dir not found : {batchOutDir}')

    resultDir = os.path.join(histDir, 'results') # NEW
    if os.path.isdir(resultDir):
        logging.info(f'{resultDir} exits !!!')
    else:
        os.mkdir(resultDir)

    skimDir = os.path.join(histDir, 'skims') # NEW
    if os.path.isdir(skimDir):
        logging.info(f'{skimDir} exits !!!')
    else:
        os.mkdir(skimDir)

    yieldDir = os.path.join(histDir, 'yields') # NEW
    if os.path.isdir(yieldDir):
        logging.info(f'{yieldDir} exits !!!')
    else:
        os.mkdir(yieldDir)

    '''
    logging.basicConfig (filename=os.path.join(histDir,'postprocess.log'),
                         level=logging.INFO,
                         format='%(asctime)s - %(levelname)s - %(message)s',
                         datefmt='%m/%d/%Y %H:%M:%S')
    # set up logging to console 
    console = logging.StreamHandler() 
    console.setLevel(logging.DEBUG) 
    # set a format which is simpler for console use 
    #formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    #console.setFormatter(formatter) 
    # add the handler to the root logger 
    logging.getLogger('').addHandler(console) 
    logger = logging.getLogger(__name__)
    '''
    pwd = os.getcwd()
    logger.info('present working dir : {}'.format(pwd))
    

    keyList = [str(key) for key in configDict.keys()]
    hasskim = args.hasSkim
    nofake  = args.noFake 

    era            = configDict.get('era')
    lumi           = configDict.get('lumi')

    logger.info('era  : {}'.format(era))
    logger.info('lumi : {} pb-1'.format(lumi))

    jobdir = os.path.join(os.getcwd(), 'JobCards_'+str(era)+'_'+dirKey_condor)

    samplesDict       = configDict.get('samplesDict')
    dataTypes         = [str(sample) for sample in samplesDict.keys()]

    if not os.path.isdir(histDir):
        logger.info('{} : output directory doesnt exist !!! '.format(histDir))

    #if args.produceplotyaml:
    #    dictToDump = dict()
    #    dictToDump.update(getCommonInfoDict())
    #    print(dictToDump)
    
    legendPosDict = configDict.get('legendPos')
    resultDict = defaultdict(list)
    xsecDict   = dict()
    legendDict = dict()
    # mc samples
    logger.info('Doing hadd ===>')
    #groupLegendDict = dict()
    #groupLegendInfo = dict()
    evtCutFlowDict   = defaultdict(list)
    evtCutFlowLabels = []
    combinedDictForYaml = dict()
    resubmitList = []
    cutflowTableFile   = open(os.path.join(yieldDir, f'CutFlow_{era}.txt'),'w')
    yieldTableFile     = open(os.path.join(yieldDir, f'Yields_{era}.txt'),'w')

    for dataType, valDict in samplesDict.items():
        '''
        e.g. dataType : MC
        valDict:
          TTToSemiLeptonic:
          filedirs: ['/eos/user/g/gsaha3/Exotic/MC_UL2017/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/']
          genEvtWtSum: 'genEventSumw'
          xsec: 365.52
          filesPerJob: 2
          group: 'TT'
        ...
        
        '''
        ismc     = False
        isdata   = False
        issignal = False
        if dataType == 'MC' : 
            ismc = True
            hasskim = True
        elif dataType == 'DATA':
            isdata = True
            hasskim = False
        elif dataType == 'SIGNAL':
            issignal = True
            hasskim = True
        else:
            raise RuntimeError('Unknown datatype. Please mention MC, DATA or SIGNAL in the main yaml')
            
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
                logger.info(' Sample : {}'.format(key))
                xsec = val.get('xsec')
                filePathList = val.get('filedirs')
                filesPerJob  = int(val.get('filesPerJob'))
                files        = []
                condorDir = os.path.join(jobdir, key+'_condorJobs')
                for item in filePathList:
                    logger.info('\t {}'.format(item))
                    files += [os.path.join(item,rfile) for rfile in os.listdir(item) if '.root' in rfile]

                infileListPerJob = [files[i:i+filesPerJob] for i in range(0, len(files), filesPerJob)]
                logger.info('\t nJobs : {}'.format(len(infileListPerJob)))
                nJobs = len(infileListPerJob)
                
                batchHistDir = os.path.join(batchOutDir, key) # ..........New
                if not os.path.isdir(batchHistDir):
                    raise RuntimeError(f'{batchHistDir} not found !')

                tobehadd         = []
                #posthaddfile     = os.path.join(histDir, str(key)+'_hist.root')
                posthaddfile     = os.path.join(batchHistDir, str(key)+'_hist.root') # NEW
                haddcmd_         = ['hadd', posthaddfile]
                if hasskim:
                    tobehaddskim         = []
                    posthaddfileskim     = os.path.join(skimDir, str(key)+'_skim.root')
                    haddcmdskim_         = ['hadd', posthaddfileskim]

                fileabsent = False                    
                #if args.dohadd:
                haddCond = os.path.exists(posthaddfile) and os.path.exists(posthaddfileskim) if hasskim else os.path.exists(posthaddfile)
                #if os.path.exists(posthaddfile) :
                if haddCond :
                    logger.info(f'hadded files ---> {posthaddfile} & {posthaddfileskim} already exist!') if hasskim else logger.info(f'hadded file ---> {posthaddfile} already exist!')
                else:
                    with alive_bar(title='h-adding', length=60, enrich_print=True) as bar:
                        for i in range(nJobs) :
                            #rootfile  = os.path.join(histDir, str(key)+'_'+str(i)+'_hist.root')
                            rootfile  = os.path.join(batchHistDir, str(key)+'_'+str(i)+'_hist.root') # NEW
                            tfile    = ROOT.TFile(rootfile,"READ")
                            if not os.path.exists(rootfile): 
                                logger.warning(f'{rootfile} >>-x-x-x-> Missing')
                                resubmitList.append(['condor_submit', os.path.join(condorDir,str(key)+'_'+str(i)+'.sub')])
                                fileabsent = True
                            elif tfile.IsZombie():
                                logger.warning(f'{rootfile} is a Zombie! Please produce this file again.')
                                resubmitList.append(['condor_submit', os.path.join(condorDir,str(key)+'_'+str(i)+'.sub')])
                                tfile.Close()
                                fileabsent = True
                            else:
                                sleep(0.03)
                                tobehadd.append(rootfile)
                    
                            if hasskim:
                                #rootfileS  = os.path.join(histDir, str(key)+'_'+str(i)+'_skim.root')
                                rootfileS  = os.path.join(batchHistDir, str(key)+'_'+str(i)+'_skim.root')
                                tfileS    = ROOT.TFile(rootfileS,"READ")
                                if not os.path.exists(rootfileS): 
                                    logger.warning(f'{rootfileS} >>-x-x-x-> Missing')
                                elif tfileS.IsZombie():
                                    logger.warning(f'{rootfileS} is a Zombie! Please produce this file again.')
                                    tfileS.Close()
                                else:
                                    sleep(0.03)
                                    tobehaddskim.append(rootfileS)

                        if len(tobehadd) == 0:
                            logger.error(f'{key} job output root files are not present')
                            sys.exit(f'TERMINATED!!! no post hadd files for {key}. Check the input files and then check the scripts :( ')  
                        elif nJobs == len(tobehadd):
                            haddcmd = haddcmd_ + tobehadd
                            process = Popen(haddcmd, stdout=PIPE)
                            process.communicate()
                            #realTimeLogging(process)
                            '''
                            # removing the job root files because
                            # we have the hadded root files now
                            if len(tobehadd) == nJobs:
                            rmcmd = ['rm'] + tobehadd
                            process2 = Popen(rmcmd, stdout=PIPE)
                            process2.communicate()
                            else:
                            logging.info('There are some missing or zombie files ... Resolve it before removing !')
                            '''
                        else:
                            logger.warning('Skipping hadd hist files because files are found missing')

                        if hasskim:
                            if len(tobehaddskim) == 0:
                                logger.error(f'{key} job output skim root files are not present')
                                sys.exit(f'TERMINATED!!! no post hadd files for {key}. Check the input files and then check the scripts :( ')  
                            elif nJobs == len(tobehaddskim):
                                haddcmdskim = haddcmdskim_ + tobehaddskim
                                #print(haddcmdskim)
                                processskim = Popen(haddcmdskim, stdout=PIPE)
                                processskim.communicate()
                                #realTimeLogging(process)
                                '''
                                # removing the job root files because
                                # we have the hadded root files now
                                if len(tobehadd) == nJobs:
                                rmcmd = ['rm'] + tobehadd
                                process2 = Popen(rmcmd, stdout=PIPE)
                                process2.communicate()
                                else:
                                logging.info('There are some missing or zombie files ... Resolve it before removing !')
                                '''
                            else:
                                logger.warning('Skipping hadd skim files because files are found missing')

                        bar()
                # hadding done .................
                if not fileabsent:
                    # getting all the evtCutFlow histograms
                    #print(dataType , key)
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
                        #print(cutFlowHist.Integral(), cutFlowWtHist.Integral(), cutFlowWtHist_ls.Integral())
                        cutFlowHist.SetDirectory(0)
                        cutFlowWtHist.SetDirectory(0)
                        cutFlowWtHist_ls.SetDirectory(0)

                    elif dataType == 'DATA':
                        evtCutFlowDict[dataType+'_'+key].append(cutFlowHist)
                        #print(cutFlowHist.Integral())
                        cutFlowHist.SetDirectory(0)
                        
                    else:
                        raise RuntimeError('please mention correct datatype : MC or DATA or SIGNAL ...')

                    SR_File, Fake_File = makeFakeFileSeparate.getSRandFakeRootFiles(posthaddfile, resultDir, isSignal=issignal)
                    
                    group = val.get('group')
                    if os.path.exists(SR_File):
                        resultDict[str(group)].append(SR_File)
                        legendDict[str(group)] = val.get('legend')
                        xsecDict[str(SR_File)]=[xsec, dataType]
                    else:
                        logging.warning('hadded file |{}| is absent'.format(SR_File))
                    if not issignal and os.path.exists(Fake_File):
                        resultDict['Fake'].append(Fake_File)
                        legendDict['Fake'] = 'Fake'
                        
    # Re-submission
    if len(resubmitList) > 0:
        logger.info('Resubmit the following jobs')
        for cmd in resubmitList:
            logger.info(' '.join(cmd)+'\n')
        logger.info('After re-running these failed jobs, do proceed again with the post-processing.')
        sys.exit()
    else:
        logger.info('All h-added files are present. So moving further --->')
     
        
    #logger.debug(resultDict)
    yamlObj = yamlMaker(era,lumi,resultDict,xsecDict,legendDict,os.path.join(histDir,'plots.yml'),histDir,args.norm)
    histograms, histTitles = yamlObj.getListOfHistograms()
    #logger.info(f'List of histograms : \n {histograms} \n')
    commonInfoDict = yamlObj.getCommonInfoDict()
    logger.info('commonInfo >>----> plots.yml')
    #print(yaml.dump(commonInfoDict,default_flow_style=False))
    combinedDictForYaml.update(commonInfoDict)

    fileInfoDict = yamlObj.getFileInfoDict(lumi,nofake)
    logger.info('fileInfo >>----> plots.yml')
    combinedDictForYaml.update(fileInfoDict)

    logger.info('legendInfo >>----> plots.yml')
    #logging.info(yaml.dump(groupLegendDict,default_flow_style=False))
    groupLegendDict = yamlObj.getLegendInfoDict(nofake)
    combinedDictForYaml.update(groupLegendDict)

    logger.info('legendPosInfo >>----> plots.yml')
    #logging.info(yaml.dump(legendPosDict,default_flow_style=False))
    combinedDictForYaml.update(legendPosDict)

    plotInfoDict = yamlObj.getHistogramDict(histograms, histTitles)
    logger.info('plotsInfo >>----> plots.yml')
    combinedDictForYaml.update(plotInfoDict)

    #print(combinedDictForYaml)
    #logger.info('preparing plots.json')
    #json_object = json.dumps(combinedDictForYaml, indent = 4)
    #with open(os.path.join(histDir,'plots.json'), 'w') as outfile:
    #    outfile.write(json_object)

    logger.info('preparing plot yml')
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
        proc       = Popen(plotItCmds, stdout=PIPE)
        # Poll proc.stdout to show stdout live
        while True:
            output = proc.stdout.readline()
            if proc.poll() is not None:
                break
            if output:
                logging.info(output.strip().decode("utf-8"))
        rc = proc.poll()


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


