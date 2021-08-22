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
from subprocess import Popen, PIPE
from collections import defaultdict
from alive_progress import alive_bar
from time import sleep
import logging
from yamlMaker import yamlMaker
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
    parser.add_argument('--configname', action='store', required=True, type=str, help='Name of the config')
    parser.add_argument('--histdir', action='store', required=True, type=str, help='just name of the output directory')
    parser.add_argument('--runplotit',action='store_true',required=False,default=False,help="")
    parser.add_argument('--norm',action='store_true',required=False,default=False,help="To produce normalised plots")

    args = parser.parse_args()
    
    with open(args.configname, 'r') as config:
        configDict = yaml.safe_load(config)

    outdir  = configDict.get('outDir')
    histDir = os.path.join(outdir, args.histdir)
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
    
    era            = configDict.get('era')
    lumi           = configDict.get('lumi')

    logger.info('era  : {}'.format(era))
    logger.info('lumi : {} pb-1'.format(lumi))
    
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
    combinedDictForYaml = dict()
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
        elif dataType == 'DATA':
            isdata = True
        elif dataType == 'SIGNAL':
            issignal = True
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
                for item in filePathList:
                    logger.info('\t {}'.format(item))
                    files += [os.path.join(item,rfile) for rfile in os.listdir(item) if '.root' in rfile]

                infileListPerJob = [files[i:i+filesPerJob] for i in range(0, len(files), filesPerJob)]
                logger.info('\t nJobs : {}'.format(len(infileListPerJob)))
                nJobs = len(infileListPerJob)
                
                tobehadd         = []
                posthaddfile     = os.path.join(histDir, str(key)+'_hist.root')
                haddcmd_         = ['hadd', posthaddfile]

                #if args.dohadd:
                if os.path.exists(posthaddfile):
                    logger.info('hadded file ---> {} : already exists!'.format(posthaddfile))
                else:
                    with alive_bar(title='h-adding', length=60, enrich_print=True) as bar:
                        for i in range(nJobs) :
                            rootfile  = os.path.join(histDir, str(key)+'_'+str(i)+'_hist.root')
                            tfile    = ROOT.TFile(rootfile,"READ")
                            if not os.path.exists(rootfile): 
                                logger.warning(f'{rootfile} >>-x-x-x-> Missing')
                            elif tfile.IsZombie():
                                logger.warning(f'{rootfile} is a Zombie! Please produce this file again.')
                                tfile.Close()
                            else:
                                #logging.info(f'to be h-added >>-----> : {rootfile}')
                                sleep(0.03)
                                tobehadd.append(rootfile)
                            #bar()
                    
                        if len(tobehadd) == 0:
                            logger.error(f'{key} job output root files are not present')
                            sys.exit(f'TERMINATED!!! no post hadd files for {key}. Check the input files and then check the scripts :( ')  
                        else:
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
                        bar()
                # hadding done .................
                SR_File, Fake_File = makeFakeFileSeparate.getSRandFakeRootFiles(posthaddfile, isSignal=issignal)
    
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
                    
    logger.debug(resultDict)
    yamlObj = yamlMaker(era,lumi,resultDict,xsecDict,legendDict,os.path.join(histDir,'plots.yml'),histDir,args.norm)
    histograms, histTitles = yamlObj.getListOfHistograms()
    logger.info(f'List of histograms : \n {histograms} \n')
    commonInfoDict = yamlObj.getCommonInfoDict()
    logger.info('commonInfo >>----> plots.yml')
    #print(yaml.dump(commonInfoDict,default_flow_style=False))
    combinedDictForYaml.update(commonInfoDict)

    fileInfoDict = yamlObj.getFileInfoDict()
    logger.info('fileInfo >>----> plots.yml')
    combinedDictForYaml.update(fileInfoDict)

    logger.info('legendInfo >>----> plots.yml')
    #logging.info(yaml.dump(groupLegendDict,default_flow_style=False))
    groupLegendDict = yamlObj.getLegendInfoDict()
    combinedDictForYaml.update(groupLegendDict)

    logger.info('legendPosInfo >>----> plots.yml')
    #logging.info(yaml.dump(legendPosDict,default_flow_style=False))
    combinedDictForYaml.update(legendPosDict)

    plotInfoDict = yamlObj.getHistogramDict(histograms, histTitles)
    logger.info('plotsInfo >>----> plots.yml')
    combinedDictForYaml.update(plotInfoDict)

    #print(combinedDictForYaml)
    logger.info('preparing plots.json')
    json_object = json.dumps(combinedDictForYaml, indent = 4)
    with open(os.path.join(histDir,'plots.json'), 'w') as outfile:
        outfile.write(json_object)

    logger.info('preparing plots.yml')
    with open(os.path.join(histDir, 'plots.yml'), 'w') as ymlfile:
        documents = yaml.dump(combinedDictForYaml, ymlfile, default_flow_style=False)

    if args.runplotit:
        logger.info('running PlotIt ... \n')
        plotDir = os.path.join(histDir, 'plots') if not args.norm else os.path.join(histDir, 'norm_plots')
        if not os.path.exists(plotDir):
            os.mkdir(plotDir)
        
        if not os.path.exists(os.path.join(histDir, 'plots.yml')):
            raise RuntimeError('plots.yml doesnt exist !!!')
        # run plotIt
        plotItCmds = ['plotIt','-y','-o',plotDir,os.path.join(histDir, 'plots.yml')]
        proc       = Popen(plotItCmds, stdout=PIPE)
        # Poll proc.stdout to show stdout live
        while True:
            output = proc.stdout.readline()
            if proc.poll() is not None:
                break
            if output:
                logging.info(output.strip().decode("utf-8"))
        rc = proc.poll()

if __name__ == "__main__":
    main()


