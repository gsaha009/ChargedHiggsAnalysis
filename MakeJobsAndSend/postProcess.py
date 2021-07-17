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
import logging
logging.basicConfig(level=logging.DEBUG,format='%(asctime)s - %(levelname)s - %(message)s',datefmt='%m/%d/%Y %H:%M:%S')
from plotter import plotter

def main():
    parser = argparse.ArgumentParser(description='Make Jobs and Send')
    
    parser.add_argument('--configName', action='store', required=True, type=str, help='Name of the config')
    parser.add_argument('--histdir', action='store', required=True, type=str, help='just name of the output directory')
    parser.add_argument('--dohadd', action='store_true',required=False,default=False,help="")
    parser.add_argument('--producePlotYaml',action='store_true',required=False,default=False,help="")

    args = parser.parse_args()

    pwd = os.getcwd()
    logging.info('present working dir : {}'.format(pwd))
    
    with open(args.configName, 'r') as config:
        configDict = yaml.safe_load(config)

    keyList = [str(key) for key in configDict.keys()]
    
    era            = configDict.get('era')
    lumi           = configDict.get('lumi')

    logging.info('era  : {}'.format(era))
    logging.info('lumi : {} pb-1'.format(lumi))
    
    outdir            = configDict.get('outDir')
    samplesDict       = configDict.get('samplesDict')
    dataTypes         = [str(sample) for sample in samplesDict.keys()]

    histDir = os.path.join(outdir, args.histdir)

    if not os.path.isdir(histDir):
        logging.info('{} : output directory doesnt exist !!! '.format(histDir))

    if args.producePlotYaml:
        dictToDump = dict()
        dictToDump.update(getCommonInfoDict())
        print(dictToDump)
    
    legendPosDict = configDict.get('legendPos')
    resultDict = defaultdict(list)
    xsecDict   = dict()
    # mc samples
    logging.info('Doing hadd ===>')
    groupLegendDict = dict()
    groupLegendInfo = dict()
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
            raise RuntimeError('Unknown datatype. Please mention MC, DATA or SIGNAL')

        if valDict == None or len(valDict) == 0:
            logging.info('Sorry! No {} samples are present in yaml'.format(str(dataType)))
        else:
            logging.info('Looking for {} - output root files >>------>'.format(str(dataType)))
            logging.info('{}_Samples : {}'.format(str(dataType), [sample for sample in valDict.keys()]))
            # Loop over valDict
            for key, val in valDict.items():
                '''
                key == sample name e.g. TTToSemiLeptonic
                val == a dictionary with keys and values 
                filedirs: ['/eos/user/g/gsaha3/Exotic/MC_UL2017/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/']
                genEvtWtSum: 'genEventSumw'
                xsec: 365.52
                filesPerJob: 2
                group: 'TT'
                fill-color: '#1a83a1'
                legend: 'tt+jets'
                order: 2
                '''
                logging.info(' Sample : {}'.format(key))
                xsec = val.get('xsec')
                filePathList = val.get('filedirs')
                filesPerJob  = int(val.get('filesPerJob'))
                files        = []
                for item in filePathList:
                    logging.info('\t {}'.format(item))
                    files += [os.path.join(item,rfile) for rfile in os.listdir(item) if '.root' in rfile]

                infileListPerJob = [files[i:i+filesPerJob] for i in range(0, len(files), filesPerJob)]
                logging.info('\t nJobs : {}'.format(len(infileListPerJob)))
                nJobs = len(infileListPerJob)
                
                tobehadd     = []
                posthaddfile = os.path.join(histDir, str(key)+'_hist.root')
                haddcmd_     = ['hadd', posthaddfile]

                if args.dohadd:
                    if os.path.exists(posthaddfile):
                        print('hadded file ---> {} : already exists!'.format(posthaddfile))
                    else:
                        for i in range(nJobs) :
                            rootfile = os.path.join(histDir, str(key)+'_'+str(i)+'_hist.root')
                            if os.path.exists(rootfile):
                                tobehadd.append(rootfile)   
                            else:
                                print('{} >>-----> Missing'.format(rootfile))

                        if len(tobehadd) == 0:
                            print('{} job output root files are not present'.format(key))
                        else:
                            haddcmd = haddcmd_ + tobehadd
                            process = Popen(haddcmd, stdout=PIPE)
                            print(process.communicate()[0])                    
                            # removing the job root files because
                            # we have the hadded root files now
                            rmcmd = ['rm'] + tobehadd
                            process2 = Popen(rmcmd, stdout=PIPE)
                            print(process2.communicate()[0])

                group = val.get('group')
                if os.path.exists(posthaddfile):
                    resultDict[str(group)].append(posthaddfile)
                    xsecDict[str(posthaddfile)]=[xsec, dataType]
                else:
                    print('hadded file |{}| is absent'.format(posthaddfile))
                if not group == 'data':
                    groupLegendInfo[group] = {
                        "fill-color": val.get('fill-color'),
                        "legend": val.get('legend'),
                        "order": val.get('order')
                    }
                else:
                    groupLegendInfo[group] = {
                        "legend":'data'
                    }
                groupLegendDict['groups'] = groupLegendInfo
                    
    print(resultDict)
    plotterObj = plotter(era,lumi,resultDict,xsecDict,os.path.join(histDir,'plots.yml'))
    histograms = plotterObj.getListOfHistograms()
    print(histograms)
    commonInfoDict = plotterObj.getCommonInfoDict()
    print('commonInfo :: to be dumped inside plots.yml')
    print(commonInfoDict)
    print(type(commonInfoDict))
    combinedDictForYaml.update(commonInfoDict)
    fileInfoDict = plotterObj.getFileInfoDict()
    print('fileInfo :: to be dumped inside plots.yml')
    print(fileInfoDict)
    print(type(fileInfoDict))
    combinedDictForYaml.update(fileInfoDict)
    print('legendInfo :: to be dumped inside plots.yml')
    print(groupLegendDict)
    combinedDictForYaml.update(groupLegendDict)
    print('legendPosInfo :: to be dumped inside plots.yml')
    print(legendPosDict)
    combinedDictForYaml.update(legendPosDict)

    print(combinedDictForYaml)
    json_object = json.dumps(combinedDictForYaml, indent = 4)
    with open(os.path.join(histDir,'plots.json'), 'w') as outfile:
        outfile.write(json_object)

    with open(os.path.join(histDir, 'plots.yml'), 'w') as ymlfile:
        documents = yaml.dump(combinedDictForYaml, ymlfile, default_flow_style=False)

if __name__ == "__main__":
    main()
