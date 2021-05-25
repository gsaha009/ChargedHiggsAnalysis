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
import logging
logging.basicConfig(level=logging.DEBUG,format='%(asctime)s - %(levelname)s - %(message)s',datefmt='%m/%d/%Y %H:%M:%S')

def getCommonInfoDict():
    '''
    dict_ = {
        "configuration": {
            "blinded-range-fill-color": "#FDFBFB",
            "blinded-range-fill-style": "4050",
            "eras": ["2017"],
            "experiment": "CMS",
            "extra-label": "Preliminary results",
            "height": '600',
            "luminosity": {
                '2018': "59740.565201546"
            },
            "luminosity-label": "%1$.2f fb^{-1} (13 TeV)",
            "margin-bottom": "0.13",
            "margin-left": "0.15",
            "margin-right": "0.03",
            "margin-top": "0.05",
            "root": "results",
            "show-overflow": "true",
            "width": "800",
            "yields-table-align": "v"
        }
    }
    '''
    dict_ = {
        "blinded-range-fill-color": "#FDFBFB",
        "blinded-range-fill-style": 4050
    }
    
    return dict_

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
    # mc samples
    logging.info('Doing hadd ===>')
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
                '''
                logging.info(' Sample : {}'.format(key))
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
                for i in range(nJobs) :
                    rootfile = os.path.join(histDir, str(key)+'_'+str(i)+'_hist.root')
                    tobehadd.append(rootfile)
                
                if args.dohadd:
                    haddcmd = haddcmd_ + tobehadd
                    process = Popen(haddcmd, stdout=PIPE)
                    print process.communicate()[0]
                    
                    # removing the job root files because
                    # we have the hadded root files now
                    rmcmd = ['rm'] + tobehadd
                    process2 = Popen(rmcmd, stdout=PIPE)
                    print process2.communicate()[0]
                '''
                group = val.get('group')
                if ismc or issignal:
                    xsec    = val.get('xsec') 
                    histFile  = ROOT.TFile.Open(posthaddfile, 'READ')
                    evWtSum_hist = histFile.Get('pt__PassTightId')
                    evWtSum = evWtSum_hist.Integral()
                elif isdata:
                    xsec    = 1.0
                    evWtSum = 1.0
                
    if args.producePlotYaml:
        print(dictToDump)
        json_object = json.dumps(dictToDump, indent = 4)
        plotJson = os.path.join(histDir, 'plot_'+str(era)+'.json')
        with open(plotJson, 'w') as file:
            plotJson.write(json_object)
                '''

if __name__ == "__main__":
    main()
