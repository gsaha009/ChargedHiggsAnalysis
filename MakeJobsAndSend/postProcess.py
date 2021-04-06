#####################################################
### A script to postprocess the output root files ###
###                 Gourab Saha                   ###
#####################################################
import os
import yaml
import argparse
import shutil
import fileinput
import sys
import stat
from subprocess import Popen, PIPE
import logging
logging.basicConfig(level=logging.DEBUG,format='%(asctime)s - %(levelname)s - %(message)s',datefmt='%m/%d/%Y %H:%M:%S')

def main():
    parser = argparse.ArgumentParser(description='Make Jobs and Send')
    
    parser.add_argument('--configName', action='store', required=True, type=str, help='Name of the config')
    parser.add_argument('--histdir', action='store', required=True, type=str, help='output directory')

    
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

    # mc samples
    logging.info('Doing hadd ===>')
    for dataType, valDict in samplesDict.items():
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
            for key, val in valDict.items():
                logging.info(' Sample : {}'.format(key))
                filePathList = val.get('filedirs')
                xsec         = val.get('xsec') if not isdata else -999.9
                evtWtSum     = val.get('genEvtWtSum') if not isdata else 'null'
                filesPerJob  = int(val.get('filesPerJob'))
                files        = []
                for item in filePathList:
                    logging.info('\t {}'.format(item))
                    files += [os.path.join(item,rfile) for rfile in os.listdir(item) if '.root' in rfile]

                infileListPerJob = [files[i:i+filesPerJob] for i in range(0, len(files), filesPerJob)]
                logging.info('\t nJobs : {}'.format(len(infileListPerJob)))
                nJobs = len(infileListPerJob)
                
                tobehadd     = []
                posthaddfile = os.path.join(histdir, str(key)+'_hist.root')
                haddcmd_     = ['hadd', 'posthaddfile']
                fileabsent   = False
                for i in range(nJobs) :
                    rootfile = os.path.join(histdir, str(key)+'_'+str(i)+'_hist.root')
                    if not os.path.isfile(rootfile):
                        fileabsent = True
                        logging.info('{}  :: N O T  F O U N D !!!'.format(rootfile))
                    tobehadd.append(rootfile)
                    
                haddcmd = haddcmd_ + tobehadd
                process = Popen(haddcmd, stdout=PIPE)
                print process.communicate()[0]

                if not fileabsent:
                    rmcmd = ['rm'] + tobehadd
                    process2 = Popen(rmcmd, stdout=PIPE)
                    print process2.communicate()[0]


if __name__ == "__main__":
    main()
