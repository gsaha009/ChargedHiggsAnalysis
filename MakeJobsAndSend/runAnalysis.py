#################################################
### A script to make the job submissions easy ###
###               Gourab Saha                 ###
#################################################
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

def replaceAll(file,searchExp,replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)

def dumpInJobCard(jfile, era, commonInfoList, evtWtSum, mvaInfoList, histDir, key, xsec, lumi, cutLists, 
                  HLT_SingleMuon, HLT_DoubleMuon, HLT_SingleElectron, HLT_DoubleEG, HLT_MuonEG, 
                  SFInfo, dataset, files, endInfoList, ismc, isdata, issignal):
    if ismc:
        jfile.write('START\n'+'era '+str(era)+'\n'+'dataType mc\n')
    if issignal:
        jfile.write('START\n'+'era '+str(era)+'\n'+'dataType mc#signal\n')
    if isdata:
        jfile.write('START\n'+'era '+str(era)+'\n'+'dataType data\n')
    for item in commonInfoList:
        jfile.write(item+'\n') 
    jfile.write('evtWtSum '+evtWtSum+'\n')
    jfile.write('############ MVA Info ###############'+'\n')
    for item in mvaInfoList:
        jfile.write(item+'\n')
    jfile.write('########### xsec,lumi,hist ###########\n')
    if not isdata:
        jfile.write('lumiWtList xsec='+str(xsec)+' intLumi='+str(lumi)+' nevents=100000'+'\n')
    jfile.write('histFile '+histDir+'/'+str(key)+'_hist.root'+'\n')
    #jfile.write('logFile '+histDir+'/'+str(key)+'_dump.log'+'\n')
    jfile.write('############ Cut lists ###############'+'\n')
    for item in cutLists:
        jfile.write(item+'\n')
    jfile.write('############ HLT lists ###############'+'\n')
    for item in HLT_SingleMuon:
        jfile.write('SingleMuon '+'HLT_'+item+'\n')
    for item in HLT_DoubleMuon:
        jfile.write('DoubleMuon '+'HLT_'+item+'\n')
    for item in HLT_SingleElectron:
        jfile.write('SingleElectron '+'HLT_'+item+'\n')
    for item in HLT_DoubleEG:
        jfile.write('DoubleEG '+'HLT_'+item+'\n')
    for item in HLT_MuonEG:
        jfile.write('MuonEG '+'HLT_'+item+'\n')
    jfile.write('################ SF ###################'+'\n')
    for item in SFInfo:
        jfile.write(item+'\n')
    jfile.write('############ Input Files ##############'+'\n')
    if isdata:
        jfile.write('dataset '+dataset+'\n')
    for item in files:
        jfile.write('inputFile '+item+'\n')
    jfile.write('########################################\n')
    for item in endInfoList:
        jfile.write(item+'\n')
    jfile.write('END')


def main():
    parser = argparse.ArgumentParser(description='Make Jobs and Send')
    
    parser.add_argument('--configName', action='store', required=True, type=str, help='Name of the config')
    parser.add_argument('--suffix', action='store', required=True, type=str, help='set jobdir name')
    parser.add_argument('--send', action='store_true', required=False, help='send jobs to HT-Condor')
    
    args = parser.parse_args()

    pwd = os.getcwd()
    logging.info(' >>>--------------| JobCard Producer |-------------->')
    logging.info('present working dir : {}'.format(pwd))
    
    with open(args.configName, 'r') as config:
        configDict = yaml.safe_load(config)

    keyList = [str(key) for key in configDict.keys()]
    
    era                = configDict.get('era')
    lumi               = configDict.get('lumi')
    tree               = configDict.get('tree')
    commonInfoList     = configDict.get('commonInfo')
    mvaInfoList        = configDict.get('mvaInfo')
    endInfoList        = configDict.get('endInfo')
    cutLists           = configDict.get('cutLists')    
    HLT_DoubleMuon     = configDict.get('HLT_DoubleMuon')
    HLT_DoubleEG       = configDict.get('HLT_DoubleEG')
    HLT_MuonEG         = configDict.get('HLT_MuonEG')
    HLT_SingleMuon     = configDict.get('HLT_SingleMuon')
    HLT_SingleElectron = configDict.get('HLT_SingleElectron')
    SFInfo             = configDict.get('ScaleFactorsInfo')

    logging.info('era  : {}'.format(era))
    logging.info('lumi : {} pb-1'.format(lumi))
    logging.info('tree : {}'.format(tree))
    
    appdir            = configDict.get('appDir')
    jobdir            = configDict.get('jobDir')
    exeToRun          = configDict.get('exeToRun')
    outdir            = configDict.get('outDir')
    samplesDict       = configDict.get('samplesDict')

    jobDir  = os.path.join(pwd,'JobCards_'+str(era)+'_'+args.suffix) 
    if os.path.isdir(jobDir):
        logging.info('{} : job directory overwritten !'.format(jobDir))
    else:
        os.mkdir(jobDir)
        logging.info('{} : job directory created'.format(jobDir))

    histDir = os.path.join(outdir,'ChargedHiggsAnalysis_'+str(era)+'_JobOutput_'+args.suffix)
    if os.path.isdir(histDir):
        raise RuntimeError('Output directory exists. Please change the suffix .')
    else:
        os.mkdir(histDir)
        logging.info('{} : out directory created'.format(histDir))

    if args.send :
        logging.info('{} : list of JobIds'.format(os.path.join(histDir, 'JobIds.txt')))
        JobIdList   = open(os.path.join(histDir, 'JobIds.txt'), 'w')
    
    logging.info('Start making job cards ............... ')
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
            logging.info('Start making job cards for {} samples >>------>'.format(str(dataType)))
            logging.info('{}_Samples : {}'.format(str(dataType), [sample for sample in valDict.keys()]))
            for key, val in valDict.items():
                logging.info(' Sample : {}'.format(key))
                dataset      = str(key.split('_')[0]) if isdata else 'null'
                filePathList = val.get('filedirs')
                xsec         = val.get('xsec') if not isdata else -999.9
                evtWtSum     = val.get('genEvtWtSum') if not isdata else 'null'
                filesPerJob  = int(val.get('filesPerJob'))
                files        = []
                for item in filePathList:
                    logging.info('\t {}'.format(item))
                    files += [os.path.join(item,rfile) for rfile in os.listdir(item) if '.root' in rfile]
                jobFile = str(key)+'.job'
                with open(os.path.join(jobDir,jobFile), 'w') as jfile:
                    dumpInJobCard(jfile, era, commonInfoList, evtWtSum, mvaInfoList, histDir, key, xsec, lumi, cutLists,
                                  HLT_SingleMuon, HLT_DoubleMuon, HLT_SingleElectron, HLT_DoubleEG, HLT_MuonEG,
                                  SFInfo, dataset, files, endInfoList, ismc, isdata, issignal)
                jfile.close()
                # jobfile production finished

                # Now prepare to send jobs to condor
                conDir = os.path.join(jobDir, str(key)+'_condorJobs')
                if not os.path.isdir(conDir):
                    os.mkdir(conDir)
                tmplsubFile = os.path.join(conDir,str(key)+'.sub.tmpl')
                tmplshFile  = os.path.join(conDir,str(key)+'.sh.tmpl')
                shutil.copy(os.path.join('CondorTemplateScripts','SAMPLE.sub.tmpl'), tmplsubFile)
                shutil.copy(os.path.join('CondorTemplateScripts','SAMPLE.sh.tmpl'), tmplshFile)
                inputFile = str(key)+'_infiles.list'
                with open(os.path.join(conDir, inputFile), 'w') as infile:
                    for file in files:
                        infile.write(file+'\n')
                infile.close()
                if not os.path.isdir(os.path.join(conDir,'log')):
                    os.mkdir(os.path.join(conDir,'log'))
                if not os.path.isdir(os.path.join(conDir,'output')):
                    os.mkdir(os.path.join(conDir,'output'))
                if not os.path.isdir(os.path.join(conDir,'error')):
                    os.mkdir(os.path.join(conDir,'error'))
                if not os.path.isdir(os.path.join(conDir,'runlogs')):
                    os.mkdir(os.path.join(conDir,'runlogs'))
                infileListPerJob = [files[i:i+filesPerJob] for i in range(0, len(files), filesPerJob)]
                logging.info('\t nFiles : {} || nJobs : {}'.format(len(files),len(infileListPerJob)))
                for i, filesList in enumerate(infileListPerJob):
                    jobkey = os.path.join(conDir,str(key)+'_'+str(i)+'.job')
                    subkey = os.path.join(conDir,str(key)+'_'+str(i)+'.sub')
                    shkey  = os.path.join(conDir,str(key)+'_'+str(i)+'.sh')
                    with open(jobkey, 'w') as tmpl:
                        tmplKey = key+'_'+str(i)
                        dumpInJobCard(tmpl, era, commonInfoList, evtWtSum, mvaInfoList, histDir, tmplKey, xsec, lumi, cutLists,
                                      HLT_SingleMuon, HLT_DoubleMuon, HLT_SingleElectron, HLT_DoubleEG, HLT_MuonEG,
                                      SFInfo, dataset, filesList, endInfoList, ismc, isdata, issignal)
                    tmpl.close()
                    
                    # Make the condor scripts ready
                    shutil.copy(tmplsubFile, subkey)
                    shutil.copy(tmplshFile, shkey)
                    # edit the sub.tmpl file
                    replaceAll(subkey, 'executable   = sample_index.sh', 'executable   = '+shkey)
                    replaceAll(subkey, 'output       = output/sample_INDEX.$(ClusterId).$(ProcId).out', 
                               'output = '+os.path.join(conDir,'output',str(key))+'_'+str(i)+'.$(ClusterId).$(ProcId).out')
                    replaceAll(subkey, 'error        = error/sample_INDEX.$(ClusterId).$(ProcId).err',  
                               'error  = '+os.path.join(conDir,'error',str(key))+'_'+str(i)+'.$(ClusterId).$(ProcId).err')
                    replaceAll(subkey, 'log          = log/sample_INDEX.$(ClusterId).log', 
                               'log    = '+os.path.join(conDir,'log',str(key))+'_'+str(i)+'.$(ClusterId).log')
                    # edit the sh.tmpl file
                    replaceAll(shkey, 'JOBDIR=NameOfJobDirGivenInYaml', 'JOBDIR='+jobdir)
                    replaceAll(shkey, 'APPDIR=NameOfAppDirGivenInYaml', 'APPDIR='+appdir)
                    replaceAll(shkey, 'cd $JOBDIR/condor_runlog_dir', 'cd '+os.path.join(conDir,'runlogs'))
                    replaceAll(shkey, 'uname -a > ./sample_INDEX.runlog 2>&1', 'uname -a > ./'+str(key)+'_'+str(i)+'.runlog 2>&1')
                    replaceAll(shkey, '$JOBDIR/EXE $JOBDIR/PathToJobFile/sample_index.job >> ./sample_index.runlog 2>&1', 
                               exeToRun+' '+jobkey+' >> ./'+str(key)+'_'+str(i)+'.runlog 2>&1')
                    
                    # All job files, sub and sh files are ready
                    # Now SHOOT !
                    if args.send :
                        condorJobCommandList = ['condor_submit',subkey]
                        st = os.stat(shkey)
                        os.chmod(shkey, st.st_mode | stat.S_IEXEC)
                        process = Popen(condorJobCommandList, stdout=PIPE)
                        report = process.communicate()[0]
                        print report
                        JobIdList.write(report.split('cluster ')[-1])


    if args.send:
        JobIdList.close()


if __name__ == "__main__":
    main()
