import os
from subprocess import Popen, PIPE
import logging

basepath  = '/eos/user/g/gsaha3/Exotic/DATA_UL2017'
workdirs  = ['DoubleMuon/RunB','DoubleMuon/RunC','DoubleMuon/RunD','DoubleMuon/RunE','DoubleMuon/RunF']
nMerge    = 100 

logging.basicConfig(filename=os.path.join(basepath,workdirs[0].split('/')[0]+'_hadd.log'), filemode='w',level=logging.DEBUG,format='%(asctime)s - %(levelname)s - %(message)s',datefmt='%m/%d/%Y %H:%M:%S')

for workdir in workdirs:
    logging.info('hadd the root files in >>>------> {}'.format(workdir))
    infiles  = [os.path.join(basepath,workdir,file) for file in os.listdir(os.path.join(basepath,workdir)) if 'root' in file]
    for i in range(1000):
        check = os.path.join(basepath,workdir,'tree_'+str(i+1)+'.root')
        if not os.path.isfile(check):
            logging.warning('{} missing '.format(check))

    nFiles   = int(len(infiles)/nMerge)
    if len(infiles)%nMerge > 0 :
        nFiles = nFiles+1

    logging.info('{} files are going to be produced'.format(nFiles)) 
    slicedinfiles = [infiles[i*nMerge:i*nMerge+nMerge] for i in range(nFiles)]
    for i, filestoadd in enumerate(slicedinfiles):
        outfile    = os.path.join(basepath,workdir, workdir.replace('/','_')+'_'+str(i+1)+'.root')
        cmdlist    = ['hadd',outfile]
        for file in filestoadd:
            cmdlist.append(file)
        logging.info('Merging {} files and producing {}'.format(len(filestoadd), outfile))
        logging.info('hadding : >>>----->')
        for item in filestoadd:
            logging.info('--- {}'.format(item))
        #logging.info('nFiles : {}'.format(len(filestoadd)))
        logging.info('Starting hadding >>>------> {}'.format(cmdlist))
        process    = Popen(cmdlist, stdout=PIPE)                                                                                  
        report     = process.communicate()[0]                                                                                                     
        logging.info(report)    
