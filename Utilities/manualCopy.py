#gfal-ls srm://se01.indiacms.res.in:8446/srm/managerv2?SFN=/dpm/indiacms.res.in/home/cms/store/user/gsaha/NanoAOD2017-UL17Data/MuonEG/
#gfal-copy 

import os
from subprocess import Popen, PIPE
import logging
import argparse
parser = argparse.ArgumentParser(description='Make Jobs and Send')
parser.add_argument('--source', action='store', required=True, type=str, help='Name of the source dir')

args         = parser.parse_args()
pwd          = os.getcwd()
basepath     = 'srm://se01.indiacms.res.in:8446/srm/managerv2?SFN=/dpm/indiacms.res.in'
files_source = os.path.join('home/cms/store/user/gsaha/NanoAOD2017-UL17Data',args.source) # change it as needed
files_target = os.path.join(pwd, files_source.split('/')[-1])
if not os.path.exists(files_target):
    os.mkdir(files_target)

pathtosearch  = os.path.join(basepath, files_source)
gfallsReport  = Popen(['gfal-ls', pathtosearch], stdout=PIPE).communicate()[0]
print(gfallsReport)
dirList = gfallsReport.split('\n')[0:-1]
print(dirList)

for item in dirList:
    finalpath = os.path.join(files_target, item)
    print('targetpath : ', finalpath)
    if not os.path.exists(finalpath):
        os.mkdir(finalpath)
    else:
        #raise RuntimeError('{} : exists !!!'.format(finalpath))
        print('{} : exists !!!'.format(finalpath))
        continue
    semipath   = os.path.join(pathtosearch, item)
    report     = Popen(['gfal-ls', semipath], stdout=PIPE).communicate()[0]
    targetpath = os.path.join(semipath, report).split('\n')[0]
    endpaths   = Popen(['gfal-ls', targetpath], stdout=PIPE).communicate()[0].split('\n')[0:-1]

    for ipath in endpaths:
        ipath  = os.path.join(targetpath, ipath)
        pathfiles = [x for x in Popen(['gfal-ls', ipath], stdout=PIPE).communicate()[0].split('\n')[0:-1] if '.root' in x]
        print('nFiles in {} : {}'.format(ipath, len(pathfiles)))
        for file in pathfiles:
            cpstatus  = Popen(['gfal-copy', os.path.join(ipath, file), finalpath], stdout=PIPE).communicate()[0]
            print(cpstatus)

