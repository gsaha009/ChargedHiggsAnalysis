import os
import sys
import time
from alive_progress import alive_bar
import ROOT
import copy
import logging
logger = logging.getLogger(__name__)

def getSRandFakeRootFiles(postHaddFilePath, isSignal=False):
    haddfile = ROOT.TFile (postHaddFilePath,"READ")
    if (haddfile.IsZombie()):
        logger.warning(f'{postHaddFilePath} is a zombie')
    else:
        SR_fileName = postHaddFilePath.replace('_hist','')
        SB_fileName = postHaddFilePath.replace('_hist','_FakeExtrapolation')
        if os.path.isfile(SR_fileName) and os.path.isfile(SB_fileName):
            logger.info("files already exist !!!")
        else:
            SR_file  = ROOT.TFile.Open(SR_fileName, 'RECREATE')
            SR_histList = [key.GetName() for key in haddfile.GetListOfKeys() if not any(['_SB_' in key.GetName(),'ObjectSelection' in key.GetName()])]
            SB_file  = ROOT.TFile.Open(SB_fileName, 'RECREATE')
            if not isSignal:
                SB_histList = [key.GetName() for key in haddfile.GetListOfKeys() if not any(['_SR_' in key.GetName(),'ObjectSelection' in key.GetName()])]
                
                if not len(SR_histList) == len(SB_histList):
                    raise RuntimeError('no of histograms is not same for SR and FakeRegion !!!')
                    
            logger.info('Producing root file for Signal & FakeExtrapolation regions ...')

            with alive_bar(len(SR_histList), title='Producing file with SR histograms | ', length=30, enrich_print=True, bar='circles') as SR_bar:
                for histName in SR_histList:
                    time.sleep(0.03)
                    hist = copy.deepcopy(haddfile.Get(histName))
                    hist.SetName(histName.replace('_SR_',''))
                    SR_file.cd()
                    hist.Write()
                    SR_bar()
            SR_file.Close()

            if not isSignal:
                with alive_bar(len(SR_histList), title='Producing file with FR histograms | ', length=30, enrich_print=True, bar='circles') as SB_bar:
                    for histName in SB_histList:
                        time.sleep(0.03)
                        hist = copy.deepcopy(haddfile.Get(histName))
                        hist.SetName(histName.replace('_SB_',''))
                        SB_file.cd()
                        hist.Write()
                        SB_bar()
            SB_file.Close()

    haddfile.Close()
    return SR_fileName, SB_fileName
