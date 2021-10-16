import os
import sys
import time
from alive_progress import alive_bar
import ROOT
import copy
import logging
logger = logging.getLogger("joblog")

def getSRandFakeRootFiles(postHaddFilePath, resultdir, noFake=False, isSignal=False):
    haddfile = ROOT.TFile (postHaddFilePath,"READ")
    if (haddfile.IsZombie()):
        logger.error(f'{postHaddFilePath} is a zombie')
        raise RuntimeError('Abort... Check your postprocessor')
    else:
        SR_fileName = os.path.join(resultdir, os.path.basename(postHaddFilePath).replace('_hist',''))
        SB_fileName = os.path.join(resultdir, os.path.basename(postHaddFilePath).replace('_hist','_FakeExtrapolation'))

        if os.path.isfile(SR_fileName):
            logger.info("SR file already exists !!!")
        else:
            SR_file     = ROOT.TFile.Open(SR_fileName, 'RECREATE')
            SR_histList = [key.GetName() for key in haddfile.GetListOfKeys() if not any(['_SB_' in key.GetName(),'ObjectSelection' in key.GetName()])]
            with alive_bar(len(SR_histList), title='Producing file with SR histograms | ', length=30, enrich_print=True, bar='circles') as SR_bar:
                for histName in SR_histList:
                    time.sleep(0.03)
                    hist = copy.deepcopy(haddfile.Get(histName))
                    hist.SetName(histName.replace('_SR_',''))
                    SR_file.cd()
                    hist.Write()
                    SR_bar()
            SR_file.Close()

        if os.path.isfile(SB_fileName):
            logger.info("SB file already exists !!!")
        else:
            if not noFake:
                if not isSignal:
                    SB_file  = ROOT.TFile.Open(SB_fileName, 'RECREATE')
                    SB_histList = [key.GetName() for key in haddfile.GetListOfKeys() if not any(['_SR_' in key.GetName(),'ObjectSelection' in key.GetName()])]
                    with alive_bar(len(SB_histList), title='Producing file with FR histograms | ', length=30, enrich_print=True, bar='circles') as SB_bar:
                        for histName in SB_histList:
                            time.sleep(0.03)
                            hist = copy.deepcopy(haddfile.Get(histName))
                            hist.SetName(histName.replace('_SB_',''))
                            SB_file.cd()
                            hist.Write()
                            SB_bar()
                    SB_file.Close()
                else:
                    logger.info('No FakeExtrapolation for Signal Processes')
            else:
                logger.info('[--nofake] ---> Skipping the production of FakeExtrapolation root file')

        haddfile.Close()
        return SR_fileName, SB_fileName
