## ChargedHiggs Analysis

This is a C++ based analysis framework. Basically it uses C++14 and ROOT. For the time being, a few stuffs like `unique_pointer` etc need CMSSW_9_4_9 [ToDo : check the CMSSW dependency]. So before running the package inside `Analysis`, one should get the CMSSW. \n

 - cmsrel CMSSW_9_4_9
 - cd CMSSW_9_4_9/src; cmsenv; cd -

`MakeJobsAndSend` contains a python script which automatize Job Card production and HTCondor submission.