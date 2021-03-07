# ChargedHiggs Analysis

This is a C++ based analysis framework. Here ROOT is used for plotting and MVA stuff.\
###### cloning the repo:
```
git clone https://github.com/gsaha009/ChargedHiggsAnalysis.git
````
Right now, this package can be divided into 2 parts. 

 - Analysis
 - JobCard production

For the Analysis part, the following instucions should be follwed.
###### Getting CMSSW_9_4_9
```
cd ChargedHiggsAnalysis
cmsrel CMSSW_9_4_9
cd CMSSW_9_4_9/src; cmsenv; cd -
```
###### Compile
```
make clean -f *MakefileName*
make cling -f *MakefileName*
make -f *MakefileName*
```
If the codes get compiled successfully, an executable will be produced. Now one would need the JobCards to run the executable i.e. the analysis.
Here comes the 2nd part `MakeJobsAndSend`. It contains a python script which automatize Job Card production and HTCondor submission.\
There is also an `yaml` with all relevant information to produce the jobs.