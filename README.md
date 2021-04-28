# ChargedHiggs Analysis

This is a C++ and ROOT based analysis framework works with NanoAOD. The old MiniAOD version is [here](https://github.com/subirsarkar/HZZ4lAnalysis).
## Codes
All the headers are in `interface` and `src` contains all the source codes.
 - AnaBase : contains kind of all basic functions. All the Branches used in this analysis are accessed here. Basic JobCard parsing is also done in `readJob`.
 - AnaUtil : Basic utilities e.g. histogram filling, cutFlow and efficiency charts.
 - PhysicsObjects : Classes for different objects to make the flat structure into object structure on the fly.
 - PhysicsObjSelector : Object selection
 - MultiLeptonMVAna : main analysis class, contains SR selections.
 - multiLeptonMVAna : main function.

## cloning the repo:
```
git clone https://github.com/gsaha009/ChargedHiggsAnalysis.git
````
Right now, this package can be divided into 2 parts. 

 - Analysis
 - JobCard production

For the Analysis part, the following instucions should be follwed.
###### Getting CMSSW_11_1_X [`To access ROOT and g++`]
```
cd ChargedHiggsAnalysis
cmsrel CMSSW_11_1_X
cd CMSSW_11_1_X/src; cmsenv; cd -
```
###### Compile
```
make clean -f *MakefileName*
make cling -f *MakefileName*
make -f *MakefileName*
```
If the codes get compiled successfully, an executable will be produced. Now one would need the JobCards to run the executable i.e. the analysis.
Here comes the 2nd part `MakeJobsAndSend`. It contains a python script which automatize Job Card production and HTCondor submission.
There is also an `yaml` with all relevant information to produce the jobs. Be careful with the names of the directories i.e. jobDir, outDir etc.
One can produce the job files and the condor scripts by using the following command :
```
cd MakeJobsAndSend/
python runAnalysis.py --configName analysisBLA.yml --suffix Foo [--send]
```
N.B. `--send` is used to send jobs to condor directly.