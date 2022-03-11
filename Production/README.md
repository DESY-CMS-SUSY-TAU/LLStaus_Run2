# NANOAOD production



## NanoAOD production

0. To source the environment:
```sh
> cd $CMSSW_BASE/src
> cmsenv
> source /cvmfs/cms.cern.ch/common/crab-setup.sh
```

1. To test the cms config file for NanoAOD production (for 1000 events of the signal SUSY MC):
```sh
> cmsRun LLStaus_Run2/Production/python/myNanoProdMc2018_NANO.py inputFiles=/store/user/myshched/mc/UL2018-pythia-v4/SUS-RunIISummer20UL18GEN-stau250_lsp1_ctau1000mm_v4/MiniAOD/220129_215847/0001/SUS-RunIISummer20UL18MiniAODv2-LLStau_1251.root fileNamePrefix=root://cms-xrd-global.cern.ch/ maxEvents=1000
```

2. To enable VOMS proxy:
```sh
> voms-proxy-init -rfc -voms cms -valid 192:00
> export X509_USER_PROXY=`voms-proxy-info -path`
```

3. To submit all datasets in configuration files use crab_submit.py:
```sh
> crab_submit.py --inputDBS phys03 --workArea <working/area/folder> --cfg LLStaus_Run2/Production/python/myNanoProdMc2018_NANO.py --site T2_DE_DESY --output <path/to/dcache/folder> ./LLStaus_Run2/Production/configs/crab/UL2018/STauSignal.txt .... <other/confo>
```
- For more command line options use crab_submit.py --help.
- For big dataset file-based splitting should be used.
- For private production `phys03` DBS should be specified, for official production nothing is needed.

4. Regularly check task status using crab_cmd.py:
```sh
crab_cmd.py --workArea work-area --cmd status
```

5. If some jobs are failed: try to understand the reason and use standard crab tools to solve the problem (e.g. crab resubmit with additional arguments). In very problematic cases a recovery task could be created.
