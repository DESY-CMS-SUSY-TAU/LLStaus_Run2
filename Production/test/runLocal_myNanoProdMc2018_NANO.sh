mkdir -p logs


#nohup cmsRun LLStaus_Run2/Production/python/myNanoProdMc2018_NANO.py \
#sourceFile=LLStaus_Run2/Production/configs/sourceFiles/SUS-RunIISummer20UL18GEN-stau100_lsp1_ctau1000mm_v4_myshched-MiniAOD-e2b541b9a2622f5ad47ca02585af3f7c_USER/SUS-RunIISummer20UL18GEN-stau100_lsp1_ctau1000mm_v4_myshched-MiniAOD-e2b541b9a2622f5ad47ca02585af3f7c_USER.txt \
#disTauTagOutputOpt=1 \
#maxEvents=30000 \
#outFile=output/nanoaod_stau100_lsp1_ctau1000mm.root \
#> logs/myNanoProdMc2018_NANO_stau100_lsp1_ctau1000mm.log &
#
#
#
#nohup cmsRun LLStaus_Run2/Production/python/myNanoProdMc2018_NANO.py \
#sourceFile=LLStaus_Run2/Production/configs/sourceFiles/SUS-RunIISummer20UL18GEN-stau250_lsp1_ctau1000mm_v4_myshched-MiniAOD-e2b541b9a2622f5ad47ca02585af3f7c_USER/SUS-RunIISummer20UL18GEN-stau250_lsp1_ctau1000mm_v4_myshched-MiniAOD-e2b541b9a2622f5ad47ca02585af3f7c_USER.txt \
#disTauTagOutputOpt=1 \
#maxEvents=30000 \
#outFile=output/nanoaod_stau250_lsp1_ctau1000mm.root \
#> logs/myNanoProdMc2018_NANO_stau250_lsp1_ctau1000mm.log &
#
#
#
#nohup cmsRun LLStaus_Run2/Production/python/myNanoProdMc2018_NANO.py \
#sourceFile=LLStaus_Run2/Production/configs/sourceFiles/SUS-RunIISummer20UL18GEN-stau400_lsp1_ctau1000mm_v4_myshched-MiniAOD-e2b541b9a2622f5ad47ca02585af3f7c_USER/SUS-RunIISummer20UL18GEN-stau400_lsp1_ctau1000mm_v4_myshched-MiniAOD-e2b541b9a2622f5ad47ca02585af3f7c_USER.txt \
#disTauTagOutputOpt=1 \
#maxEvents=30000 \
#outFile=output/nanoaod_stau400_lsp1_ctau1000mm.root \
#> logs/myNanoProdMc2018_NANO_stau400_lsp1_ctau1000mm.log &
#
#
#
#nohup cmsRun LLStaus_Run2/Production/python/myNanoProdMc2018_NANO.py \
#sourceFile=LLStaus_Run2/Production/configs/sourceFiles/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1_MINIAODSIM/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1_MINIAODSIM.txt \
#disTauTagOutputOpt=1 \
#maxEvents=30000 \
#outFile=output/nanoaod_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root \
#> logs/myNanoProdMc2018_TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.log &



nohup cmsRun LLStaus_Run2/Production/python/myNanoProdMc2018_NANO.py \
sourceFile=LLStaus_Run2/Production/configs/sourceFiles/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2_MINIAODSIM/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2_MINIAODSIM.txt \
disTauTagOutputOpt=1 \
maxEvents=30000 \
outFile=output/nanoaod_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root \
> logs/myNanoProdMc2018_TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.log &



nohup cmsRun LLStaus_Run2/Production/python/myNanoProdMc2018_NANO.py \
sourceFile=LLStaus_Run2/Production/configs/sourceFiles/TTToHadronic_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1_MINIAODSIM/TTToHadronic_TuneCP5_13TeV-powheg-pythia8_RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v1_MINIAODSIM.txt \
disTauTagOutputOpt=1 \
maxEvents=30000 \
outFile=output/nanoaod_TTToHadronic_TuneCP5_13TeV-powheg-pythia8.root \
> logs/myNanoProdMc2018_TTToHadronic_TuneCP5_13TeV-powheg-pythia8.log &

