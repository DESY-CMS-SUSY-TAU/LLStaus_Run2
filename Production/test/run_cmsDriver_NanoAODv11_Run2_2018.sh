#!/bin/bash

cmsDriver.py myNanoProdData2018 \
-s NANO \
--data \
--eventcontent NANOAOD \
--datatier NANOAOD \
--no_exec \
--conditions "auto:phase1_2018_realistic" \
--era Run2_2018,run2_nanoAOD_106Xv2 \
--customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" \
--customise JMEAnalysis/JetToolbox/nanoAOD_jetToolbox_cff.nanoJTB_customizeMC
