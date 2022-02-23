#!/bin/bash

cmsDriver.py myNanoProdMc2018 \
-s NANO \
--mc \
--eventcontent NANOAODSIM \
--datatier NANOAODSIM \
--no_exec \
--conditions 106X_upgrade2018_realistic_v16_L1v1 \
--era Run2_2018,run2_nanoAOD_106Xv2 \
--customise_commands="process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))" \
