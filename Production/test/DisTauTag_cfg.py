# coding: utf-8

import os

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing


# get the data/ directory
graph_file = 'LLStaus_Run2/Production/data/models/particlenet_v1_a27159734e304ea4b7f9e0042baa9e22.pb'
graph_file_abs = os.path.abspath(graph_file)

if not os.path.exists(graph_file_abs):
    raise Exception("Error: graph_file is not found")

# setup minimal options
options = VarParsing("python")
options.setDefault("inputFiles", "root://xrootd-cms.infn.it//store/user/myshched/mc/UL2018-pythia-v5-100cm/SUS-RunIISummer20UL18GEN-stau100_lsp1_ctau100mm_v5/MiniAOD/220525_205521/0000/SUS-RunIISummer20UL18MiniAODv2-LLStau_401.root")  # noqa
options.parseArguments()

# define the process to run
process = cms.Process("TEST")

# minimal configuration
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(1000))
process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring(options.inputFiles))

# process options
process.options = cms.untracked.PSet(
    allowUnscheduled=cms.untracked.bool(True),
    wantSummary=cms.untracked.bool(True),
)

# setup DisTauTag by loading the auto-generated cfi (see DisTauTag.fillDescriptions)
# process.load("LLStaus_Run2.Production.disTauTag_cfi")
# process.disTauTag.graphPath    = cms.string(graph_file_abs)
# process.disTauTag.jets         = cms.InputTag('slimmedJets')
# process.disTauTag.pfCandidates = cms.InputTag('packedPFCandidates')

process.disTauTag = cms.EDProducer("DisTauTag",
    graphPath = cms.string(graph_file_abs),
    jets         = cms.InputTag('slimmedJets'),
    pfCandidates = cms.InputTag('packedPFCandidates'),
    save_inputs  = cms.bool(True)
)

# define what to run in the path
process.p = cms.Path(process.disTauTag)