# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line args: myNanoProdMc2018 -s NANO --mc --eventcontent NANOAODSIM --datatier NANOAODSIM --no_exec --conditions 106X_upgrade2018_realistic_v16_L1v1 --era Run2_2018,run2_nanoAOD_106Xv2 --customise_commands=process.add_(cms.Service("InitRootHandlers", EnableIMT = cms.untracked.bool(False)))
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run2_2018_cff import Run2_2018
from Configuration.Eras.Modifier_run2_nanoAOD_106Xv2_cff import run2_nanoAOD_106Xv2

process = cms.Process("NANO",Run2_2018,run2_nanoAOD_106Xv2)

# import of standard configurations
process.load("Configuration.StandardSequences.Services_cff")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.EventContent.EventContent_cff")
process.load("SimGeneral.MixingModule.mixNoPU_cfi")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("PhysicsTools.NanoAOD.nano_cff")
process.load("Configuration.StandardSequences.EndOfProcess_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

from LLStaus_Run2.Production.arg_config import *
args = get_args()

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(args.maxEvents)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(args.inputFiles),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string("myNanoProdMc2018 nevts:1"),
    name = cms.untracked.string("Applications"),
    version = cms.untracked.string("$Revision: 1.19 $")
)

# Output definition

process.NANOAODSIMoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string("LZMA"),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string("NANOAODSIM"),
        filterName = cms.untracked.string("")
    ),
    fileName = cms.untracked.string(args.outFile),
    outputCommands = process.NANOAODSIMEventContent.outputCommands
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, "106X_upgrade2018_realistic_v16_L1v1", "")

# Path and EndPath definitions
process.nanoAOD_step = cms.Path(process.nanoSequenceMC)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOAODSIMoutput_step = cms.EndPath(process.NANOAODSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.nanoAOD_step,process.endjob_step,process.NANOAODSIMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

# Automatic addition of the customisation function from PhysicsTools.NanoAOD.nano_cff
from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeMC 

#call to customisation function nanoAOD_customizeMC imported from PhysicsTools.NanoAOD.nano_cff
process = nanoAOD_customizeMC(process)

from LLStaus_Run2.Production.customize_nanoaod_eventcontent_cff import *
customize_process_and_associate(process)

# End of customisation functions

# Customisation from command line

process.add_(cms.Service("InitRootHandlers", EnableIMT = cms.untracked.bool(False)))
# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion


# Debug EDM
if (args.debugEDM) :
    
    process.out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string("debugEDM.root")
    )
    
    process.output_step = cms.EndPath(process.out)
    process.schedule.extend([process.output_step])
