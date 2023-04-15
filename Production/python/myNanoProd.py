# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line args: myNanoProdMc2018 -s NANO --mc --eventcontent NANOAOD --datatier NANOAOD --no_exec --conditions 106X_upgrade2018_realistic_v16_L1v1 --era Run2_2018,run2_nanoAOD_106Xv2 --customise_commands=process.add_(cms.Service("InitRootHandlers", EnableIMT = cms.untracked.bool(False)))
import importlib

import FWCore.ParameterSet.Config as cms

from LLStaus_Run2.Production.arg_config import *
args = get_args()

d_procConfig = {
    "Data": {
        "2016": {
            "condition": "auto:run2_data",
            "era": "Run2_2016",
            "eramodifier": "run2_nanoAOD_106Xv2",
        },
        "2017":{
            "condition": "auto:run2_data",
            "era": "Run2_2017",
            "eramodifier": "run2_nanoAOD_106Xv2",
        },
        "2018":{
            "condition": "auto:run2_data",
            "era": "Run2_2018",
            "eramodifier": "run2_nanoAOD_106Xv2",
        },
    },
    
    "MC": {
        # 2016 conditions not checked yet; just a placeholder for now
        "2016": {
            "condition": "auto:phase1_2016_realistic",
            "era": "Run2_2016",
            "eramodifier": "run2_nanoAOD_106Xv2",
        },
        "2017":{
            "condition": "auto:phase1_2017_realistic",
            "era": "Run2_2017",
            "eramodifier": "run2_nanoAOD_106Xv2",
        },
        "2018":{
            "condition": "auto:phase1_2018_realistic",
            "era": "Run2_2018",
            "eramodifier": "run2_nanoAOD_106Xv2",
        },
    }
}

isMC = (args.sampleType == "MC")
condition_str = d_procConfig[args.sampleType][args.era]["condition"]
era_str = d_procConfig[args.sampleType][args.era]["era"]
eramodifier_str = d_procConfig[args.sampleType][args.era]["eramodifier"]

era_cff = importlib.import_module(f"Configuration.Eras.Era_{era_str}_cff")
era = getattr(era_cff, era_str)

eramodifier_cff = importlib.import_module(f"Configuration.Eras.Modifier_{eramodifier_str}_cff")
eramodifier = getattr(eramodifier_cff, eramodifier_str)

process = cms.Process("NANO", era, eramodifier)

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

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, condition_str, "")

process.MessageLogger.cerr.enableStatistics = True

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(args.maxEvents)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(),
    secondaryFileNames = cms.untracked.vstring()
    )
from LLStaus_Run2.Production.readFileList import *
if len(args.inputFiles) > 0:
    addList(process.source.fileNames, args.inputFiles, fileNamePrefix=args.fileNamePrefix)
elif len(args.sourceFile) > 0:
    addList(process.source.fileNames, args.inputFiles, fileNamePrefix=None)

if len(args.lumiFile) > 0:
    import FWCore.PythonUtilities.LumiList as LumiList
    process.source.lumisToProcess = LumiList.LumiList(filename = args.lumiFile).getVLuminosityBlockRange()

if args.eventRange != "":
    process.source.eventsToProcess = cms.untracked.VEventRange(re.split(",", args.eventRange))

if args.maxEvents > 0:
    process.maxEvents.input = args.maxEvents

process.options = cms.untracked.PSet()

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string("myNanoProd{args.sampleType}{args.era}"),
    name = cms.untracked.string("Applications"),
    version = cms.untracked.string("$Revision: 1.19 $")
)

# Output definition
assert(args.disTauTagOutputOpt in [0, 1, 2])

if isMC :
    outputCommands = process.NANOAODSIMEventContent.outputCommands
else :
    outputCommands = process.NANOAODEventContent.outputCommands

if args.disTauTagOutputOpt == 1 :
    
    args.outFile = args.outFile.replace(".root", "_with-disTauTagScore.root")

elif args.disTauTagOutputOpt == 2 :
    
    outputCommands = cms.untracked.vstring(
        "drop *",
        #"keep *_*_*disTauTag*_*",
        #"keep nanoaodFlatTable_*Table_*_*",
        #"keep nanoaodFlatTable_*Table_*_*",
        "keep nanoaodFlatTable_jetTable_*_*",
    )
    
    args.outFile = args.outFile.replace(".root", "_only-disTauTagScore.root")

process.NANOAODoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string("LZMA"),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string("NANOAODSIM") if isMC else cms.untracked.string("NANOAOD"),
        filterName = cms.untracked.string("")
    ),
    fileName = cms.untracked.string(args.outFile),
    outputCommands = outputCommands,
)

# Additional output definition

# Other statements

# Path and EndPath definitions
if isMC :
    process.nanoAOD_step = cms.Path(process.nanoSequenceMC)
else :
    process.nanoAOD_step = cms.Path(process.nanoSequence)

process.endjob_step = cms.EndPath(process.endOfProcess)
process.NANOAODoutput_step = cms.EndPath(process.NANOAODoutput)

# Schedule definition
process.schedule = cms.Schedule(process.nanoAOD_step,process.endjob_step,process.NANOAODoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# customisation of the process.

if isMC :
    from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeMC 
    process = nanoAOD_customizeMC(process)

else :
    from PhysicsTools.NanoAOD.nano_cff import nanoAOD_customizeData
    process = nanoAOD_customizeData(process)

from LLStaus_Run2.Production.customize_nanoaod_eventcontent_cff import *
customize_process_and_associate(process, isMC = isMC, disTauTagOutputOpt = args.disTauTagOutputOpt)

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
