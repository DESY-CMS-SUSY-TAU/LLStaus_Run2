import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *

def customize_process_and_associate(process) :
    
    # Lost tracks
    process.lostTrackTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src = cms.InputTag("lostTracks"),
        #cut = cms.string(""),
        cut = cms.string("pt > 1"),
        name= cms.string("LostTrack"),
        doc = cms.string("Lost tracks"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table
        variables = cms.PSet(
            CandVars,
        )
    )
    
    
    # PF candidates
    process.isFromTauForPfCand = cms.EDProducer("IsFromPatTauMapProducer",
        packedPFCandidates = cms.InputTag("packedPFCandidates"),
        #patTaus = cms.InputTag("slimmedTaus"),
        patTaus = cms.InputTag("linkedObjects", "taus"),
        #patTaus = cms.InputTag("selectedPatTaus"),
    )
    
    process.pfCandTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src = cms.InputTag("packedPFCandidates"),
        #cut = cms.string(""),
        cut = cms.string("pt > 1"),
        name= cms.string("PFCandidate"),
        doc = cms.string("PF candidates"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table
        variables = cms.PSet(
            CandVars,
        ),
        externalVariables = cms.PSet(
            isTauSignalCand     = ExtVar("isFromTauForPfCand:isTauSignalCand"       , int, doc = "Belongs to pat::Tau::signalCands()"),
            isTauIsoCand        = ExtVar("isFromTauForPfCand:isTauIsoCand"          , int, doc = "Belongs to pat::Tau::isolationCands()"),
            isTauLeadChHadCand  = ExtVar("isFromTauForPfCand:isTauLeadChHadCand"    , int, doc = "Is pat::Tau::leadChargedHadrCand()"),
        )
    )
    
    
    # Unfiltered taus
    from PhysicsTools.NanoAOD.taus_cff import *
    
    myfinalTaus = finalTaus.clone(
        src = cms.InputTag("slimmedTausUpdated"),
        cut = cms.string("pt > 18"),
    )
    
    process.globalReplace("finalTaus", myfinalTaus)
    
    
    # Create the task
    process.custom_nanoaod_task = cms.Task(
        myfinalTaus,
        
        #process.lostTrackTable,
        
        #process.isFromTauForPfCand,
        #process.pfCandTable,
    )
    
    # Associate the task to the associate
    process.schedule.associate(process.custom_nanoaod_task)
    
    
    
    
