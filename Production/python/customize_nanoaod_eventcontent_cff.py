import FWCore.ParameterSet.Config as cms

from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.genparticles_cff import *
from PhysicsTools.NanoAOD.taus_cff import *

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
    myfinalTaus = finalTaus.clone(
        src = cms.InputTag("slimmedTausUpdated"),
        cut = cms.string("pt > 18"),
    )
    
    process.globalReplace("finalTaus", myfinalTaus)
    
    
    # CaloJets
    process.caloJetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src = cms.InputTag("slimmedCaloJets"),
        #cut = cms.string(""),
        cut = cms.string("pt > 15"),
        name= cms.string("CaloJet"),
        doc = cms.string("AK4 calo jets"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table
        variables = cms.PSet(
            P4Vars,
            emEnergyFraction                = Var("emEnergyFraction"                , float),
            emEnergyInEB                    = Var("emEnergyInEB"                    , float),
            emEnergyInEE                    = Var("emEnergyInEE"                    , float),
            emEnergyInHF                    = Var("emEnergyInHF"                    , float),
            energyFractionHadronic          = Var("energyFractionHadronic"          , float),
            hadEnergyInHB                   = Var("hadEnergyInHB"                   , float),
            hadEnergyInHE                   = Var("hadEnergyInHE"                   , float),
            hadEnergyInHF                   = Var("hadEnergyInHF"                   , float),
            hadEnergyInHO                   = Var("hadEnergyInHO"                   , float),
            #maxEInEmTowers                  = Var("maxEInEmTowers"                  , float),
            #maxEInHadTowers                 = Var("maxEInHadTowers"                 , float),
            #towersArea                      = Var("towersArea"                      , float),
            detectorP4pt                    = Var("detectorP4.Pt"                   , float),
            detectorP4eta                   = Var("detectorP4.Eta"                  , float),
            detectorP4phi                   = Var("detectorP4.Phi"                  , float),
            detectorP4mass                  = Var("detectorP4.M"                    , float),
            detectorP4energy                = Var("detectorP4.E"                    , float),
        ),
    )
    
    
    # GenParticles
    myGenParticleTable = genParticleTable.clone()
    myGenParticleTable.variables.vertexX        = Var("vertex.X"      , float)
    myGenParticleTable.variables.vertexY        = Var("vertex.Y"      , float)
    myGenParticleTable.variables.vertexZ        = Var("vertex.Z"      , float)
    myGenParticleTable.variables.vertexRho      = Var("vertex.Rho"    , float)
    myGenParticleTable.variables.vertexR        = Var("vertex.R"      , float)
    
    process.globalReplace("genParticleTable", myGenParticleTable)
    
    
    # Create the task
    process.custom_nanoaod_task = cms.Task(
        process.lostTrackTable,
        
        process.isFromTauForPfCand,
        process.pfCandTable,
        
        process.caloJetTable,
    )
    
    # Associate the task to the associate
    process.schedule.associate(process.custom_nanoaod_task)
