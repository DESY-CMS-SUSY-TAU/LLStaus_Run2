import FWCore.ParameterSet.Config as cms

from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.genparticles_cff import *
from PhysicsTools.NanoAOD.taus_cff import *
from PhysicsTools.NanoAOD.jetsAK4_CHS_cff import *
#from PhysicsTools.NanoAOD.jets_cff import *

def customize_process_and_associate(process, disTauTagOutputOpt = 1) :
    
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
    
    trk_cond = "hasTrackDetails"
    
    # https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2017#Packed_ParticleFlow_Candidates
    # https://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_12_4_0/doc/html/d8/d79/classpat_1_1PackedCandidate.html
    # lostInnerHits: https://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_12_4_0/doc/html/d8/d79/classpat_1_1PackedCandidate.html#ab9ef9a12f92e02fa61653ba77ee34274
    # fromPV: https://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_12_4_0/doc/html/d8/d79/classpat_1_1PackedCandidate.html#a1e86b4e893738b7cbae410b7f106f339
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
            fromPV                  = Var("fromPV"                              , int       , doc = "isolated track comes from PV"),
            lostInnerHits           = Var("lostInnerHits"                       , int       , doc = "Lost inner hits"),
            hasTrackDetails         = Var("hasTrackDetails"                     , bool      , doc = "True if a bestTrack can be extracted from this Candidate"),
            phiAtVtx                = Var("phiAtVtx"                            , float     , doc = "Phi of the candidate's track at the vertex; this is identical to phi() for the vast majority of the particles, but the two might differ for some of them if the calorimeters had contributed significantly in defining the 4-vector of the particle"),
            dxy                     = Var(f"?{trk_cond}?dxy:-999"                , float     , doc = "dxy w.r.t. associated PV"),
            dxyError                = Var(f"?{trk_cond}?dxyError:-999"           , float     , doc = "Error on dxy"),
            dz                      = Var(f"?{trk_cond}?dzAssociatedPV:-999"     , float     , doc = "dz w.r.t. associated PV"),
            dzError                 = Var(f"?{trk_cond}?dzError:-999"            , float     , doc = "Error on dz"),
            vx                      = Var("vx"                                  , float     , doc = "Vertex x"),
            vy                      = Var("vx"                                  , float     , doc = "Vertex y"),
            vz                      = Var("vz"                                  , float     , doc = "Vertex z"),
        ),
        externalVariables = cms.PSet(
            isTauIdxSignalCand     = ExtVar("isFromTauForPfCand:isTauIdxSignalCand"       , int, doc = "Index of the tau if it belongs to pat::Tau::signalCands(); else -1"),
            isTauIdxIsoCand        = ExtVar("isFromTauForPfCand:isTauIdxIsoCand"          , int, doc = "Index of the tau if it belongs to pat::Tau::isolationCands(); else -1"),
            isTauIdxLeadChHadCand  = ExtVar("isFromTauForPfCand:isTauIdxLeadChHadCand"    , int, doc = "Index of the tau if it is pat::Tau::leadChargedHadrCand(); else -1"),
        )
    )
    
    
    # Unfiltered taus
    process.finalTaus.cut = cms.string("pt > 18")
    process.tauTable.doc = cms.string("slimmedTaus after basic selection (" + process.finalTaus.cut.value()+")")
    
    
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
            #emEnergyInEB                    = Var("emEnergyInEB"                    , float),
            #emEnergyInEE                    = Var("emEnergyInEE"                    , float),
            #emEnergyInHF                    = Var("emEnergyInHF"                    , float),
            energyFractionHadronic          = Var("energyFractionHadronic"          , float),
            #hadEnergyInHB                   = Var("hadEnergyInHB"                   , float),
            #hadEnergyInHE                   = Var("hadEnergyInHE"                   , float),
            #hadEnergyInHF                   = Var("hadEnergyInHF"                   , float),
            #hadEnergyInHO                   = Var("hadEnergyInHO"                   , float),
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
    
    
    ## GenVisTau
    #myGenVisTauTable = genVisTauTable.clone()
    #myGenVisTauTable.variables.vertexX        = Var("vertex.X"      , float)
    #myGenVisTauTable.variables.vertexY        = Var("vertex.Y"      , float)
    #myGenVisTauTable.variables.vertexZ        = Var("vertex.Z"      , float)
    #myGenVisTauTable.variables.vertexRho      = Var("vertex.Rho"    , float)
    #myGenVisTauTable.variables.vertexR        = Var("vertex.R"      , float)
    #
    #process.globalReplace("genVisTauTable", myGenVisTauTable)
    
    if (disTauTagOutputOpt > 0) :
        
        process.disTauTag = cms.EDProducer(
            "DisTauTag",
            graphPath = cms.string("data/particlenet_v1_a27159734e304ea4b7f9e0042baa9e22.pb"),
            #jets = cms.InputTag("finalJets"),
            jets = process.jetTable.src,
            pfCandidates = cms.InputTag('packedPFCandidates'),
            save_inputs  = cms.bool(False)
        )
        
        d_disTauTagVars = {
            "disTauTag_score0":     ExtVar("disTauTag:score0"       , float, doc = "Score 0"),
            "disTauTag_score1":     ExtVar("disTauTag:score1"       , float, doc = "Score 1"),
        }
    
    ##process.jetTable.externalVariables = process.jetTable.externalVariables.clone(
    ##    #disTauTag_score0         = ExtVar("disTauTag:score0"       , float, doc = "Score 0"),
    ##    #disTauTag_score1         = ExtVar("disTauTag:score1"       , float, doc = "Score 1"),
    ##    **d_disTauTagVars
    ##)
    
    
    # Create the task
    if (disTauTagOutputOpt == 0) :
        
        process.custom_nanoaod_task = cms.Task(
            process.lostTrackTable,
            
            process.isFromTauForPfCand,
            process.pfCandTable,
            
            process.caloJetTable,
        )
    
    elif (disTauTagOutputOpt == 1) :
        
        process.jetTable.externalVariables = process.jetTable.externalVariables.clone(**d_disTauTagVars)
        
        process.custom_nanoaod_task = cms.Task(
            process.lostTrackTable,
            
            process.isFromTauForPfCand,
            process.pfCandTable,
            
            process.caloJetTable,
            
            process.disTauTag,
        )
    
    elif (disTauTagOutputOpt == 2) :
        
        process.jetTable.variables = cms.PSet()
        process.jetTable.externalVariables = cms.PSet(**d_disTauTagVars)
        
        process.custom_nanoaod_task = cms.Task(process.disTauTag)
    
    
    
    # Associate the task to the associate
    process.schedule.associate(process.custom_nanoaod_task)
