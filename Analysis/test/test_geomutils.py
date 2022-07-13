import aghast
import awkward
import boost_histogram
import coffea
import coffea.hist
import coffea.processor
import dataclasses
import logging
import numba
import numpy
import os
import sortedcontainers
import uproot3

import ROOT
ROOT.gROOT.SetBatch(True)

#import CMS_lumi, tdrstyle
import geometry_utils
import utils

from coffea.nanoevents import NanoEventsFactory, NanoAODSchema


newVertex = [20.0, 20.0, 0.0]
detP4_PtEtaPhiM = [50.0, 0.45, 1.9, 10.4]
oldVertex = [0.0, 0.0, 0.0]

rt_detP4 = ROOT.Math.PtEtaPhiMVector(*detP4_PtEtaPhiM)

detP4 = [rt_detP4.Px(), rt_detP4.Py(), rt_detP4.Pz(), rt_detP4.E()]

#l_physP4 = [geometry_utils.physicsP4(
#    newVertex = newVertex,
#    inParticle = detP4,
#    oldVertex = oldVertex,
#)]

#l_physP4 = [geometry_utils.physicsP4(
#    newVertex = newVertex,
#    inParticle = detP4,
#    oldVertex = oldVertex,
#)]

#l_newVertex = [newVertex]
#l_detP4 = [detP4]
#l_oldVertex = [oldVertex]
#
#l_physP4 = geometry_utils.vphysicsP4( 
#    newVertex = l_newVertex,
#    inParticle = l_detP4,
#    oldVertex = l_oldVertex,
#)
#
##print(l_physP4)
#
#for iVec, physP4 in enumerate(l_physP4) :
#    
#    print(physP4)
#    rt_phP4 = ROOT.Math.PxPyPzEVector(*physP4)
#    
#    print("[vec %d] Old: pt %0.4f, eta %+0.4f, phi %+0.4f, mass %0.4f, energy %0.4f" %(iVec, rt_detP4.Pt(), rt_detP4.Eta(), rt_detP4.Phi(), rt_detP4.M(), rt_detP4.E()))
#    print("[vec %d] New: pt %0.4f, eta %+0.4f, phi %+0.4f, mass %0.4f, energy %0.4f" %(iVec, rt_phP4.Pt(), rt_phP4.Eta(), rt_phP4.Phi(), rt_phP4.M(), rt_phP4.E()))
#    print()


@dataclasses.dataclass
class MyProcessor(coffea.processor.ProcessorABC) :
    
    trigger_str         : str       = None
    
    def __post_init__(self) :
        
        self.dataset_axis = coffea.hist.Cat("dataset", "Dataset")
        
        self._accumulator = coffea.processor.dict_accumulator({
            "MET_pt": coffea.hist.Hist(
                "Counts",
                self.dataset_axis,
                coffea.hist.Bin("MET_pt", "MET_pt", 100, 0, 1000),
            ),
        })
    
    
    @property
    def accumulator(self) :
        
        return self._accumulator
    
    
    # we will receive a NanoEvents instead of a coffea DataFrame
    def process(self, events) :
        
        if (self.trigger_str is not None and len(self.trigger_str)) :
            
            #events = events[events.HLT.PFMETNoMu110_PFMHTNoMu110_IDTight | events.HLT.PFMETNoMu120_PFMHTNoMu120_IDTight | ]
            events = events[eval(self.trigger_str)]
        
        #events = events[0: 5]
        
        events["GenVisTau", "vertexX"] = events.GenPart.vertexX[events.GenVisTau.genPartIdxMother]
        events["GenVisTau", "vertexY"] = events.GenPart.vertexY[events.GenVisTau.genPartIdxMother]
        events["GenVisTau", "vertexZ"] = events.GenPart.vertexZ[events.GenVisTau.genPartIdxMother]
        
        output = self.accumulator.identity()
        
        #print(events.GenVisTau.x)
        
        z_caloJetP4 = awkward.zip([events.CaloJet.x, events.CaloJet.y, events.CaloJet.z, events.CaloJet.t])
        l_caloJetP4 = awkward.to_list(z_caloJetP4)
        
        #print(l_caloJetP4)
        print(l_caloJetP4[0])
        
        #l_physP4 = geometry_utils.vphysicsP4(
        #    inParticle = l_caloJetP4[0],
        #    oldVertex = (0.0, 0.0, 0.0),
        #    newVertex = (10.0, -10.0, 0.0)#, (-20.0, 5.0, 5.0)],
        #)
        #
        #print(l_physP4)
        
        print("*"*10, len(events))
        print("$"*10, len(events.GenVisTau))
        print("$"*10, len(events.CaloJet))
        print("#"*10, len(events.GenVisTau[0]))
        print("#"*10, len(events.CaloJet[0]))
        print("%"*10, len(events.GenVisTau[1]))
        print("%"*10, len(events.CaloJet[1]))
        
        myCaloJets = events.CaloJet[events.CaloJet.pt > 100]
        
        nearest = events.GenVisTau.nearest(
            #events.CaloJet,
            myCaloJets,
            
            #axis=1,
            #metric = lambda v1, v2: geometry_utils.get_deltaR_shiftVertex(v1s = v1, v2s = v2, oldVertex = (0.0, 0.0, 0.0), newVertex = (10.0, -10.0, 0.0)),
            
            #metric = lambda v1, v2: geometry_utils.get_deltaR_shiftVertex(v1s = v1, v2s = v2, oldVertex = (0.0, 0.0, 0.0), newVertex = (0.0, 0.0, 0.0)),
            #metric = lambda v1, v2: geometry_utils.coffea_nearest_metric_deltaR_shiftVertex(v1s = v1, v2s = v2, oldVertex = (0.0, 0.0, 0.0), newVertex = (0.0, 0.0, 0.0)),
            
            #metric = lambda v1, v2: geometry_utils.get_deltaR_shiftVertex(v1s = v1, v2s = v2, oldVertex = (0.0, 0.0, 0.0)),
            metric = lambda v1, v2: geometry_utils.coffea_nearest_metric_deltaR_shiftVertex(v1s = v1, v2s = v2, oldVertex = (0.0, 0.0, 0.0)),
            
            return_metric = True,
        )
        
        print(nearest[0])
        print(nearest[1])
        
        return output
    
    
    def postprocess(self, accumulator):
        
        return accumulator


def main() :
    
    datasets = sortedcontainers.SortedDict({
        "test_sample": ["tmp/nanoaod.root"],
    })
    
    #trigger_str = ""
    
    trigger_str = (
        "HLT_PFMET120_PFMHT120_IDTight"
        " | HLT_PFMET130_PFMHT130_IDTight "
        " | HLT_PFMET140_PFMHT140_IDTight "
        " | HLT_PFMETNoMu110_PFMHTNoMu110_IDTight "
        " | HLT_PFMETNoMu120_PFMHTNoMu120_IDTight "
        " | HLT_PFMETNoMu130_PFMHTNoMu130_IDTight "
        " | HLT_PFMETNoMu140_PFMHTNoMu140_IDTight "
        
        " | HLT_IsoMu24 "
        #" | HLT_TkMu100 "
        
        #" | HLT_IsoMu27_MET90 "
        
        #" | HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg "
        #" | HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg "
        #" | HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg "
        #" | HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg "
        
        " | HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100 "
        " | HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110 "
        " | HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120 "
        " | HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130 "
        " | HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140 "
        
        " | HLT_MET105_IsoTrk50 "
        " | HLT_MET120_IsoTrk50 "
        
        #" | HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1 "
        #" | HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1 "
        
        #" | HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL"
        
        #" | HLT_IsoTrackHB "
    )
    
    trigger_str = trigger_str.replace("HLT_", "events.HLT.")
    
    mySchema = NanoAODSchema
    
    mySchema.mixins.update({
        "CaloJet": "PtEtaPhiMCollection",
    })
    
    output_num = coffea.processor.run_uproot_job(
        datasets,
        "Events",
        MyProcessor(
            #trigger_str = trigger_str,
        ),
        coffea.processor.iterative_executor,
        #{"schema": NanoAODSchema},
        {"schema": mySchema},
    )
    
    
    return 0


if __name__ == "__main__" :
    
    main()

