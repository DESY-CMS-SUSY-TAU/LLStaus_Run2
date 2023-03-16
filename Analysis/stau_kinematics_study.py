#!/usr/bin/env python3

import aghast
import awkward
import boost_histogram
import coffea
import coffea.hist
import coffea.processor
import dataclasses
import hist
import logging
import numba
import numpy
import os
import sortedcontainers
import typing
import uproot
import uproot3

import ROOT
ROOT.gROOT.SetBatch(True)

#import CMS_lumi, tdrstyle
import utils

from coffea.nanoevents import NanoEventsFactory, NanoAODSchema

@dataclasses.dataclass
class MyProcessor(coffea.processor.ProcessorABC) :
    
    datasets         : typing.List[str]
    
    def __post_init__(self) :
        
        #self.dataset_axis = coffea.hist.Cat("dataset", "dataset")
        self.dataset_axis = hist.axis.StrCategory(self.datasets, growth = True, name = "dataset", label = "dataset")
        
        self._accumulator = sortedcontainers.SortedDict({
            "GenVisTauh": hist.Hist(
                self.dataset_axis,
                hist.axis.Regular(100, 0, 1000, name = "GenVisTauh_e", label = "GenVisTauh_e"),
                hist.axis.Regular(100, 0, 1000, name = "GenVisTauh_pt", label = "GenVisTauh_pt"),
                hist.axis.Regular(200, 0, 2, name = "GenVisTauh_e_by_stau_mass", label = "GenVisTauh_e_by_stau_mass"),
                storage = "weight",
                name = "Counts"
            ),
            
            "GenVisTaul": hist.Hist(
                self.dataset_axis,
                hist.axis.Regular(100, 0, 1000, name = "GenVisTaul_e", label = "GenVisTaul_e"),
                hist.axis.Regular(100, 0, 1000, name = "GenVisTaul_pt", label = "GenVisTaul_pt"),
                hist.axis.Regular(200, 0, 2, name = "GenVisTaul_e_by_stau_mass", label = "GenVisTaul_e_by_stau_mass"),
                storage = "weight",
                name = "Counts"
            ),
        })
    
    
    @property
    def accumulator(self) :
        
        return self._accumulator
    
    
    # we will receive a NanoEvents instead of a coffea DataFrame
    def process(self, events) :
        
        output = self.accumulator
        
        stau_mass = 250.0
        
        GenStau = events.GenPart[
            ((abs(events.GenPart.pdgId) == 1000015) | (abs(events.GenPart.pdgId) == 2000015))
            #& (abs(events.GenPart.pdgId[events.GenPart.genPartIdxMother]) != 1000015)
            & events.GenPart.hasFlags(["isHardProcess"])
            & events.GenPart.hasFlags(["isFirstCopy"])
            & (events.GenPart.mass == stau_mass)
        ]
        
        #print("GenStau.distinctChildren.pdgId:", GenStau.distinctChildren.pdgId)
        #print(GenStau.pdgId)
        #print(GenStau.mass)
        
        GenLsp = events.GenPart[
            (abs(events.GenPart.pdgId) == 1000022)
            #& (abs(events.GenPart.pdgId[events.GenPart.genPartIdxMother]) == 1000015)
            & events.GenPart.hasFlags(["isHardProcess"])
            & events.GenPart.hasFlags(["isFirstCopy"])
            & (events.GenPart.mass == 1)
        ]
        
        sel_idx = (awkward.num(GenStau, axis = 1) == 2) & (awkward.num(GenLsp, axis = 1) == 2)
        events = events[sel_idx]
        
        #distinctChildren
        
        GenTau = events.GenPart[
            (abs(events.GenPart.pdgId) == 15)
            & events.GenPart.hasFlags(["isHardProcess"])
            & events.GenPart.hasFlags(["isFirstCopy"])
            #& ((abs(events.GenPart.distinctParent.pdgId) == 1000015) | (abs(events.GenPart.distinctParent.pdgId) == 2000015))
        ]
        
        GenVisTaul = events.GenPart[
            ((abs(events.GenPart.pdgId) == 11) | (abs(events.GenPart.pdgId) == 13))
            & events.GenPart.hasFlags(["isFirstCopy"])
            #& events.GenPart.hasFlags(["isLastCopy"])
            #& events.GenPart.hasFlags(["isDirectPromptTauDecayProduct"])
            & events.GenPart.hasFlags(["isDirectHardProcessTauDecayProduct"])
            #& (abs(events.GenPart.distinctParent.pdgId) == 15)
            #& ((abs(events.GenPart.distinctParent.distinctParent.pdgId) == 1000015) | (abs(events.GenPart.distinctParent.distinctParent.pdgId) == 2000015))
        ]
        
        print(GenVisTaul)
        print(awkward.num(GenVisTaul, axis = 1))
        
        #sel_idx = (awkward.num(events.GenVisTau, axis = 1) >= 1) | (awkward.num(GenVisTaul, axis = 1) >= 1)
        #sel_idx = (awkward.num(GenVisTaul, axis = 1) >= 1)
        #events = events[sel_idx]
        #GenVisTaul = GenVisTaul[sel_idx]
        
        GenStaul = GenVisTaul.distinctParent.distinctParent
        GenVisTaul_stauRF = GenVisTaul.boost(-GenStaul.boostvec)
        
        #print(awkward.sum(GenStaul.mass == 0, axis = None))
        print(awkward.flatten(GenStaul.pdgId[(abs(GenStaul.pdgId) != 1000015) & (abs(GenStaul.pdgId) != 2000015)]))
        
        #print("GenVisTaul.distinctParent.pdgId:", GenVisTaul.distinctParent.pdgId)
        #print("GenVisTaul.distinctParent.distinctParent.pdgId:", GenVisTaul.distinctParent.distinctParent.pdgId)
        #print("GenVisTau.distinctParent.pdgId:", events.GenPart.pdgId[events.GenVisTau.genPartIdxMother])
        
        #print(events)
        
        #GenTauh = events[awkward.num(events.GenVisTau, axis = 1) >= 1].GenPart[events.GenVisTau.genPartIdxMother]
        GenTauh = events.GenVisTau.parent
        #awkward.drop_none(GenTauh)
        #GenTauh = GenTauh[~awkward.is_none(GenTauh, axis = 1)]
        #GenTauh = GenTauh[(abs(GenTauh.pdgId) == 15)]
        
        GenStauh = GenTauh.distinctParent
        GenStauh = GenStauh[
            ((abs(GenStauh.pdgId) == 1000015) | (abs(GenStauh.pdgId) == 2000015))
            | ((abs(GenStauh.distinctParent.pdgId) == 1000015) | (abs(GenStauh.distinctParent.pdgId) == 2000015))
        ]
        
        
        print("GenTauh.pdgId:", GenTauh.pdgId)
        print("GenStauh.pdgId:", GenStauh.pdgId)
        #print("GenStauh.distinctParent.pdgId:", GenStauh.distinctParent.pdgId)
        print("num(events.GenVisTau):", awkward.num(events.GenVisTau, axis = 1))
        
        # Skip processing as it is an EmptyArray
        if not len(events) :
            
            return output
        
        output["GenVisTauh"].fill(
            dataset = events.metadata["dataset"],
            GenVisTauh_e = awkward.flatten(events.GenVisTau.energy),
            GenVisTauh_pt = awkward.flatten(events.GenVisTau.pt),
            GenVisTauh_e_by_stau_mass = awkward.flatten(events.GenVisTau.pt),
        )
        
        output["GenVisTaul"].fill(
            dataset = events.metadata["dataset"],
            GenVisTaul_e = awkward.flatten(GenVisTaul.energy),
            GenVisTaul_pt = awkward.flatten(GenVisTaul.pt),
            GenVisTaul_e_by_stau_mass = awkward.flatten(GenVisTaul_stauRF.energy / GenStaul.mass),
        )
        
        return output
    
    
    def postprocess(self, accumulator):
        
        return accumulator





def main() :
    
    d_fnamelist = {}
    #d_fnamelist["stau_LH"] = list(numpy.loadtxt("../Production/configs/sourceFiles/SMS-TStauStau_lefthanded_mStau-225to250_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18NanoAODv7-Nano02Apr2020_GridpackScan_102X_upgrade2018_realistic_v21-v1_NANOAODSIM/SMS-TStauStau_lefthanded_mStau-225to250_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18NanoAODv7-Nano02Apr2020_GridpackScan_102X_upgrade2018_realistic_v21-v1_NANOAODSIM.txt", dtype = str))
    #d_fnamelist["stau_RH"] = list(numpy.loadtxt("../Production/configs/sourceFiles/SMS-TStauStau_righthanded_mStau-225to250_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18NanoAODv7-Nano02Apr2020_GridpackScan_102X_upgrade2018_realistic_v21-v1_NANOAODSIM/SMS-TStauStau_righthanded_mStau-225to250_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18NanoAODv7-Nano02Apr2020_GridpackScan_102X_upgrade2018_realistic_v21-v1_NANOAODSIM.txt", dtype = str))
    #d_fnamelist["stau_MM"] = list(numpy.loadtxt("../Production/configs/sourceFiles/SMS-TStauStau_ctau-0p01to10_mStau-250to500_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18NanoAODv7-Nano02Apr2020_GridpackScan_102X_upgrade2018_realistic_v21-v1_NANOAODSIM/SMS-TStauStau_ctau-0p01to10_mStau-250to500_TuneCP5_13TeV-madgraphMLM-pythia8_RunIIAutumn18NanoAODv7-Nano02Apr2020_GridpackScan_102X_upgrade2018_realistic_v21-v1_NANOAODSIM.txt", dtype = str))
    
    #d_fnamelist["stau_LH"] = ["/pnfs/desy.de/cms/tier2/store/user/sobhatta/LongLivedStaus/NanoAOD/SUS-RunIISummer20UL18GEN-stau250_lsp1_ctau100mm_v6/crab_stau250_lsp1_ctau100mm/230311_014339/0000/nanoaod_with-disTauTagScore_1.root"]
    
    d_fnamelist["stau_LH"] = ["SMS-TStauStau_lefthanded_mStau-225to250_TuneCP5_13TeV-madgraphMLM-pythia8.root"]
    d_fnamelist["stau_RH"] = ["SMS-TStauStau_righthanded_mStau-225to250_TuneCP5_13TeV-madgraphMLM-pythia8.root"]
    d_fnamelist["stau_MM"] = ["SMS-TStauStau_ctau-0p01to10_mStau-250to500_TuneCP5_13TeV-madgraphMLM-pythia8.root"]
    #print(d_fnamelist["stau_LH"].shape)
    #print(d_fnamelist["stau_LH"][0: 1])
    
    datasets = sortedcontainers.SortedDict({
        "stau_LH": d_fnamelist["stau_LH"],
        "stau_RH": d_fnamelist["stau_RH"],
        "stau_MM": d_fnamelist["stau_MM"],
    })
    
    
    output = coffea.processor.run_uproot_job(
        datasets,
        "Events",
        MyProcessor(
            datasets = list(datasets.keys())
        ),
        #executor = coffea.processor.iterative_executor,
        executor = coffea.processor.futures_executor,
        executor_args = {
            "schema": NanoAODSchema,
            #"skipbadfiles": True,
            "xrootdtimeout": 600, #sec
            "workers": 10
        },
    )
    
    print(output)
    
    with uproot.recreate("output_stau_kinematics_study.root") as fout:
        
        for key in output:
            
            histo = output[key]
            print(key)
            print(histo.axes)
            #print(histo.axes["dataset"])
            print(histo.__dict__)
            print(histo.axes.__dict__)
            
            for s in histo.axes["dataset"]:
                
                for ax in histo.axes :
                    
                    if ax.name == "dataset" :
                        
                        continue
                        
                    print(s, ax)
                    fout[f"{s}/{ax.name}"] = histo[{"dataset": s}].project(ax.name)
        
    
    return 0


if __name__ == "__main__" :
    
    main()



#/SMS-TStauStau_lefthanded_mStau-225to250_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_GridpackScan_102X_upgrade2018_realistic_v21-v1/NANOAODSIM
#/SMS-TStauStau_righthanded_mStau-225to250_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_GridpackScan_102X_upgrade2018_realistic_v21-v1/NANOAODSIM
#/SMS-TStauStau_ctau-0p01to10_mStau-250to500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_GridpackScan_102X_upgrade2018_realistic_v21-v1/NANOAODSIM

#/SMS-TStauStau_ctau-0p01to10_mStau-90_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIAutumn18NanoAODv7-Nano02Apr2020_GridpackScan_102X_upgrade2018_realistic_v21-v1/NANOAODSIM