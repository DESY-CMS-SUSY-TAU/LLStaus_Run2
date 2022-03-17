import aghast
import awkward as ak
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

import glob
from omegaconf import DictConfig, OmegaConf
import hydra

# import ROOT
# ROOT.gROOT.SetBatch(True)

import utils

from coffea.nanoevents import NanoEventsFactory, NanoAODSchema

import matplotlib.pyplot as plt
# import matplotlib as mpl

# import mplhep
# plt.style.use(mplhep.style.ROOT)

class JetMatching(coffea.processor.ProcessorABC):
    def __init__(self):
        self.hists = {}
        
        self.hists["dRJetSTau"] = coffea.hist.Hist(
                "Events",
                coffea.hist.Cat("dataset", "Dataset"),
                coffea.hist.Bin("dRTauSTau", "dR(Tau,STau)", 100, -1.1, 2.5),
                coffea.hist.Bin("dRJetSTau", "dR(Jet,STau)", 100, -1.1, 2.5),
            )
        
        self.hists["dRJetTau"] = coffea.hist.Hist(
                "Events",
                coffea.hist.Cat("dataset", "Dataset"),
                coffea.hist.Bin("dRTauSTau", "dR(Tau,STau)", 100, -1.1, 2.5),
                coffea.hist.Bin("dRJetTau", "dR(Jet,Tau)", 100, -1.1, 2.5),
            )

        self._accumulator = coffea.processor.dict_accumulator(self.hists)

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        out = self.accumulator.identity()

        objects = {}

        # First way of finding pairs
        # objects["Taus_susy"] = events.GenVisTau[
        #     abs(events.GenVisTau.parent.parent.parent.pdgId) == 1000015
        # ]
        # objects["STaus_susy"] = objects["Taus_susy"].parent.parent.parent

        events["GenPart", "motherPdgId"] = events.GenPart.pdgId[events.GenPart.genPartIdxMother]

        print(events.GenPart.motherPdgId)
        print(len(events.GenPart.motherPdgId))
        print(len(events.GenPart.pdgId[0]))
        print(len(events.GenPart.motherPdgId[0]))
        print(len(ak.flatten(events.GenPart.pdgId)))
        print(len(ak.flatten(events.GenPart.motherPdgId)))
        #print( abs(events.GenPart.pdgId) == 1000015 & abs(events.GenPart.motherPdgId) != 1000015)
        print(events.GenPart.hasFlags(["isHardProcess"]))
        print(len(events.GenPart.hasFlags(["isHardProcess"])))
        print(len(events.GenPart.hasFlags(["isHardProcess"])[0]))

        # Second way of calculating pairs
        objects["Taus_susy"] = events.GenVisTau

        objects["STaus_susy"] = events.GenPart[
            abs(events.GenPart.pdgId) == 1000015
            # & events.GenPart.hasFlags(["isHardProcess"])
            #& ak.all(abs(events.GenPart.parent.pdgId) != 1000015, axis=-1)
            #& abs(events.GenPart.pdgId[events.GenPart.genPartIdxMother]) != 1000015
            & abs(events.GenPart.motherPdgId) != 1000015
        ]
        #print(len(ak.all(abs(events.GenPart.children.pdgId) != 1000015, axis=-1)))
        print(len(objects["STaus_susy"]))
        print(len(ak.flatten(objects["STaus_susy"])))
        exit()


        # objects["dRTauSTau"] = objects["STaus_susy"].delta_r(objects["Taus_susy"])
        objects["dRTauSTau"] = objects["Taus_susy"].nearest(objects["STaus_susy"], return_metric=True, threshold=None)[1]
        objects["dRSTauTau"] = objects["STaus_susy"].nearest(objects["Taus_susy"], return_metric=True, threshold=None)[1]

        objects["dRJetSTau"] = objects["STaus_susy"].nearest(events.Jet, return_metric=True, threshold=None)[1]
        objects["dRJetTau"] = objects["Taus_susy"].nearest(events.Jet, return_metric=True, threshold=None)[1]


        objects["dRTauSTau"] = ak.flatten(objects["dRTauSTau"])
        objects["dRSTauTau"] = ak.flatten(objects["dRSTauTau"])
        objects["dRJetSTau"] = ak.flatten(objects["dRJetSTau"])
        objects["dRJetTau"] = ak.flatten(objects["dRJetTau"])


        objects["dRTauSTau"] = ak.fill_none(objects["dRTauSTau"], -1)
        objects["dRSTauTau"] = ak.fill_none(objects["dRSTauTau"], -1)
        objects["dRJetSTau"] = ak.fill_none(objects["dRJetSTau"], -1)
        objects["dRJetTau"] = ak.fill_none(objects["dRJetTau"], -1)

        # print(len(objects["dRSTauTau"]), len(objects["dRJetSTau"]), len(objects["dRTauSTau"]), len(objects["dRJetTau"]))
        # exit()

        out["dRJetSTau"].fill(
            dataset = events.metadata["dataset"],
            dRJetSTau = objects["dRJetSTau"],
            dRTauSTau = objects["dRSTauTau"]
            )
        
        out["dRJetTau"].fill(
            dataset = events.metadata["dataset"],
            dRJetTau = objects["dRJetTau"],
            dRTauSTau = objects["dRTauSTau"]
            )

        return out

    def postprocess(self, accumulator):
        return accumulator


@hydra.main(config_path="../configs", config_name="plot_regionStudy.yaml")
def regionStudy(cfg: DictConfig) -> None:
    '''
    The following script is performing the region study
    for the stau signal nanoAOD sample
    '''

    samples = {}
    for name in cfg.input:
        if cfg.n_files > 0:
            samples[name] = glob.glob(f'{cfg.input[name]}/**/*.root', recursive=True)[:cfg.n_files]
        else:
            samples[name] = glob.glob(f'{cfg.input[name]}/**/*.root', recursive=True)

    os.makedirs(cfg.output, exist_ok=True)

    result_JetMatching = coffea.processor.run_uproot_job(
        samples,
        "Events",
        JetMatching(),
        coffea.processor.iterative_executor,
        {"schema": NanoAODSchema},
    )

    import matplotlib.colors as colors
    # norm=colors.LogNorm(vmin=Z.min(), vmax=Z.max())

    ax = coffea.hist.plot2d(result_JetMatching["dRJetSTau"].sum("dataset"), xaxis='dRJetSTau', patch_opts={"norm":colors.LogNorm(vmin=1, vmax=10000)})
    plt.show()
    plt.savefig(cfg.output+'/dRJetSTau.png')

    ax = coffea.hist.plot2d(result_JetMatching["dRJetTau"].sum("dataset"), xaxis='dRJetTau', patch_opts={"norm":colors.LogNorm(vmin=1, vmax=10000)})
    plt.show()
    plt.savefig(cfg.output+'/dRJetTau.png')

if __name__ == "__main__":
    regionStudy()
    
