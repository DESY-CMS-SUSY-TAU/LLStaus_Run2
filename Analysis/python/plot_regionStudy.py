import awkward as ak
import coffea
import coffea.hist
import coffea.processor
import os

import glob
from omegaconf import DictConfig, OmegaConf
import hydra

# import ROOT
# ROOT.gROOT.SetBatch(True)

import utils

from coffea.nanoevents import NanoEventsFactory, NanoAODSchema

class JetMatching(coffea.processor.ProcessorABC):
    def __init__(self):
        self.hists = {}
        
        self.hists["dRJetSTau"] = coffea.hist.Hist(
                "Events",
                coffea.hist.Cat("dataset", "Dataset"),
                coffea.hist.Bin("dRTauSTau", "dR(STau,Tau)", 100, 0.0, 3.5),
                coffea.hist.Bin("dRJetSTau", "dR(STau,Jet)", 100, 0.0, 3.5),
            )
        
        self.hists["dRJetTau"] = coffea.hist.Hist(
                "Events",
                coffea.hist.Cat("dataset", "Dataset"),
                coffea.hist.Bin("dRTauSTau", "dR(Tau,STau)", 100, 0.0, 3.5),
                coffea.hist.Bin("dRJetTau", "dR(Tau,Jet)", 100, 0.0, 3.5),
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

        # Second way of calculating pairs
        objects["Taus_susy"] = events.GenVisTau
        objects["STaus_susy"] = events.GenPart[
            (abs(events.GenPart.pdgId) == 1000015)
            & (events.GenPart.hasFlags(["isLastCopy"]))
        ]

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

    import matplotlib.pyplot as plt
    import matplotlib.colors as colors

    for dataset in samples:

        hist_dRJetSTau = result_JetMatching["dRJetSTau"].integrate("dataset", dataset)
        ax = coffea.hist.plot2d(hist_dRJetSTau, xaxis='dRJetSTau', patch_opts={"norm":colors.LogNorm()})
        plt.show()
        plt.savefig(cfg.output+f'/dRJetSTau_{dataset}.png')

        hist_dRJetTau = result_JetMatching["dRJetTau"].integrate("dataset", dataset)
        ax = coffea.hist.plot2d(hist_dRJetTau, xaxis='dRJetTau', patch_opts={"norm":colors.LogNorm()})
        plt.show()
        plt.savefig(cfg.output+f'/dRJetTau_{dataset}.png')

if __name__ == "__main__":
    regionStudy()
    
