import awkward as ak
import coffea
import coffea.hist
import coffea.processor
import matplotlib.pyplot
import os

#import uproot
#uproot.open.defaults["xrootd_handler"] = uproot.source.xrootd.MultithreadedXRootDSource

import glob
from omegaconf import DictConfig, OmegaConf
import hydra

import ROOT
ROOT.gROOT.SetBatch(True)

import utils
import geometry_utils_jit
import geometry_utils

from coffea.nanoevents import NanoEventsFactory, NanoAODSchema

class JetMatching(coffea.processor.ProcessorABC):
    def __init__(self, collections, setups={}):
        
        self.collections = collections
        #self.jet_collection = collections.Jets.name
        self.hists = {}

        self.hists["dR_STau_Jet"] = coffea.hist.Hist(
                "Events",
                coffea.hist.Cat("dataset", "Dataset"),
                coffea.hist.Bin("dR_STau_Tau", "dR(STau,Tau)", *setups.dR_bins),
                coffea.hist.Bin("dR_STau_Jet", "dR(STau,Jet)", *setups.dR_bins),
                coffea.hist.Bin("Lxy","Lxy", *setups.Lxy_bins)
                # coffea.hist.Bin("Lxy","Lxy", setups["rho_bins"])
            )

        self.hists["dR_STau_Track"] = coffea.hist.Hist(
                "Events",
                coffea.hist.Cat("dataset", "Dataset"),
                coffea.hist.Bin("dR_STau_Jet", "dR(STau,Jet)", *setups.dR_bins),
                coffea.hist.Bin("dR_STau_lostTrack","dR(STau,lostTrack)", *setups.dR_bins),
                coffea.hist.Bin("Lxy","Lxy", *setups.Lxy_bins)
                # coffea.hist.Bin("Lxy","Lxy", setups["rho_bins"])
            )
        
        self.hists["dR_STau_pfCand"] = coffea.hist.Hist(
                "Events",
                coffea.hist.Cat("dataset", "Dataset"),
                coffea.hist.Bin("dR_STau_Jet", "dR(STau,Jet)", *setups.dR_bins),
                coffea.hist.Bin("dR_STau_pfCand","dR(STau,pfCand)", *setups.dR_bins),
                coffea.hist.Bin("Lxy","Lxy", *setups.Lxy_bins)
                # coffea.hist.Bin("Lxy","Lxy", setups["rho_bins"])
            )
        
        self.hists["dR_Tau_Jet"] = coffea.hist.Hist(
                "Events",
                coffea.hist.Cat("dataset", "Dataset"),
                coffea.hist.Bin("dR_Tau_STau", "dR(Tau,STau)", *setups.dR_bins),
                coffea.hist.Bin("dR_Tau_Jet", "dR(Tau,Jet)", *setups.dR_bins),
                coffea.hist.Bin("Lxy","Lxy", *setups.Lxy_bins)
                # coffea.hist.Bin("Lxy","Lxy", setups["rho_bins"])
            )

        self._accumulator = coffea.processor.dict_accumulator(self.hists)

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        out = self.accumulator.identity()

        objects = {}
        # events = events[0: 1000]

        events[self.collections.GenVisTaus.name, "vertexX"] = events.GenPart.vertexX[events[self.collections.GenVisTaus.name].genPartIdxMother]
        events[self.collections.GenVisTaus.name, "vertexY"] = events.GenPart.vertexY[events[self.collections.GenVisTaus.name].genPartIdxMother]
        events[self.collections.GenVisTaus.name, "vertexZ"] = events.GenPart.vertexZ[events[self.collections.GenVisTaus.name].genPartIdxMother]

        # First way of finding pairs
        # objects["Taus_susy"] = events.GenVisTau[
        #     abs(events.GenVisTau.parent.parent.parent.pdgId) == 1000015
        # ]
        # objects["STaus_susy"] = objects["Taus_susy"].parent.parent.parent

        # Second way of calculating pairs
        objects["Taus_susy"] = events[self.collections.GenVisTaus.name][eval(self.collections.GenVisTaus.cut.format(name = "events.%s" %(self.collections.GenVisTaus.name)))]
        objects["STaus_susy"] = events.GenPart[
            (abs(events.GenPart.pdgId) == 1000015)
            & (events.GenPart.hasFlags(["isLastCopy"]))
        ]

        objects["Taus_susy","disp"] = objects["Taus_susy"].parent.vertexRho
        objects["STaus_susy","disp"] = objects["STaus_susy"].children[:,:,0].vertexRho

        events[self.collections.Jets.name,"px"] = events[self.collections.Jets.name].x
        events[self.collections.Jets.name,"py"] = events[self.collections.Jets.name].y
        events[self.collections.Jets.name,"pz"] = events[self.collections.Jets.name].z
        events[self.collections.Jets.name,"E"]  = events[self.collections.Jets.name].t
        
        objects["jets"] = events[self.collections.Jets.name][eval(self.collections.Jets.cut.format(name = "events.%s" %(self.collections.Jets.name)))]

        # objects["dRTauSTau"] = objects["STaus_susy"].delta_r(objects["Taus_susy"])
        objects["dR_STau_Tau"] = objects["STaus_susy"].nearest(objects["Taus_susy"], return_metric=True, threshold=None)[1]
        objects["dR_STau_Jet"] = objects["STaus_susy"].nearest(objects["jets"], 
            #metric = lambda v1, v2: geometry_utils_jit.coffea_nearest_metric_deltaR_shiftVertex(v1s = v1, v2s = v2),
            return_metric=True, threshold=None)[1]
        objects['dR_STau_lostTrack'] = objects["STaus_susy"].nearest(events.LostTrack, return_metric=True, threshold=None)[1]
        objects['dR_STau_pfCand'] = objects["STaus_susy"].nearest(events.PFCandidate, return_metric=True, threshold=None)[1]

        objects["dR_Tau_STau"] = objects["Taus_susy"].nearest(objects["STaus_susy"], return_metric=True, threshold=None)[1]
        objects["dR_Tau_Jet"]  = objects["Taus_susy"].nearest(objects["jets"], 
            #metric = lambda v1, v2: geometry_utils_jit.coffea_nearest_metric_deltaR_shiftVertex(v1s = v1, v2s = v2),
            return_metric=True, threshold=None)[1]

        objects["dR_STau_Tau"] = ak.flatten(objects["dR_STau_Tau"])
        objects["dR_STau_Jet"] = ak.flatten(objects["dR_STau_Jet"])
        objects['dR_STau_lostTrack'] = ak.flatten(objects['dR_STau_lostTrack'])
        objects['dR_STau_pfCand'] = ak.flatten(objects['dR_STau_pfCand'])
        objects['stau_disp'] = ak.flatten(objects["STaus_susy","disp"])

        objects["dR_Tau_STau"] = ak.flatten(objects["dR_Tau_STau"])
        objects["dR_Tau_Jet"]  = ak.flatten(objects["dR_Tau_Jet"])
        objects['tau_disp'] = ak.flatten(objects["Taus_susy","disp"])

        objects["dR_STau_Tau"] = ak.fill_none(objects["dR_STau_Tau"], -1)
        objects["dR_STau_Jet"] = ak.fill_none(objects["dR_STau_Jet"], -1)
        objects['dR_STau_lostTrack'] = ak.fill_none(objects["dR_STau_lostTrack"], -1)
        objects['dR_STau_pfCand'] = ak.fill_none(objects["dR_STau_pfCand"], -1)
        objects['stau_disp'] = ak.fill_none(objects['stau_disp'],-1)

        objects["dR_Tau_STau"] = ak.fill_none(objects["dR_Tau_STau"], -1)
        objects["dR_Tau_Jet"]  = ak.fill_none(objects["dR_Tau_Jet"], -1)
        objects['tau_disp'] = ak.fill_none(objects['tau_disp'],-1)

        out["dR_STau_Jet"].fill(
            dataset = events.metadata["dataset"],
            dR_STau_Tau = objects['dR_STau_Tau'],
            dR_STau_Jet = objects['dR_STau_Jet'],
            Lxy = objects['stau_disp'],
            )

        out["dR_STau_Track"].fill(
            dataset = events.metadata["dataset"],
            dR_STau_Jet = objects['dR_STau_Jet'],
            dR_STau_lostTrack = objects['dR_STau_lostTrack'],
            Lxy = objects['stau_disp'],
            )
        
        out["dR_STau_pfCand"].fill(
            dataset = events.metadata["dataset"],
            dR_STau_Jet = objects['dR_STau_Jet'],
            dR_STau_pfCand = objects['dR_STau_pfCand'],
            Lxy = objects['stau_disp'],
            )
   
        out["dR_Tau_Jet"].fill(
            dataset = events.metadata["dataset"],
            dR_Tau_STau = objects["dR_Tau_STau"],
            dR_Tau_Jet = objects["dR_Tau_Jet"],
            Lxy = objects['tau_disp'],
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
        
        #samples[name] = [_file.replace("/pnfs", "root://dcache-cms-xrootd.desy.de://pnfs") for _file in samples[name]]
        #print(samples[name])
    
    os.makedirs(cfg.output, exist_ok=True)

    mySchema = NanoAODSchema

    mySchema.mixins.update({
        "CaloJet": "PtEtaPhiMCollection",
    })

    result_JetMatching = coffea.processor.run_uproot_job(
        samples,
        "Events",
        JetMatching(
            collections = cfg.collections,
            setups = cfg.bin_setups
        ),
        
        #coffea.processor.iterative_executor,
        #{"schema": mySchema},
        
        executor = coffea.processor.futures_executor,
        executor_args = {"schema": mySchema, "workers": 20},
    )

    import matplotlib.pyplot as plt
    import matplotlib.colors as colors

    path_jet = cfg.output+"/match_"+cfg.collections.Jets.name
    os.makedirs(path_jet, exist_ok=True)

    path_pfCand = cfg.output+"/match_pfcand"
    os.makedirs(path_pfCand, exist_ok=True)  

    for dataset in samples:
        
        for Lxy_slice in cfg.bin_setups.Lxy_slice:

            hist_dRJetSTau = result_JetMatching["dR_STau_Jet"].integrate("dataset",dataset).integrate("Lxy",int_range=slice(*Lxy_slice))
            ax = coffea.hist.plot2d(hist_dRJetSTau, xaxis='dR_STau_Jet', patch_opts={"norm":colors.LogNorm()})
            ax.set_title(f'dR(STau, {cfg.collections.Jets.name}) ({dataset}) Lxy({Lxy_slice[0]}-{Lxy_slice[1]})')
            ax.figure.set_dpi(72)
            ax.figure.tight_layout()
            ax.figure.savefig(path_jet+f'/dR_STau_Jet_dataset_{dataset}_Lxy_{Lxy_slice[0]}-{Lxy_slice[1]}.png')
            ax.figure.savefig(path_jet+f'/dR_STau_Jet_dataset_{dataset}_Lxy_{Lxy_slice[0]}-{Lxy_slice[1]}.pdf')
            matplotlib.pyplot.close(ax.figure)

            hist_dRJetTau = result_JetMatching["dR_Tau_Jet"].integrate("dataset",dataset).integrate("Lxy",int_range=slice(*Lxy_slice))
            ax = coffea.hist.plot2d(hist_dRJetTau, xaxis='dR_Tau_Jet', patch_opts={"norm":colors.LogNorm()})
            ax.set_title(f'dR(Tau, {cfg.collections.Jets.name}) ({dataset}) Lxy({Lxy_slice[0]}-{Lxy_slice[1]})')
            ax.figure.set_dpi(72)
            ax.figure.tight_layout()
            ax.figure.savefig(path_jet+f'/dR_Tau_Jet_dataset_{dataset}_Lxy_{Lxy_slice[0]}-{Lxy_slice[1]}.png')
            ax.figure.savefig(path_jet+f'/dR_Tau_Jet_dataset_{dataset}_Lxy_{Lxy_slice[0]}-{Lxy_slice[1]}.pdf')
            matplotlib.pyplot.close(ax.figure)

            hist_dRJetSTau = result_JetMatching["dR_STau_Track"].integrate("dataset",dataset).integrate("Lxy",int_range=slice(*Lxy_slice))
            ax = coffea.hist.plot2d(hist_dRJetSTau, xaxis='dR_STau_lostTrack', patch_opts={"norm":colors.LogNorm()})
            ax.set_title(f'dR(Tau, lostTrack) ({dataset}) Lxy({Lxy_slice[0]}-{Lxy_slice[1]})')
            ax.figure.set_dpi(72)
            ax.figure.tight_layout()
            ax.figure.savefig(path_pfCand+f'/dR_STau_lostTrack_dataset_{dataset})_Lxy_{Lxy_slice[0]}-{Lxy_slice[1]}.png')
            ax.figure.savefig(path_pfCand+f'/dR_STau_lostTrack_dataset_{dataset})_Lxy_{Lxy_slice[0]}-{Lxy_slice[1]}.pdf')
            matplotlib.pyplot.close(ax.figure)

            hist_dRJetSTau = result_JetMatching["dR_STau_pfCand"].integrate("dataset",dataset).integrate("Lxy",int_range=slice(*Lxy_slice))
            ax = coffea.hist.plot2d(hist_dRJetSTau, xaxis='dR_STau_pfCand', patch_opts={"norm":colors.LogNorm()})
            ax.set_title(f'dR(Tau, pfCand) ({dataset}) Lxy({Lxy_slice[0]}-{Lxy_slice[1]})')
            ax.figure.set_dpi(72)
            ax.figure.tight_layout()
            ax.figure.savefig(path_pfCand+f'/dR_STau_pfCand_dataset_{dataset})_Lxy_{Lxy_slice[0]}-{Lxy_slice[1]}.png')
            ax.figure.savefig(path_pfCand+f'/dR_STau_pfCand_dataset_{dataset})_Lxy_{Lxy_slice[0]}-{Lxy_slice[1]}.pdf')
            matplotlib.pyplot.close(ax.figure)

if __name__ == "__main__":
    regionStudy()
    
