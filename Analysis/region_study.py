import awkward as ak
import coffea
import coffea.hist
import coffea.processor
import matplotlib.pyplot
import os
import aghast

#import uproot
#uproot.open.defaults["xrootd_handler"] = uproot.source.xrootd.MultithreadedXRootDSource

import glob
from omegaconf import DictConfig, OmegaConf
import hydra
import sys

import ROOT
ROOT.gROOT.SetBatch(True)

import numpy as np

import utils.utils as utils
# import utils.geometry_utils as geometry_utils_jit
import utils.geometry_utils

from coffea.nanoevents import NanoEventsFactory, NanoAODSchema

def conv_to_root(coffea_hist, name):
    w, data = coffea_hist.to_boost().to_numpy()
    return aghast.to_root(
            aghast.from_numpy((w.astype(int), data)),
            name
           )

def stand_arr(array):
    return ak.fill_none( ak.flatten(array), -1)


class JetdR2D(object):

    @staticmethod
    def __init__(hists, setups={}):

        hists["dR_STau_Jet"] = coffea.hist.Hist(
            "Events",
            coffea.hist.Cat("dataset", "Dataset"),
            coffea.hist.Bin("dR_STau_Tau", "dR(STau,Tau)", *setups.dR_bins),
            coffea.hist.Bin("dR_STau_Jet", "dR(STau,Jet)", *setups.dR_bins),
            coffea.hist.Bin("Lxy","Lxy", np.array(setups.Lxy_bins))
            # coffea.hist.Bin("Lxy","Lxy", setups["rho_bins"])
        )
        hists["dR_STau_Track"] = coffea.hist.Hist(
            "Events",
            coffea.hist.Cat("dataset", "Dataset"),
            coffea.hist.Bin("dR_STau_Jet", "dR(STau,Jet)", *setups.dR_bins),
            coffea.hist.Bin("dR_STau_lostTrack","dR(STau,lostTrack)", *setups.dR_bins),
            coffea.hist.Bin("Lxy","Lxy", np.array(setups.Lxy_bins))
            # coffea.hist.Bin("Lxy","Lxy", setups["rho_bins"])
        )
        
        hists["dR_STau_pfCand"] = coffea.hist.Hist(
            "Events",
            coffea.hist.Cat("dataset", "Dataset"),
            coffea.hist.Bin("dR_STau_Jet", "dR(STau,Jet)", *setups.dR_bins),
            coffea.hist.Bin("dR_STau_pfCand","dR(STau,pfCand)", *setups.dR_bins),
            coffea.hist.Bin("Lxy","Lxy", np.array(setups.Lxy_bins))
            # coffea.hist.Bin("Lxy","Lxy", setups["rho_bins"])
        )
        
        hists["dR_Tau_Jet"] = coffea.hist.Hist(
            "Events",
            coffea.hist.Cat("dataset", "Dataset"),
            coffea.hist.Bin("dR_Tau_STau", "dR(Tau,STau)", *setups.dR_bins),
            coffea.hist.Bin("dR_Tau_Jet", "dR(Tau,Jet)", *setups.dR_bins),
            coffea.hist.Bin("Lxy","Lxy", np.array(setups.Lxy_bins))
            # coffea.hist.Bin("Lxy","Lxy", setups["rho_bins"])
        )
    
    @staticmethod
    def process_events(out, objects, events):

        # objects["dRTauSTau"] = objects["genSUSYTaus"].delta_r(objects["genTau"])
        objects["dR_STau_Tau"] = objects["genSUSYTaus"].nearest(objects["genTau"], return_metric=True, threshold=None)[1]
        objects["dR_STau_Jet"] = objects["genSUSYTaus"].nearest(objects["jets"], 
            #metric = lambda v1, v2: geometry_utils_jit.coffea_nearest_metric_deltaR_shiftVertex(v1s = v1, v2s = v2),
            return_metric=True, threshold=None)[1]
        
        objects['dR_STau_lostTrack'] = objects["genSUSYTaus"].nearest(events.LostTrack, return_metric=True, threshold=None)[1]
        objects['dR_STau_pfCand'] = objects["genSUSYTaus"].nearest(events.PFCandidate, return_metric=True, threshold=None)[1]

        objects["dR_Tau_STau"] = objects["genTau"].nearest(objects["genSUSYTaus"], return_metric=True, threshold=None)[1]
        objects["dR_Tau_Jet"]  = objects["genTau"].nearest(objects["jets"], 
            #metric = lambda v1, v2: geometry_utils_jit.coffea_nearest_metric_deltaR_shiftVertex(v1s = v1, v2s = v2),
            return_metric=True, threshold=None)[1]

        objects["dR_STau_Tau"] = ak.flatten(objects["dR_STau_Tau"])
        objects["dR_STau_Jet"] = ak.flatten(objects["dR_STau_Jet"])
        objects['dR_STau_lostTrack'] = ak.flatten(objects['dR_STau_lostTrack'])
        objects['dR_STau_pfCand'] = ak.flatten(objects['dR_STau_pfCand'])
        objects['stau_disp'] = ak.flatten(objects["genSUSYTaus","disp"])

        objects["dR_Tau_STau"] = ak.flatten(objects["dR_Tau_STau"])
        objects["dR_Tau_Jet"]  = ak.flatten(objects["dR_Tau_Jet"])
        objects['tau_disp'] = ak.flatten(objects["genTau","disp"])

        objects["dR_STau_Tau"] = ak.fill_none(objects["dR_STau_Tau"], -1)
        objects["dR_STau_Jet"] = ak.fill_none(objects["dR_STau_Jet"], -1)
        objects['dR_STau_lostTrack'] = ak.fill_none(objects["dR_STau_lostTrack"], -1)
        objects['dR_STau_pfCand'] = ak.fill_none(objects["dR_STau_pfCand"], -1)
        objects['stau_disp'] = ak.fill_none(objects['stau_disp'],-1)

        objects["dR_Tau_STau"] = ak.fill_none(objects["dR_Tau_STau"], -1)
        objects["dR_Tau_Jet"]  = ak.fill_none(objects["dR_Tau_Jet"], -1)
        objects['tau_disp'] = ak.fill_none(objects['tau_disp'],-1)

        # print(objects['dR_STau_Jet'])
        # print(objects['stau_disp'])
        # exit()

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

    

class EfficiencyStudy(object):

    @staticmethod
    def __init__(hists, setups={}):

        hists["gen_tau_pt"] = coffea.hist.Hist(
            "Events",
            coffea.hist.Cat("dataset", "Dataset"),
            coffea.hist.Bin("gen_tau_pt", "efficiency", np.array(setups.gen_jet_pt)),
            # coffea.hist.Bin("Lxy","Lxy", *setups.Lxy_bins)
            coffea.hist.Bin("Lxy","Lxy", np.array(setups.Lxy_bins))
        )
        
        hists["gen_tau_pt_reco"] = coffea.hist.Hist(
            "Events",
            coffea.hist.Cat("dataset", "Dataset"),
            coffea.hist.Bin("gen_tau_pt", "efficiency", np.array(setups.gen_jet_pt)),
            # coffea.hist.Bin("Lxy","Lxy", *setups.Lxy_bins)
            coffea.hist.Bin("Lxy","Lxy", np.array(setups.Lxy_bins))
        )

        hists["gen_jet_pt"] = coffea.hist.Hist(
            "Events",
            coffea.hist.Cat("dataset", "Dataset"),
            coffea.hist.Bin("gen_tau_pt", "efficiency", np.array(setups.gen_jet_pt)),
            # coffea.hist.Bin("Lxy","Lxy", *setups.Lxy_bins)
            coffea.hist.Bin("Lxy","Lxy", np.array(setups.Lxy_bins))
        )

        hists["gen_jet_pt_reco"] = coffea.hist.Hist(
            "Events",
            coffea.hist.Cat("dataset", "Dataset"),
            coffea.hist.Bin("gen_tau_pt", "efficiency", np.array(setups.gen_jet_pt)),
            # coffea.hist.Bin("Lxy","Lxy", *setups.Lxy_bins)
            coffea.hist.Bin("Lxy","Lxy", np.array(setups.Lxy_bins))
        )

        # hists["gen_energy_ratio"] = coffea.hist.Hist(
        #     "Events",
        #     coffea.hist.Cat("dataset", "Dataset"),
        #     coffea.hist.Bin("energy_ratio", "energy_ratio", *setups.energy_ratio),
        #     coffea.hist.Bin("Lxy","Lxy", *setups.Lxy_bins)
        # )

        # hists["tau_childr_energy_ratio"] = coffea.hist.Hist(
        #     "Events",
        #     coffea.hist.Cat("dataset", "Dataset"),
        #     coffea.hist.Bin("energy_ratio", "energy_ratio", *setups.energy_ratio_children),
        #     coffea.hist.Bin("Lxy","Lxy", *setups.Lxy_bins)
        # )

        # hists["tau_childr_visTau_energy_ratio"] = coffea.hist.Hist(
        #     "Events",
        #     coffea.hist.Cat("dataset", "Dataset"),
        #     coffea.hist.Bin("energy_ratio", "energy_ratio", *setups.energy_ratio_children),
        #     coffea.hist.Bin("Lxy","Lxy", *setups.Lxy_bins)
        # )

        # hists["pt_ratio_genJet"] = coffea.hist.Hist(
        #     "Events",
        #     coffea.hist.Cat("dataset", "Dataset"),
        #     coffea.hist.Bin("pt_ratio", "pt_reco/pt_gen", *(setups.pt_ratio[:3])),
        #     coffea.hist.Bin("pt_gen", "pt_gen", *(setups.pt_ratio[-3:])),
        #     coffea.hist.Bin("Lxy","Lxy", *setups.Lxy_bins)
        # )

        hists["pt_ratio_visTau"] = coffea.hist.Hist(
            "Events",
            coffea.hist.Cat("dataset", "Dataset"),
            coffea.hist.Bin("pt_ratio", "pt_reco/pt_gen", *(setups.pt_ratio[:3])),
            coffea.hist.Bin("pt_gen", "pt_gen", *(setups.pt_ratio[-3:])),
            coffea.hist.Bin("Lxy","Lxy", np.array(setups.Lxy_bins))
        )

        hists["pt_resolution_visTau"] = coffea.hist.Hist(
            "Events",
            coffea.hist.Cat("dataset", "Dataset"),
            coffea.hist.Bin("pt_resolution", "pt_reco-pt_gen/pt_gen", *(setups.pt_resolution[:3])),
            coffea.hist.Bin("pt_gen", "pt_gen", *(setups.pt_resolution[-3:])),
            coffea.hist.Bin("Lxy","Lxy", np.array(setups.Lxy_bins))
        )

    @staticmethod
    def process_events(out, objects, events):
        
        # genTaus that have matching reco tau

        hpsTaus = objects["genTau"].nearest(objects["hpsTaus"], return_metric=False, threshold=0.4)
        objects["match_genTau_to_recoTau"] = objects["genTau"][(~ak.is_none(hpsTaus.energy,-1))]

        out["gen_tau_pt"].fill(
            dataset = events.metadata["dataset"],
            gen_tau_pt = stand_arr( objects["genTau"].pt ),
            Lxy = stand_arr( objects["genTau"].parent.vertexRho ),
            )
        
        out["gen_tau_pt_reco"].fill(
            dataset = events.metadata["dataset"],
            gen_tau_pt = stand_arr( objects["match_genTau_to_recoTau"].pt ),
            Lxy = stand_arr( objects["match_genTau_to_recoTau"].parent.vertexRho ),
            )
        
        genJets = objects["genTau"].nearest(objects["genJets"], return_metric=False, threshold=0.4)
        objects["match_genTau_to_genJet"] = objects["genTau"][(~ak.is_none(genJets.energy,-1))]

        out["gen_jet_pt"].fill(
            dataset = events.metadata["dataset"],
            gen_tau_pt = stand_arr( objects["match_genTau_to_genJet"].pt ),
            Lxy = stand_arr( objects["match_genTau_to_genJet"].parent.vertexRho ),
            )
        
        jets = objects["genTau"].nearest(objects["jets"], return_metric=False, threshold=0.4)
        objects["match_genTau_to_recoJet"] = objects["genTau"][(~ak.is_none(jets.energy,-1))]
   
        out["gen_jet_pt_reco"].fill(
            dataset = events.metadata["dataset"],
            gen_tau_pt = stand_arr( objects["match_genTau_to_recoJet"].pt ),
            Lxy = stand_arr( objects["match_genTau_to_recoJet"].parent.vertexRho ),
            )
        
        objects["match_recoJet_to_genTau"] = jets[(~ak.is_none(jets.energy,-1))]
        
        out["pt_ratio_visTau"].fill(
            dataset = events.metadata["dataset"],
            pt_ratio = stand_arr( objects["match_recoJet_to_genTau"].pt /objects["match_genTau_to_recoJet"].pt  ),
            pt_gen = stand_arr( objects["match_genTau_to_recoJet"].pt ),
            Lxy = stand_arr( objects["match_genTau_to_recoJet"].parent.vertexRho ),
            )

        out["pt_resolution_visTau"].fill(
            dataset = events.metadata["dataset"],
            pt_resolution = stand_arr( (objects["match_recoJet_to_genTau"].pt - objects["match_genTau_to_recoJet"].pt) /objects["match_genTau_to_recoJet"].pt  ),
            pt_gen = stand_arr( objects["match_genTau_to_recoJet"].pt ),
            Lxy = stand_arr( objects["match_genTau_to_recoJet"].parent.vertexRho ),
        )
        
        
        '''

        # genJets matched to the visible tau - matched_genJet, visible taus that mathced to this jets match_genTau_to_genJet
        objects["match_genTau_to_genJet"] = objects["genJets"].nearest(objects["genTau"], return_metric=False, threshold=0.4) # genTaus matched to jets
        objects["match_genJet_to_genTau"] = objects["genJets"][(objects["match_genTau_to_genJet"].energy != None)] # genJets that have matched genTaus

        # recoJets matched to the visible tau - match_recoJet, visible taus that mathced to this jets match_genTau_to_recoJet
        objects["match_genTau_to_recoJet"] = objects["jets"].nearest(objects["genTau"], return_metric=False, threshold=0.4) # genTaus matched to reco jets
        objects["match_recoJet"] = objects["jets"][(objects["match_genTau_to_recoJet"].energy != None)] # Jets that have matched genTaus

        # recoJet matched to the genJet - match_jets, genJet that matched to this jets match_genTau_to_recoJet
        objects["match_recoJet_to_genJet"] = objects["match_genJet_to_genTau"].nearest(objects["jets"], return_metric=False, threshold=0.4)
        objects["match_genJet_to_recoJet_and_genTau"] = objects["match_genJet_to_genTau"][(objects["match_recoJet_to_genJet"].energy != None)]
        
        objects["match_genTau_to_genJet_and_recoJet"] = objects["match_genTau_to_genJet"][(objects["match_recoJet_to_genJet"].energy != None)]

        objects["clone_Taus_susy"] = objects["genTau"]


        objects["match_genJet_to_genTau", "jet_E"] = objects["match_genJet_to_genTau"].energy
        objects["match_genJet_to_genTau", "tau_E"] = objects["match_genTau_to_genJet"].energy[(objects["match_genTau_to_genJet"].energy != None)]
        objects["match_genJet_to_genTau", "tau_disp"] = objects["match_genTau_to_genJet"].parent.vertexRho[(objects["match_genTau_to_genJet"].energy != None)]
        objects["match_genJet_to_genTau", "ratio_E"]  = objects["match_genJet_to_genTau", "jet_E"] / objects["match_genJet_to_genTau", "tau_E"]
        
        # Select tau (that mathed to the jet children) that not None 
        tau = objects["match_genTau_to_genJet"].parent
        tau_children = tau[(~ak.is_none(tau.pdgId,-1))].children
        # Select only visible tau children (drop neutrino)
        tau_children = tau_children[(abs(tau_children.pdgId)!=16)]

        # matched genJet without jets that does not match to anything
        match_genJet = objects["match_genJet_to_genTau"][(~ak.is_none(objects["match_genJet_to_genTau"].energy,-1))]

        # print(match_genJet.energy.to_list())
        # print(tau_children.pdgId.to_list())
        # print(len(match_genJet), len(tau_children))

        dR = tau_children.delta_r(match_genJet)
        # visible tau children that within 0.4 cone to the jet
        tau_children_in_cone = tau_children[(dR<0.4)]
        tau_children_in_cone_E = ak.sum(tau_children_in_cone.energy, axis=2)
        tau_children_E = ak.sum(tau_children.energy, axis=2)

        tau_vis_E = objects["match_genTau_to_genJet"][(~ak.is_none(tau.energy,-1))].energy
        tau_vis_disp = objects["match_genTau_to_genJet"][(~ak.is_none(tau.energy,-1))].parent.vertexRho

        # print(tau_tau_children_in_cone_E.to_list())
        # print(tau_tau_children_E.to_list())

        # print(len(tau_children_E), len(tau_children_in_cone_E), len(tau_vis_E))

        tau_childr_energy_ratio = tau_children_E/tau_children_in_cone_E 
        tau_childr_visTau_energy_ratio = tau_vis_E/tau_children_in_cone_E

        # for E1, E2 in zip(tau_vis_E.to_list(), tau_children_in_cone_E.to_list()):
        #     print(E1, E2)

        # exit()
        # print(len(tau_childr_energy_ratio), len(tau_childr_visTau_energy_ratio))
        # print(len(objects["match_genJet_to_genTau", "tau_disp"]))
        
        
        out["gen_tau_pt"].fill(
            dataset = events.metadata["dataset"],
            gen_tau_pt = stand_arr( objects["genTau"].pt ),
            Lxy = stand_arr( objects["genTau"].parent.vertexRho ),
            )

        out["gen_jet_pt"].fill(
            dataset = events.metadata["dataset"],
            gen_tau_pt = stand_arr( objects["match_genTau_to_genJet"].pt ),
            Lxy = stand_arr( objects["match_genTau_to_genJet"].parent.vertexRho ),
            )
   
        out["gen_jet_pt_reco"].fill(
            dataset = events.metadata["dataset"],
            gen_tau_pt = stand_arr( objects["match_genTau_to_recoJet"].pt ),
            Lxy = stand_arr( objects["match_genTau_to_recoJet"].parent.vertexRho ),
            )

        # out["gen_energy_ratio"].fill(
        #     dataset = events.metadata["dataset"],
        #     energy_ratio = stand_arr( objects["match_genJet_to_genTau", "ratio_E"] ),
        #     Lxy = stand_arr( objects["match_genJet_to_genTau", "tau_disp"] ),
        #     )

        # out["tau_childr_energy_ratio"].fill(
        #     dataset = events.metadata["dataset"],
        #     energy_ratio = stand_arr( tau_childr_energy_ratio ),
        #     Lxy = stand_arr( tau_vis_disp ),
        #     )

        # out["tau_childr_visTau_energy_ratio"].fill(
        #     dataset = events.metadata["dataset"],
        #     energy_ratio = stand_arr( tau_childr_visTau_energy_ratio ),
        #     Lxy = stand_arr( tau_vis_disp ),
        #     )

        # ######

        # out["pt_ratio_genJet"].fill(
        #     dataset = events.metadata["dataset"],
        #     pt_ratio = stand_arr( objects["match_recoJet_to_genJet"].pt / objects["match_genJet_to_recoJet_and_genTau"].pt  ),
        #     pt_gen = stand_arr( objects["match_genJet_to_recoJet_and_genTau"].pt ),
        #     Lxy = stand_arr( objects["match_genTau_to_genJet_and_recoJet"].parent.vertexRho ),
        #     )

        # out["pt_ratio_visTau"].fill(
        #     dataset = events.metadata["dataset"],
        #     pt_ratio = stand_arr( objects["match_recoJet"].pt /objects["match_genTau_to_recoJet"].pt  ),
        #     pt_gen = stand_arr( objects["match_genTau_to_recoJet"].pt ),
        #     Lxy = stand_arr( objects["match_genTau_to_recoJet"].parent.vertexRho ),
        #     )

        # out["pt_resolution_visTau"].fill(
        #     dataset = events.metadata["dataset"],
        #     pt_resolution = stand_arr( (objects["match_recoJet"].pt - objects["match_genTau_to_recoJet"].pt) /objects["match_genTau_to_recoJet"].pt  ),
        #     pt_gen = stand_arr( objects["match_genTau_to_recoJet"].pt ),
        #     Lxy = stand_arr( objects["match_genTau_to_recoJet"].parent.vertexRho ),
        # )
        
        '''

class JetMatching(coffea.processor.ProcessorABC):
    def __init__(self, cfg, tag_ids_files):
        
        self.mode = cfg.mode
        self.collections = cfg.collections
        self.tag_ids_files = tag_ids_files
        self.hists = {}

        if 'jet_dR_matching' in self.mode:
            JetdR2D(self.hists, cfg.bin_setups)
        if 'efficiency_study' in self.mode:
            EfficiencyStudy(self.hists, cfg.eff_setups)

        self._accumulator = coffea.processor.dict_accumulator(self.hists)

    @property
    def accumulator(self):
        return self._accumulator

    def get_jets_score(self, _metadata): 
        '''
        parse distautag score from separate files
        '''
       
        # file_path = self.tag_ids_files[_metadata['dataset']][_metadata['filename']]
        # print(file_path)
        # print(_metadata['entrystart'])
        # print(_metadata['entrystop'])
        # exit()

        pass
        
    def process(self, events):
        out = self.accumulator.identity()

        if self.tag_ids_files:
            self.get_jets_score(events.metadata)

        objects = {}
        # events = events[:20]

        events[self.collections.GenVisTaus.name, "vertexX"] = events.GenPart.vertexX[events[self.collections.GenVisTaus.name].genPartIdxMother]
        events[self.collections.GenVisTaus.name, "vertexY"] = events.GenPart.vertexY[events[self.collections.GenVisTaus.name].genPartIdxMother]
        events[self.collections.GenVisTaus.name, "vertexZ"] = events.GenPart.vertexZ[events[self.collections.GenVisTaus.name].genPartIdxMother]

        # First way of finding pairs
        # objects["genTau"] = events.GenVisTau[
        #     abs(events.GenVisTau.parent.parent.parent.pdgId) == 1000015
        # ]
        # objects["genSUSYTaus"] = objects["genTau"].parent.parent.parent

        # Second way of calculating pairs
        objects["genTau"] = events[self.collections.GenVisTaus.name][eval(self.collections.GenVisTaus.cut.format(name = "events.%s" %(self.collections.GenVisTaus.name)))]
        objects["genSUSYTaus"] = events.GenPart[ (abs(events.GenPart.pdgId) == 1000015) & (events.GenPart.hasFlags(["isLastCopy"])) ]

        objects["genTau","disp"] = objects["genTau"].parent.vertexRho
        objects["genSUSYTaus","disp"] = objects["genSUSYTaus"].children[:,:,0].vertexRho

        events[self.collections.Jets.name,"px"] = events[self.collections.Jets.name].x
        events[self.collections.Jets.name,"py"] = events[self.collections.Jets.name].y
        events[self.collections.Jets.name,"pz"] = events[self.collections.Jets.name].z
        events[self.collections.Jets.name,"E"]  = events[self.collections.Jets.name].t
        
        objects["jets"] = events[self.collections.Jets.name][eval(self.collections.Jets.cut.format(name = "events.%s" %(self.collections.Jets.name)))]
        objects["genJets"] = events[self.collections.GenJets.name][eval(self.collections.GenJets.cut.format(name = "events.%s" %(self.collections.GenJets.name)))]
        objects["hpsTaus"] = events[self.collections.Tau.name][eval(self.collections.Tau.cut.format(name = "events.%s" %(self.collections.Tau.name)))]

        # print(objects["genTau"].energy)

        if "jet_dR_matching" in self.mode:
            JetdR2D.process_events(out, objects, events)

        if 'efficiency_study' in self.mode:
            EfficiencyStudy.process_events(out, objects, events)

        return out

    def postprocess(self, accumulator):
        return accumulator


@hydra.main(config_path="./configs", config_name="region_study.yaml")
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
        
    # to apply id from different file:
    id_scores = {}
    if cfg.input_disID:
        print("Apply DisTauTag scrore from file!")
        for name in samples.keys():
            score_files = glob.glob(f'{cfg.input_disID[name]}/**/*.root', recursive=True)
            id_scores[name] = {}
            for file in samples[name]:
                base = os.path.basename(file)
                file_id = base.replace("nanoaod_with-disTauTagScore", "nanoaod_only-disTauTagScore")
                matching = [s for s in score_files if file_id in s]
                if len(matching) != 1:
                    raise RuntimeError("Matching error (no or many matches): ", matching)
                id_scores[name][file] = matching[0]

    os.makedirs(cfg.output, exist_ok=True)

    mySchema = NanoAODSchema

    mySchema.mixins.update({
        "CaloJet": "PtEtaPhiMCollection",
    })

    result_JetMatching = coffea.processor.run_uproot_job(
        samples,
        "Events",
        JetMatching(
            cfg = cfg,
            tag_ids_files = id_scores if cfg.input_disID else False,
        ),
        
        executor = coffea.processor.futures_executor,
        executor_args = {"schema": mySchema, "workers": 20},
    )

    import matplotlib.pyplot as plt
    import matplotlib.colors as colors

    if 'jet_dR_matching' in cfg.mode:

        path_jet = cfg.output+"/match_"+cfg.collections.Jets.name
        os.makedirs(path_jet, exist_ok=True)

        path_pfCand = cfg.output+"/match_pfcand"
        os.makedirs(path_pfCand, exist_ok=True)  

        for dataset in samples:
            
            for Lxy_slice in cfg.bin_setups.Lxy_slice:
         
                dR_STau_Jet = result_JetMatching["dR_STau_Jet"].integrate("dataset",dataset).integrate("Lxy",int_range=slice(*Lxy_slice))
                ax = coffea.hist.plot2d(dR_STau_Jet, xaxis='dR_STau_Jet', patch_opts={"norm":colors.LogNorm()})
                ax.set_title(f'dR(STau, {cfg.collections.Jets.name}) ({dataset}) Lxy({Lxy_slice[0]}-{Lxy_slice[1]})')
                ax.figure.set_dpi(72)
                ax.figure.tight_layout()
                ax.figure.savefig(path_jet+f'/dR_STau_Jet_Tau_dataset_{dataset}_Lxy={Lxy_slice[0]}_{Lxy_slice[1]}.png')
                # ax.figure.savefig(path_jet+f'/dR_STau_Jet_dataset_{dataset}_Lxy_{Lxy_slice[0]}-{Lxy_slice[1]}.pdf')
                matplotlib.pyplot.close(ax.figure)

                dR_Tau_Jet = result_JetMatching["dR_Tau_Jet"].integrate("dataset",dataset).integrate("Lxy",int_range=slice(*Lxy_slice))
                ax = coffea.hist.plot2d(dR_Tau_Jet, xaxis='dR_Tau_Jet', patch_opts={"norm":colors.LogNorm()})
                ax.set_title(f'dR(Tau, {cfg.collections.Jets.name}) ({dataset}) Lxy={Lxy_slice[0]}_{Lxy_slice[1]}')
                ax.figure.set_dpi(72)
                ax.figure.tight_layout()
                ax.figure.savefig(path_jet+f'/dR_Tau_Jet_STau_dataset_{dataset}_Lxy={Lxy_slice[0]}_{Lxy_slice[1]}.png')
                # ax.figure.savefig(path_jet+f'/dR_Tau_Jet_dataset_{dataset}_Lxy_{Lxy_slice[0]}-{Lxy_slice[1]}.pdf')
                matplotlib.pyplot.close(ax.figure)

                dR_STau_Track = result_JetMatching["dR_STau_Track"].integrate("dataset",dataset).integrate("Lxy",int_range=slice(*Lxy_slice))
                ax = coffea.hist.plot2d(dR_STau_Track, xaxis='dR_STau_lostTrack', patch_opts={"norm":colors.LogNorm()})
                ax.set_title(f'dR(Tau, lostTrack) ({dataset}) Lxy={Lxy_slice[0]}_{Lxy_slice[1]}')
                ax.figure.set_dpi(72)
                ax.figure.tight_layout()
                ax.figure.savefig(path_pfCand+f'/dR_STau_lostTrack_dataset_{dataset})_Lxy={Lxy_slice[0]}_{Lxy_slice[1]}.png')
                # ax.figure.savefig(path_pfCand+f'/dR_STau_lostTrack_dataset_{dataset})_Lxy_{Lxy_slice[0]}-{Lxy_slice[1]}.pdf')
                matplotlib.pyplot.close(ax.figure)

                dR_STau_pfCand = result_JetMatching["dR_STau_pfCand"].integrate("dataset",dataset).integrate("Lxy",int_range=slice(*Lxy_slice))
                ax = coffea.hist.plot2d(dR_STau_pfCand, xaxis='dR_STau_pfCand', patch_opts={"norm":colors.LogNorm()})
                ax.set_title(f'dR(Tau, pfCand) ({dataset}) Lxy={Lxy_slice[0]}_{Lxy_slice[1]}')
                ax.figure.set_dpi(72)
                ax.figure.tight_layout()
                ax.figure.savefig(path_pfCand+f'/dR_STau_pfCand_dataset_{dataset})_Lxy={Lxy_slice[0]}_{Lxy_slice[1]}.png')
                # ax.figure.savefig(path_pfCand+f'/dR_STau_pfCand_dataset_{dataset})_Lxy_{Lxy_slice[0]}-{Lxy_slice[1]}.pdf')
                matplotlib.pyplot.close(ax.figure)

                # Additional projection plots:
                dR_Tau_Jet_project = dR_Tau_Jet.integrate("dR_Tau_STau")
                dR_Tau_Jet_project = conv_to_root(dR_Tau_Jet_project, f'dR_Tau_Jet_project{dataset})_Lxy={Lxy_slice[0]}_{Lxy_slice[1]}')
                
                dR_Tau_Jet_project.SetLineColor(46)
                dR_Tau_Jet_project.SetLineWidth(2)
                dR_Tau_Jet_project.SetMarkerColor(3)
                dR_Tau_Jet_project.SetMarkerSize(0)
                dR_Tau_Jet_project.SetTitle("dR(#tau,Jet)")
                
                outfile = path_jet+f'/dR_Tau_Jet_{dataset}_Lxy_{Lxy_slice[0]}-{Lxy_slice[1]}.png'
                
                utils.utils.root_plot1D(
                    l_hist = [dR_Tau_Jet_project],
                    outfile = outfile,
                    xrange = [0, 0.5],
                    yrange = (0, dR_Tau_Jet_project.GetMaximum()+0.5*dR_Tau_Jet_project.GetMaximum()),
                    logx = False, logy = False,
                    ytitle = "Efficiency",
                    xtitle = "dR (#tau, jet)",
                    centertitlex = True, centertitley = True,
                    centerlabelx = False, centerlabely = False,
                    gridx = True, gridy = True,
                    ndivisionsx = None,
                    stackdrawopt = "nostack",
                    legendpos = "UL",
                    legendtitle = f"Lxy: {Lxy_slice[0]}-{Lxy_slice[1]} cm",
                    legendncol = 1,
                    #legendtextsize = 0.04,
                    legendwidthscale = 1.0,
                    legendheightscale = 1.0,
                    lumiText = "2018 (13 TeV)"
                )
                
                dR_STau_Jet_project = dR_STau_Jet.integrate("dR_STau_Tau")
                dR_STau_Jet_project = conv_to_root(dR_STau_Jet_project, f'dR_STau_Jet_project{dataset})_Lxy={Lxy_slice[0]}_{Lxy_slice[1]}')
                
                dR_STau_Jet_project.SetLineColor(46)
                dR_STau_Jet_project.SetLineWidth(2)
                dR_STau_Jet_project.SetMarkerColor(3)
                dR_STau_Jet_project.SetMarkerSize(0)
                dR_STau_Jet_project.SetTitle("dR(s#tau,Jet)")
                
                
                outfile = path_jet+f'/dR_STau_Jet_{dataset}__Lxy={Lxy_slice[0]}_{Lxy_slice[1]}.png'
                
                utils.utils.root_plot1D(
                    l_hist = [dR_STau_Jet_project],
                    outfile = outfile,
                    xrange = [0, 0.5],
                    yrange = (0, dR_STau_Jet_project.GetMaximum()+0.5*dR_STau_Jet_project.GetMaximum()),
                    logx = False, logy = False,
                    ytitle = "Efficiency",
                    xtitle = "dR (s#tau, jet)",
                    centertitlex = True, centertitley = True,
                    centerlabelx = False, centerlabely = False,
                    gridx = True, gridy = True,
                    ndivisionsx = None,
                    stackdrawopt = "nostack",
                    legendpos = "UL",
                    legendtitle = f"Lxy: {Lxy_slice[0]}-{Lxy_slice[1]} cm",
                    legendncol = 1,
                    #legendtextsize = 0.04,
                    legendwidthscale = 1.0,
                    legendheightscale = 1.0,
                    lumiText = "2018 (13 TeV)"
                )
                

    if 'efficiency_study' in cfg.mode:

        path_jet = cfg.output+"/eff_study/pt"
        os.makedirs(path_jet, exist_ok=True)

        for dataset in samples:
            
            for Lxy_slice in cfg.eff_setups.Lxy_slice:
                
                hist_gen_jet_pt = result_JetMatching["gen_jet_pt"].integrate("dataset",dataset).integrate("Lxy",int_range=slice(*Lxy_slice))
                hist_gen_jet_pt_reco = result_JetMatching["gen_jet_pt_reco"].integrate("dataset",dataset).integrate("Lxy",int_range=slice(*Lxy_slice))
                hist_gen_tau_pt = result_JetMatching["gen_tau_pt"].integrate("dataset",dataset).integrate("Lxy",int_range=slice(*Lxy_slice))
                hist_gen_tau_pt_reco = result_JetMatching["gen_tau_pt_reco"].integrate("dataset",dataset).integrate("Lxy",int_range=slice(*Lxy_slice))
                
                hist_gen_jet_pt = conv_to_root(hist_gen_jet_pt, f'gen_jet_pt_{dataset})_Lxy={Lxy_slice[0]}_{Lxy_slice[1]}')
                hist_gen_jet_pt_reco = conv_to_root(hist_gen_jet_pt_reco, f'gen_jet_pt_reco_{dataset})_Lxy={Lxy_slice[0]}_{Lxy_slice[1]}')
                hist_gen_tau_pt = conv_to_root(hist_gen_tau_pt, f'gen_tau_pt_{dataset})_Lxy={Lxy_slice[0]}_{Lxy_slice[1]}')
                hist_gen_tau_pt_reco = conv_to_root(hist_gen_tau_pt_reco, f'gen_tau_pt_reco_{dataset})_Lxy={Lxy_slice[0]}_{Lxy_slice[1]}')

                hist_gen_jet_pt_ratio = hist_gen_jet_pt.Clone()
                hist_gen_jet_pt_reco_ratio = hist_gen_jet_pt_reco.Clone()
                hist_gen_tau_pt_reco_ratio= hist_gen_tau_pt_reco.Clone()
                
                hist_gen_jet_pt_ratio.Divide(hist_gen_tau_pt)
                hist_gen_jet_pt_reco_ratio.Divide(hist_gen_tau_pt)
                hist_gen_tau_pt_reco_ratio.Divide(hist_gen_tau_pt)

                hist_gen_jet_pt_ratio.SetLineColor(9)
                hist_gen_jet_pt_ratio.SetLineWidth(2)
                hist_gen_jet_pt_ratio.SetMarkerColor(2)
                hist_gen_jet_pt_ratio.SetMarkerSize(0)
                hist_gen_jet_pt_ratio.SetTitle("genJet/genTau (dR<0.4)")

                hist_gen_jet_pt_reco_ratio.SetLineColor(8)
                hist_gen_jet_pt_reco_ratio.SetLineWidth(2)
                hist_gen_jet_pt_reco_ratio.SetMarkerColor(3)
                hist_gen_jet_pt_reco_ratio.SetMarkerSize(0)
                hist_gen_jet_pt_reco_ratio.SetTitle("recoJet/genTau (dR<0.4)")
                
                hist_gen_tau_pt_reco_ratio.SetLineColor(46)
                hist_gen_tau_pt_reco_ratio.SetLineWidth(2)
                hist_gen_tau_pt_reco_ratio.SetMarkerColor(3)
                hist_gen_tau_pt_reco_ratio.SetMarkerSize(0)
                hist_gen_tau_pt_reco_ratio.SetTitle("recoTau/genTau (dR<0.4)")
                
                outfile = path_jet+f'/ratio_pt_{dataset}_Lxy_{Lxy_slice[0]}-{Lxy_slice[1]}.png'
                
                utils.utils.root_plot1D(
                    # l_hist = [hist_gen_jet_pt_ratio, hist_gen_jet_pt_reco_ratio, hist_gen_tau_pt_reco_ratio],
                    l_hist = [hist_gen_jet_pt_ratio, hist_gen_jet_pt_reco_ratio ],
                    outfile = outfile,
                    xrange = [10, 1010],
                    yrange = (0, hist_gen_jet_pt_ratio.GetMaximum()+0.5*hist_gen_jet_pt_ratio.GetMaximum()),
                    logx = True, logy = False,
                    ytitle = "Efficiency",
                    xtitle = "p_{T}(#tau_{vis}) GeV",
                    centertitlex = True, centertitley = True,
                    centerlabelx = False, centerlabely = False,
                    gridx = True, gridy = True,
                    ndivisionsx = None,
                    stackdrawopt = "nostack",
                    legendpos = "UL",
                    legendtitle = f"Lxy: {Lxy_slice[0]}-{Lxy_slice[1]} cm",
                    legendncol = 1,
                    #legendtextsize = 0.04,
                    legendwidthscale = 1.0,
                    legendheightscale = 1.0,
                    lumiText = "2018 (13 TeV)"
                )

        path_jet = cfg.output+"/eff_study/Lxy"
        os.makedirs(path_jet, exist_ok=True)

        for dataset in samples:

            for dataset in samples:

                hist_gen_jet_pt = result_JetMatching["gen_jet_pt"].integrate("dataset",dataset).integrate("gen_tau_pt")
                hist_gen_jet_pt_reco = result_JetMatching["gen_jet_pt_reco"].integrate("dataset",dataset).integrate("gen_tau_pt")
                hist_gen_tau_pt = result_JetMatching["gen_tau_pt"].integrate("dataset",dataset).integrate("gen_tau_pt")
                hist_gen_tau_pt_reco = result_JetMatching["gen_tau_pt_reco"].integrate("dataset",dataset).integrate("gen_tau_pt")
                
                hist_gen_jet_pt = conv_to_root(hist_gen_jet_pt, f'gen_jet_Lxy')
                hist_gen_jet_pt.Rebin(1)
                hist_gen_jet_pt_reco = conv_to_root(hist_gen_jet_pt_reco, f'gen_jet_reco_Lxy')
                hist_gen_jet_pt_reco.Rebin(1)
                hist_gen_tau_pt = conv_to_root(hist_gen_tau_pt, f'gen_tau_Lxy')
                hist_gen_tau_pt.Rebin(1)
                hist_gen_tau_pt_reco = conv_to_root(hist_gen_tau_pt_reco, f'gen_tau_reco_Lxy')
                hist_gen_tau_pt_reco.Rebin(1)

                hist_gen_jet_pt_ratio = hist_gen_jet_pt.Clone()
                hist_gen_jet_pt_reco_ratio = hist_gen_jet_pt_reco.Clone()
                hist_gen_tau_pt_reco_ratio = hist_gen_tau_pt_reco.Clone()

                hist_gen_jet_pt_ratio.Divide(hist_gen_tau_pt)
                hist_gen_jet_pt_reco_ratio.Divide(hist_gen_tau_pt)
                hist_gen_tau_pt_reco_ratio.Divide(hist_gen_tau_pt)

                hist_gen_jet_pt_ratio.SetLineColor(9)
                hist_gen_jet_pt_ratio.SetLineWidth(2)
                hist_gen_jet_pt_ratio.SetMarkerColor(2)
                hist_gen_jet_pt_ratio.SetMarkerSize(0)
                hist_gen_jet_pt_ratio.SetTitle("genJet/genTau (dR<0.4)")

                hist_gen_jet_pt_reco_ratio.SetLineColor(8)
                hist_gen_jet_pt_reco_ratio.SetLineWidth(2)
                hist_gen_jet_pt_reco_ratio.SetMarkerColor(3)
                hist_gen_jet_pt_reco_ratio.SetMarkerSize(0)
                hist_gen_jet_pt_reco_ratio.SetTitle("recoJet/genTau (dR<0.4)")
                
                hist_gen_tau_pt_reco_ratio.SetLineColor(46)
                hist_gen_tau_pt_reco_ratio.SetLineWidth(2)
                hist_gen_tau_pt_reco_ratio.SetMarkerColor(3)
                hist_gen_tau_pt_reco_ratio.SetMarkerSize(0)
                hist_gen_tau_pt_reco_ratio.SetTitle("recoTau/genTau (dR<0.4)")
                
                outfile = path_jet+f'/ratio_Lxy_{dataset}.png'
                
                utils.utils.root_plot1D(
                    l_hist = [hist_gen_jet_pt_ratio, hist_gen_jet_pt_reco_ratio, hist_gen_tau_pt_reco_ratio],
                    outfile = outfile,
                    xrange = [0.0001, 110],
                    yrange = (0, hist_gen_jet_pt_ratio.GetMaximum()+0.5*hist_gen_jet_pt_ratio.GetMaximum()),
                    logx = True, logy = False,
                    ytitle = "Efficiency",
                    xtitle = "Lxy(#tau_{vis}) cm",
                    centertitlex = True, centertitley = True,
                    centerlabelx = False, centerlabely = False,
                    gridx = True, gridy = True,
                    ndivisionsx = None,
                    stackdrawopt = "nostack",
                    legendpos = "UL",
                    legendtitle = f"",
                    legendncol = 1,
                    #legendtextsize = 0.04,
                    legendwidthscale = 1.0,
                    legendheightscale = 1.0,
                    lumiText = "2018 (13 TeV)"
                )


        # path_jet = cfg.output+"/eff_study/ratio_E"
        # os.makedirs(path_jet, exist_ok=True)

        # for dataset in samples:
            
        #     for Lxy_slice in cfg.eff_setups.Lxy_slice:

        #         # gen_energy_ratio
        #         # energy_ratio

        #         hist_gen_jet_ratio = result_JetMatching["gen_energy_ratio"].integrate("dataset",dataset).integrate("Lxy",int_range=slice(*Lxy_slice))
        #         hist_tau_childr_energy_ratio = result_JetMatching["tau_childr_energy_ratio"].integrate("dataset",dataset).integrate("Lxy",int_range=slice(*Lxy_slice))
        #         hist_tau_childr_visTau_energy_ratio = result_JetMatching["tau_childr_visTau_energy_ratio"].integrate("dataset",dataset).integrate("Lxy",int_range=slice(*Lxy_slice))  

                
        #         hist_gen_jet_ratio = conv_to_root(hist_gen_jet_ratio, f'energy_ratio{dataset})_Lxy={Lxy_slice[0]}_{Lxy_slice[1]}')
        #         hist_tau_childr_energy_ratio = conv_to_root(hist_tau_childr_energy_ratio, f'children/mat_children_ratio{dataset})_Lxy={Lxy_slice[0]}_{Lxy_slice[1]}')
        #         hist_tau_childr_visTau_energy_ratio = conv_to_root(hist_tau_childr_visTau_energy_ratio, f'vis_tau/mat_children_ratio{dataset})_Lxy={Lxy_slice[0]}_{Lxy_slice[1]}')

        #         hist_gen_jet_ratio.SetLineColor(2)
        #         hist_gen_jet_ratio.SetLineWidth(2)
        #         hist_gen_jet_ratio.SetMarkerColor(2)
        #         hist_gen_jet_ratio.SetMarkerSize(0)
        #         hist_gen_jet_ratio.SetTitle("genJet_{E}/#tau^{vis}_{E} (dR<0.4)")
                
        #         hist_tau_childr_energy_ratio.SetLineColor(4)
        #         hist_tau_childr_energy_ratio.SetLineWidth(2)
        #         hist_tau_childr_energy_ratio.SetMarkerColor(4)
        #         hist_tau_childr_energy_ratio.SetMarkerSize(0)
        #         hist_tau_childr_energy_ratio.SetTitle("#sum(#tau_{child_E})/#sum(#tau_{child_E}, dR_{genJet}<0.4)")

        #         hist_tau_childr_visTau_energy_ratio.SetLineColor(6)
        #         hist_tau_childr_visTau_energy_ratio.SetLineWidth(2)
        #         hist_tau_childr_visTau_energy_ratio.SetMarkerColor(6)
        #         hist_tau_childr_visTau_energy_ratio.SetMarkerSize(0)
        #         hist_tau_childr_visTau_energy_ratio.SetTitle("#sum(#tau^{vis}_{E})/#sum(#tau_{child_E}, dR_{genJet}<0.4)")

        #         outfile = path_jet+f'/ratio_E_{dataset}_Lxy_{Lxy_slice[0]}-{Lxy_slice[1]}.png'
                
        #         utils.utils.root_plot1D(
        #             l_hist = [hist_gen_jet_ratio],
        #             outfile = outfile,
        #             xrange = [cfg.eff_setups.energy_ratio[1], cfg.eff_setups.energy_ratio[2]],
        #             yrange = (0, hist_gen_jet_ratio.GetMaximum()+0.05*hist_gen_jet_ratio.GetMaximum()),
        #             logx = False, logy = False,
        #             ytitle = "Freq.",
        #             xtitle = "E_{obj_1}/E_{obj_2}",
        #             centertitlex = True, centertitley = True,
        #             centerlabelx = False, centerlabely = False,
        #             gridx = True, gridy = True,
        #             ndivisionsx = None,
        #             stackdrawopt = "nostack",
        #             legendpos = "UR",
        #             legendtitle = f"Lxy=({Lxy_slice[0]}-{Lxy_slice[1]})",
        #             legendncol = 1,
        #             #legendtextsize = 0.04,
        #             legendwidthscale = 1.3,
        #             legendheightscale = 1.5,
        #             lumiText = "2018 (13 TeV)"
        #         )

        #         outfile2 = path_jet+f'/ratio_E_child_{dataset}_Lxy_{Lxy_slice[0]}-{Lxy_slice[1]}.png'

        #         utils.utils.root_plot1D(
        #             l_hist = [hist_tau_childr_energy_ratio, hist_tau_childr_visTau_energy_ratio],
        #             outfile = outfile2,
        #             xrange = [cfg.eff_setups.energy_ratio_children[1], cfg.eff_setups.energy_ratio_children[2]],
        #             yrange = (0, hist_tau_childr_energy_ratio.GetMaximum()+0.05*hist_tau_childr_energy_ratio.GetMaximum()),
        #             logx = False, logy = False,
        #             ytitle = "Freq.",
        #             xtitle = "E_{obj_1}/E_{obj_2}",
        #             centertitlex = True, centertitley = True,
        #             centerlabelx = False, centerlabely = False,
        #             gridx = True, gridy = True,
        #             ndivisionsx = None,
        #             stackdrawopt = "nostack",
        #             legendpos = "UL",
        #             legendtitle = f"Lxy=({Lxy_slice[0]}-{Lxy_slice[1]})",
        #             legendncol = 1,
        #             #legendtextsize = 0.04,
        #             legendwidthscale = 1.3,
        #             legendheightscale = 1.5,
        #             lumiText = "2018 (13 TeV)"
        #         )


        path_res = cfg.output+"/eff_study/resolution_pt"
        os.makedirs(path_res, exist_ok=True)

        for dataset in samples:
            
            for Lxy_slice in cfg.eff_setups.Lxy_slice:

                # hist_pt_ratio_genJet = result_JetMatching["pt_ratio_genJet"].integrate("dataset",dataset).integrate("Lxy",int_range=slice(*Lxy_slice))
                # ax = coffea.hist.plot2d(hist_pt_ratio_genJet, xaxis='pt_gen', patch_opts={"norm":colors.LogNorm()})
                # ax.set_title(f'pt_reco/pt_genJet ({dataset}) Lxy({Lxy_slice[0]}-{Lxy_slice[1]})')
                # ax.figure.set_dpi(72)
                # ax.figure.tight_layout()
                # ax.figure.savefig(path_res+f'/pt_ratio_genJet_{dataset})_Lxy_{Lxy_slice[0]}-{Lxy_slice[1]}.png')
                # matplotlib.pyplot.close(ax.figure)


                hist_pt_ratio_visTau = result_JetMatching["pt_ratio_visTau"].integrate("dataset",dataset).integrate("Lxy",int_range=slice(*Lxy_slice))
                ax = coffea.hist.plot2d(hist_pt_ratio_visTau, xaxis='pt_gen', patch_opts={"norm":colors.LogNorm()})
                ax.set_title(f'pt_reco/pt_visTau ({dataset}) Lxy({Lxy_slice[0]}-{Lxy_slice[1]})')
                ax.figure.set_dpi(72)
                ax.figure.tight_layout()
                ax.figure.savefig(path_res+f'/pt_ratio_visTau_{dataset})_Lxy={Lxy_slice[0]}_{Lxy_slice[1]}.png')
                matplotlib.pyplot.close(ax.figure)


                hist_visTau_resolution = result_JetMatching["pt_resolution_visTau"].integrate("dataset",dataset).integrate("pt_gen").integrate("Lxy",int_range=slice(*Lxy_slice))
                hist_visTau_resolution = conv_to_root(hist_visTau_resolution, f'visTau_resolution_pt_{dataset})_Lxy={Lxy_slice[0]}_{Lxy_slice[1]}')

                hist_visTau_resolution.SetLineColor(2)
                hist_visTau_resolution.SetLineWidth(2)
                hist_visTau_resolution.SetMarkerColor(2)
                hist_visTau_resolution.SetMarkerSize(0)
                hist_visTau_resolution.SetTitle("pt_resolution")

                outfile3 = path_res+f'/pt_resolution_visTau_{dataset}_Lxy={Lxy_slice[0]}_{Lxy_slice[1]}.png'

                utils.utils.root_plot1D(
                    l_hist = [hist_visTau_resolution],
                    outfile = outfile3,
                    xrange = [cfg.eff_setups.pt_resolution[1], cfg.eff_setups.pt_resolution[2]],
                    yrange = (1, hist_visTau_resolution.GetMaximum()+0.05*hist_visTau_resolution.GetMaximum()),
                    logx = False, logy = True,
                    ytitle = "arb. units",
                    xtitle = "pt(recoJet)-pt(#tau_{vis}) / pt(#tau_{vis})",
                    centertitlex = True, centertitley = True,
                    centerlabelx = False, centerlabely = False,
                    gridx = True, gridy = True,
                    ndivisionsx = None,
                    stackdrawopt = "nostack",
                    legendpos = "UR",
                    legendtitle = f"Lxy: {Lxy_slice[0]}-{Lxy_slice[1]} cm",
                    legendncol = 1,
                    #legendtextsize = 0.04,
                    legendwidthscale = 0.7,
                    legendheightscale = 0.7,
                    lumiText = "2018 (13 TeV)"
                )

if __name__ == "__main__":

    regionStudy()
    
