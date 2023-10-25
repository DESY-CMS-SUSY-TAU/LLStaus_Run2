# This file illustrates how to implement a processor, realizing the selection
# steps and outputting histograms and a cutflow with efficiencies.
# Here we create a very simplified version of the ttbar-to-dilep processor.
# One can run this processor using
# 'python3 -m pepper.runproc --debug example_processor.py example_config.json'
# Above command probably will need a little bit of time before all cuts are
# applied once. This is because a chunk of events are processed simultaneously.
# You change adjust the number of events in a chunk and thereby the memory
# usage by using the --chunksize parameter (the default value is 500000).

from nis import match
from unittest import result
import pepper
import awkward as ak
from functools import partial
import numpy as np
import mt2
import numba as nb

from coffea.nanoevents import NanoAODSchema
# import utils.processor_stau as processor_stau
# np.set_printoptions(threshold=np.inf)


class Processor(pepper.ProcessorSTau):
    # We use the ConfigTTbarLL instead of its base Config, to use some of its
    # predefined extras
    # config_class = pepper.ConfigTTbarLL
    
    config_class = pepper.ConfigBasicPhysics
    
    def zero_handler(func):
        def _function(self, data, *args, **kwargs):
            if len(data) > 0: return func(self, data, *args, **kwargs)
            else: return ak.Array([])
        return _function

    def __init__(self, config, eventdir):
        # Initialize the class, maybe overwrite some config variables and
        # load additional files if needed
        # Can set and modify configuration here as well
        config["histogram_format"] = "root"
        # Need to call parent init to make histograms and such ready
        super().__init__(config, eventdir)
        
        if "jet_fake_rate" in self.config:
            self.jet_fake_rate = self.config["jet_fake_rate"]
        else:
            logger.warning("No jet fake rates specified")
            self.jet_fake_rate = []


        # It is not recommended to put anything as member variable into a
        # a Processor because the Processor instance is sent as raw bytes
        # between nodes when running on HTCondor.

    def process_selection(self, selector, dsname, is_mc, filler):
     
        # if not is_mc:
        #     selector.add_cut("Lumi", partial(
        #         self.good_lumimask, is_mc, dsname))

        pos_triggers, neg_triggers = pepper.misc.get_trigger_paths_for(
            dsname, is_mc, self.config["dataset_trigger_map"],
            self.config["dataset_trigger_order"])
        selector.add_cut("Trigger", partial(
            self.passing_trigger, pos_triggers, neg_triggers))

          # Select valid objects
        selector.set_column("muons_valid", self.muons_valid)
        selector.set_column("electrons_valid", self.electrons_valid)
        selector.set_column("hps_taus_valid", self.hps_taus_valid)
        selector.set_column("jets_valid", self.jets_valid)
        # selector.set_column("pfcand_valid", self.pfcand_valid)
        # selector.set_column("jets_valid_pf", partial(self.get_matched_pfCands, match_object="jets_valid", dR=0.4))

        selector.add_cut("MET_cut", self.MET_cut)
        selector.add_cut("muon_veto", self.muon_veto)
        selector.add_cut("elec_veto", self.elec_veto)
        
        selector.set_column("jet_b", self.b_tagged_jet)
        # selector.add_cut("b_tagged_0_cut", self.b_tagged_tight_cut)
        # selector.add_cut("is_two_valid_jets", self.has_jets)
        
        # selector.set_column("n_distautag_loose", partial(self.n_dis_tagged, threshold = 0.7785))
        # selector.set_column("n_distautag_medium", partial(self.n_dis_tagged, threshold = 0.9687))
        
        # signal histogram:
        selector.set_column("n_distautag_tight", partial(self.n_dis_tagged, threshold = 0.9972))
        
        selector.add_cut("apply_jet_fake_sf", self.apply_jet_fake_sf)
    
    @zero_handler
    def apply_jet_fake_sf(self, data):
        
        # Current scale factor yield does not consider fake rate dependancy over dxy of the jet

        jet_pt = data["jets_valid"].pt
        jet_eta = data["jets_valid"].eta
        # jets:    [ [j1, j2, j3, ...], [j1, j2, ...] ]
        # print(jet_pt)
        # print(jet_eta )
        
        fake_sf =  self.jet_fake_rate(pt=jet_pt, eta=jet_eta)
        # print(fake_sf)
        epsilon = fake_sf / (1 - fake_sf)
        # epsilon: [ [e1, e2, e3, ...], [e1, e2, ...] ]
        # print(epsilon)
        
        ###############
        # for 0-tight #
        ###############
        combinations = ak.combinations(epsilon, 2, axis=1)
        # combinations: [ [(e1, e2), (e2, e3), ...], [(e1, e2), ...] 
        # print(combinations)
        combinations_unzipped = ak.unzip(combinations)
        # unzip: [ [e1, e2, e3, ...], [e1, e2, ...], ...] <->  [ [e2, e3, e1, ...], [e2, e3, ...], ...]
        # print(combinations_unzipped)
        products = combinations_unzipped[0] * combinations_unzipped[1]
        # print(products)
        weight_0 = ak.sum(products, axis=1)
        # print(weight_0)
        # weight_0 per event: [ w1, w2, w3, .. ]
        
        ###############
        # for 1-tight #
        ###############
        epsilon_nontagged = epsilon[(data["jets_valid"].disTauTag_score1 < 0.9972)]
        # print(epsilon_nontagged)
        # epsilon: [ [e1, e2, e3, ...], [e1, e2, ...] ] *without tightly tagged jet"
        weight_1 = ak.sum(epsilon_nontagged, axis=1)
        # print(weight_1)
        weights = ak.where(data["n_distautag_tight"]==0, weight_0, 
                        ak.where(data["n_distautag_tight"]==1, weight_1, 1.0)
                    )
        # print(weights)
        return weights
    
    # @zero_handler
    # def apply_jet_fake_sf(self, data):
        
    #     jet_pt = data["jets_valid"].pt
    #     jet_eta = data["jets_valid"].eta
        
    #     fake_sf =  self.jet_fake_rate(pt=jet_pt, eta=jet_eta)
    #     epsilon = fake_sf / (1 - fake_sf)
    
    #     ###############
    #     # for 0-tight 
    #     ###############
    #     weight_1 = ak.sum(epsilon, axis=1)
    #     weights = ak.where(data["n_distautag_tight"]==0, weight_1, 1.0)
    #     return weights
    
    