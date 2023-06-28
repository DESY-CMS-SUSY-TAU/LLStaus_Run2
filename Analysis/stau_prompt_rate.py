from nis import match
from unittest import result
import pepper
import awkward as ak
from functools import partial
import numpy as np
import mt2
import numba as nb

from coffea.nanoevents import NanoAODSchema
# np.set_printoptions(threshold=np.inf)


class Processor(pepper.ProcessorSTau):
    # We use the ConfigTTbarLL instead of its base Config, to use some of its
    # predefined extras
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

        # It is not recommended to put anything as member variable into a
        # a Processor because the Processor instance is sent as raw bytes
        # between nodes when running on HTCondor.

    def process_selection(self, selector, dsname, is_mc, filler):
        
        # Triggers
        pos_triggers, neg_triggers = pepper.misc.get_trigger_paths_for(
            dsname, is_mc, self.config["dataset_trigger_map"],
            self.config["dataset_trigger_order"])
        # selector.add_cut("Trigger", partial(
        #     self.passing_trigger, pos_triggers, neg_triggers))
        
        # MET cut
        selector.add_cut("MET", self.MET_cut)
        
        
        # add cuts and selections on the jets
        selector.set_column("valid_jets", self.jet_selection)
        selector.set_column("n_jets", self.n_jets)
        # selector.add_cut("n_loose_bjets", self.n_loose_bjets_cut)
        # selector.add_cut("n_valid_jets", self.n_valid_jets_cut)

        # match pfcands for dxy
        selector.set_column("pfcand_valid", self.pfcand_valid)
        selector.set_column("jets_lead_pfcands", partial(self.get_matched_pfCands, match_object="valid_jets", dR=0.4))
        
        selector.set_column("jets_updated_all", self.jets_updated_all)
        selector.set_column("jets_updated_loose", self.jets_updated_loose)
        selector.set_column("jets_updated_tight", self.jets_updated_tight)
    
    @zero_handler
    def MET_cut(self, data):
        return data["MET"].pt > self.config["MET_pt"]

    @zero_handler
    def jet_selection(self, data):
        jets = data["Jet"]
        jets = jets[(
            (self.config["jet_eta_min"] < jets.eta)
            & (jets.eta < self.config["jet_eta_max"])
            & (self.config["jet_pt_min"] < jets.pt)
            )]
        # matches_h, dRlist = jets.nearest(data["muon_tag"], return_metric=True, threshold=self.config["tag_muon_veto_dR"])
        # isoJets = jets[ak.is_none(matches_h, axis=-1)]
        
        tau_vis = data.GenVisTau[ ((data.GenVisTau.pt > 30) & (abs(data.GenVisTau.eta) < 2.4) &
                                  (data.GenVisTau.parent.hasFlags(["fromHardProcess"])))
                                ]
        matches_h, dRlist = jets.nearest(tau_vis, return_metric=True, threshold=0.3)
        tau_jet = jets[~ak.is_none(matches_h, axis=1)]
        # tau_jet = jets
        return tau_jet
    
    @zero_handler
    def n_jets(self, data):
        return ak.num(data["valid_jets"])
    
    # @zero_handler
    # def n_loose_bjets_cut(self, data):
    #     return ak.num( data["valid_jets"][(data["valid_jets"].btagDeepFlavB > 0.2783)] ) == 0
    
    # @zero_handler
    # def n_valid_jets_cut(self,data):
    #     return (ak.num(data["valid_jets"]) < 3) & (ak.num(data["valid_jets"]) > 0)
    
    @zero_handler
    def jets_updated_all(self, data):
        jets = data["valid_jets"]
        jets["dxy"] = np.abs(data["jets_lead_pfcands"].dxy)
        jets["dz"] = np.abs(data["jets_lead_pfcands"].dz)
        return jets
    
    @zero_handler
    def jets_updated_loose(self, data):
        jets = data["jets_updated_all"]
        jets = jets[ (jets.disTauTag_score1 > self.config["tag_loose"]) ]
        return jets
    
    @zero_handler
    def jets_updated_tight(self, data):
        jets = data["jets_updated_all"]
        jets = jets[ (jets.disTauTag_score1 > self.config["tag_tight"]) ]
        return jets
