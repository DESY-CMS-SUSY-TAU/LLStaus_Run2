from nis import match
from unittest import result
import pepper
import awkward as ak
from functools import partial
import numpy as np
import mt2
import numba as nb

from coffea.nanoevents import NanoAODSchema
import utils.processor_stau_basic as stau_basic
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
        selector.add_cut("Trigger", partial(
            self.passing_trigger, pos_triggers, neg_triggers))
        
        # MET cut
        selector.add_cut("MET", self.MET_cut)
        
        # one tight tag muon
        selector.set_column("muon_tag", self.muons_tag) 
        selector.add_cut("single_muon_tag", self.single_muon_tag_cut)
        
        # veto all loose muons
        selector.set_column("muon_veto", self.muons_veto)
        selector.add_cut("loose_muon_veto", self.loose_muon_veto_cut)
        
        # veto all loose electrons
        selector.set_column("electron_veto", self.electron_veto)
        selector.add_cut("loose_electron_veto", self.loose_electron_veto_cut)
        
        # mt of muon and pt_miss
        selector.set_column("mt_muon", partial(self.mt, name="muon_tag"))
        selector.add_cut("mt_muon", self.mt_muon_cut)
        
        # add cuts and selections on the jets
        selector.set_column("valid_jets", self.jet_selection)
        selector.add_cut("n_loose_bjets", self.n_loose_bjets_cut)
        selector.add_cut("n_valid_jets", self.n_valid_jets_cut)

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
    def single_muon_tag_cut(self, data):
        return ak.num(data["muon_tag"])==1
    
    @zero_handler
    def loose_muon_veto_cut(self, data):
        return ak.num(data["muon_veto"])==0
    
    @zero_handler
    def loose_electron_veto_cut(self, data):
        return ak.num(data["electron_veto"])==0
    
    @zero_handler
    def mt_muon_cut(self, data):
        return ((data["mt_muon"] > self.config["mt_muon_min"]) & (data["mt_muon"] < self.config["mt_muon_max"]))
    
    @zero_handler
    def muons_tag(self, data):
        muons = data["Muon"]
        is_good = (
              (muons.pt > self.config["muon_pt_min"])
            & (muons.eta < self.config["muon_eta_max"])
            & (muons.eta > self.config["muon_eta_min"])
            & (muons[self.config["muon_ID"]] == 1)
            )
        return muons[is_good]
    
    @zero_handler
    def muons_veto(self, data):
        muons = data["Muon"]
        
        is_good = (
              (muons.pt > self.config["muon_veto_pt_min"])
            & (muons.eta < self.config["muon_veto_eta_max"])
            & (muons.eta > self.config["muon_veto_eta_min"])
            & (muons[self.config["muon_veto_ID"]] == 1)
            )
        
        is_tag =(
              (muons.pt > self.config["muon_pt_min"])
            & (muons.eta < self.config["muon_eta_max"])
            & (muons.eta > self.config["muon_eta_min"])
            & (muons[self.config["muon_ID"]] == 1)
            )
        
        return muons[(is_good & (~is_tag))]
    
    @zero_handler
    def electron_veto(self, data):
        ele = data["Electron"]
        # ele_low_eta_iso = ne.eval(self.config["elec_low_eta_iso"])
        # ele_high_eta_iso = ne.eval(self.config["elec_high_eta_iso"])
        ele_low_eta_iso =  ((np.abs(ele.eta) < 1.479) & (ele.pfRelIso03_all < (0.198+0.506/ele.pt)))
        ele_high_eta_iso = ((np.abs(ele.eta) > 1.479) & (np.abs(ele.eta) < 2.5) & (ele.pfRelIso03_all < (0.203+0.963/ele.pt)))
        isolation_cut = ( ele_low_eta_iso | ele_high_eta_iso )
        is_good = (
            isolation_cut
            & (ele.pt > self.config["elec_veto_pt"])
            & (ele.eta < self.config["elec_veto_eta_min"])
            & (ele.eta > self.config["elec_veto_eta_max"])
            & (ele[self.config["elec_Veto"]] == 1)
            )
        return ele[is_good]
    
    @zero_handler
    def jet_selection(self, data):
        jets = data["Jet"]
        jets = jets[(
            (self.config["jet_eta_min"] < jets.eta)
            & (jets.eta < self.config["jet_eta_max"])
            & (self.config["jet_pt_min"] < jets.pt)
            )]
        matches_h, dRlist = jets.nearest(data["muon_tag"], return_metric=True, threshold=self.config["tag_muon_veto_dR"])
        isoJets = jets[ak.is_none(matches_h, axis=-1)]
        return isoJets
    
    @zero_handler
    def n_loose_bjets_cut(self, data):
        return ak.num( data["valid_jets"][(data["valid_jets"].btagDeepFlavB > 0.2783)] ) == 0
    
    @zero_handler
    def n_valid_jets_cut(self,data):
        return ak.num(data["valid_jets"]) < 3
    
    @zero_handler
    def jets_updated_all(self, data):
        jets = data["valid_jets"]
        jets["dxy"] = data["jets_lead_pfcands"].dxy
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
