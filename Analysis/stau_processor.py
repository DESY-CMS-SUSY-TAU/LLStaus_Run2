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

from coffea.nanoevents import NanoAODSchema
# np.set_printoptions(threshold=np.inf)


# All processors should inherit from pepper.ProcessorBasicPhysics
class Processor(pepper.ProcessorBasicPhysics):
    # We use the ConfigTTbarLL instead of its base Config, to use some of its
    # predefined extras
    config_class = pepper.ConfigTTbarLL
    # schema_class = NanoAODSchema
    # schema_class.mixins.update({
    #     "PFCandidate": "PtEtaPhiMCollection",
    #     "SV": "PtEtaPhiMCollection",
    # })

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

    def zero_handler(func):
        def _function(self, data, *args, **kwargs):
            if len(data) > 0: return func(self, data, *args, **kwargs)
            else: return ak.Array([])
        return _function

    def process_selection(self, selector, dsname, is_mc, filler):
        # Implement the selection steps: add cuts, define objects and/or
        # compute event weights

        # Add a cut only allowing events according to the golden JSON
        # The good_lumimask method is specified in pepper.ProcessorBasicPhysics
        # It also requires a lumimask to be specified in config
        if not is_mc:
            selector.add_cut("Lumi", partial(
                self.good_lumimask, is_mc, dsname))

        # Only allow events that pass triggers specified in config
        # This also takes into account a trigger order to avoid triggering
        # the same event if it's in two different data datasets.
        pos_triggers, neg_triggers = pepper.misc.get_trigger_paths_for(
            dsname, is_mc, self.config["dataset_trigger_map"],
            self.config["dataset_trigger_order"])
        selector.add_cut("Trigger", partial(
            self.passing_trigger, pos_triggers, neg_triggers))

        # Select valid objects
        selector.set_column("muons_valid", self.muons_valid)
        selector.set_column("electrons_valid", self.electrons_valid)
        selector.set_column("jets_valid", self.jets_valid)
        selector.set_column("hps_taus_valid", self.hps_taus_valid)

        # Add basics cuts (at least 2 good jets)
        selector.add_cut("is_two_valid_jets", self.has_jets)

        # Pick up leading jet (by score)
        selector.set_column("jet_1", partial(self.leading_jet, order=0)) #1-st jet
        selector.set_column("jet_2", partial(self.leading_jet, order=1)) #2-nd jet
        selector.set_column("jet_b", self.b_tagged_jet)
        selector.set_column("sum_2jets", self.add_j1_j2)
        selector.set_column("distautag_double", self.distautag_double)
        
        selector.set_column("mt2_j1_j2_MET", partial(self.get_mt2, name_1 = "jet_1", name_2 = "jet_2"))

        selector.add_cut("b_tagged_cut", self.b_tagged_cut)
        selector.add_cut("MET_cut", self.MET_cut)

        # from jet1/2 we select only the objects that match / not match to tau
        selector.set_column("jet_1_gtau", partial(self.match_nearest, coll1="jet_1", coll2="GenVisTau", dR=0.4))
        selector.set_column("jet_2_gtau", partial(self.match_nearest, coll1="jet_2", coll2="GenVisTau", dR=0.4))

        selector.set_column("jet_1_!gtau", partial(self.match_nearest, coll1="jet_1", coll2="GenVisTau", dR=0.7, not_matched=True))
        selector.set_column("jet_2_!gtau", partial(self.match_nearest, coll1="jet_2", coll2="GenVisTau", dR=0.7, not_matched=True))

        selector.set_column("jet_1_hpstau", partial(self.match_nearest, coll1="jet_1", coll2="hps_taus_valid", dR=0.4))
        selector.set_column("jet_2_hpstau", partial(self.match_nearest, coll1="jet_2", coll2="hps_taus_valid", dR=0.4))

        selector.set_column("jet_1_!hpstau", partial(self.match_nearest, coll1="jet_1", coll2="hps_taus_valid", dR=0.7, not_matched=True))
        selector.set_column("jet_2_!hpstau", partial(self.match_nearest, coll1="jet_2", coll2="hps_taus_valid", dR=0.7, not_matched=True))

        selector.set_column("jet_1_gtau_hpstau", partial(self.match_nearest, coll1="jet_1_gtau", coll2="hps_taus_valid", dR=0.4))
        selector.set_column("jet_2_gtau_hpstau", partial(self.match_nearest, coll1="jet_2_gtau", coll2="hps_taus_valid", dR=0.4))
        
        selector.set_column("jet_1_!gtau_hpstau", partial(self.match_nearest, coll1="jet_1_!gtau", coll2="hps_taus_valid", dR=0.4))
        selector.set_column("jet_2_!gtau_hpstau", partial(self.match_nearest, coll1="jet_2_!gtau", coll2="hps_taus_valid", dR=0.4))

        selector.set_column("jet_1_gtau_!hpstau", partial(self.match_nearest, coll1="jet_1_gtau", coll2="hps_taus_valid", dR=0.4, not_matched=True))
        selector.set_column("jet_2_gtau_!hpstau", partial(self.match_nearest, coll1="jet_2_gtau", coll2="hps_taus_valid", dR=0.4, not_matched=True))
        
        selector.set_column("jet_1_!gtau_!hpstau", partial(self.match_nearest, coll1="jet_1_!gtau", coll2="hps_taus_valid", dR=0.4, not_matched=True))
        selector.set_column("jet_2_!gtau_!hpstau", partial(self.match_nearest, coll1="jet_2_!gtau", coll2="hps_taus_valid", dR=0.4, not_matched=True))

        # Pick up matched SVs
        selector.set_column("SV_jet1", partial(self.match_nearest, coll1="SV", coll2="jet_1", dR=0.4))
        selector.set_column("SV_jet2", partial(self.match_nearest, coll1="SV", coll2="jet_2", dR=0.4))

        # Pick up matched PfCand
        self.pfCand_Sequence(selector, input_jet_n = 1) # jet_1
        self.pfCand_Sequence(selector, input_jet_n = 2) # jet_2

        # study regions:
        # PP - 2 OS prompt tau_h
        # PD - 1 prompt tau_h, 1 displaced tau_h
        # DD - 2 displaced tau_h
        # selector.set_cat("categories", {"PP", "PD", "DD", "ALL"})
        # selector.set_multiple_columns(partial(self.category_masks))

    def pfCand_Sequence(self, selector, input_jet_n=None):
        if not input_jet_n:
            raise Exception("Error: input_jet_name should not be None")
        
        jet_n = str(input_jet_n)
        selector.set_column("pfCands_jet"+jet_n, partial(self.get_matched_pfCands, match_object="jet_"+jet_n, dR=0.4))
        selector.set_column("pfCands_jet"+jet_n+"_gtau", partial(self.get_matched_pfCands, match_object="jet_"+jet_n+"_gtau", dR=0.4))
        selector.set_column("pfCands_jet"+jet_n+"_!gtau", partial(self.get_matched_pfCands, match_object="jet_"+jet_n+"_!gtau", dR=0.4))

        selector.set_column("gen_stau",    self.gen_stau)
        selector.set_column("gen_stau_tau", self.gen_stau_tau)

        selector.set_column("jet_"+jet_n+"_!gtau_!allgtau", partial(self.match_nearest, coll1="jet_"+jet_n+"_!gtau", coll2="gen_stau_tau", dR=0.7, not_matched=True))
        selector.set_column("jet_"+jet_n+"_!gtau_!allgtau_!stau", partial(self.match_nearest, coll1="jet_"+jet_n+"_!gtau_!allgtau", coll2="gen_stau", dR=0.7, not_matched=True))

        # Choose PF Cands that are daughters of isolated jets
        selector.set_column("pfCands_jet"+jet_n+"_!gtau_!allgtau", partial(self.get_matched_pfCands, match_object="jet_"+jet_n+"_!gtau_!allgtau", dR=0.4))
        selector.set_column("pfCands_jet"+jet_n+"_!gtau_!allgtau_!stau", partial(self.get_matched_pfCands, match_object="jet_"+jet_n+"_!gtau_!allgtau_!stau", dR=0.4))

        # Choose PF Cands that match to HPS or not match to HPS
        selector.set_column("pfCands_jet"+jet_n+"_hpstau", partial(self.get_matched_pfCands, match_object="jet_"+jet_n+"_hpstau", dR=0.4))
        selector.set_column("pfCands_jet"+jet_n+"_!hpstau", partial(self.get_matched_pfCands, match_object="jet_"+jet_n+"_!hpstau", dR=0.4))

        selector.set_column("pfCands_jet"+jet_n+"_gtau_hpstau", partial(self.get_matched_pfCands, match_object="jet_"+jet_n+"_gtau_hpstau", dR=0.4))
        selector.set_column("pfCands_jet"+jet_n+"_gtau_!hpstau", partial(self.get_matched_pfCands, match_object="jet_"+jet_n+"_gtau_!hpstau", dR=0.4))

    @zero_handler
    def jets_valid(self, data):
        jets = data["Jet"]

        # study kinem region
        jets = jets[(
            (self.config["jet_eta_min"] < jets.eta)
            & (jets.eta < self.config["jet_eta_max"])
            & (self.config["jet_pt_min"] < jets.pt)
            )]

        # print(jets)
        print(data.metadata["dataset"])
        print(data.metadata['filename'])
        print(data["muons_valid"])
        print(data["electrons_valid"])
        print("mu:", ak.sum(ak.num(data["muons_valid"],axis=-1)), "per", ak.num(data["muons_valid"],axis=0),
              "e:", ak.sum(ak.num(data["electrons_valid"],axis=-1)),  "per", ak.num(data["electrons_valid"], axis=0))

        # muon veto
        jets_vetoed = jets.nearest(data["muons_valid"], return_metric=True, threshold=0.4)[0]
        print("1st veto", jets_vetoed)
        idx_vetoed_jets = ak.is_none(jets_vetoed, axis=-1)
        print(idx_vetoed_jets)
        jets = jets[idx_vetoed_jets]

        # electron veto
        jets_vetoed = jets.nearest(data["electrons_valid"], return_metric=True, threshold=0.4)[0]
        print("2nd veto", jets_vetoed)
        idx_vetoed_jets = ak.is_none(jets_vetoed, axis=-1)
        print(idx_vetoed_jets)
        jets = jets[idx_vetoed_jets]

        return jets

    @zero_handler
    def muons_valid(self, data):
        muons = data["Muon"]
        is_good = (
            (muons.pt > 20)
            & (muons.eta < 2.5)
            & (-2.5 < muons.eta)
            # & (muons.tightId == 1)
            & (muons.pfRelIso04_all<0.2)
            )
        return muons[is_good]

    @zero_handler
    def electrons_valid(self, data):
        ele = data["Electron"]
        low_eta_iso = ((np.abs(ele.eta) < 1.479) & (ele.pfRelIso03_all < (0.198+0.506/ele.pt)))
        high_eta_iso = ((np.abs(ele.eta) > 1.479) & (np.abs(ele.eta) < 2.5) & (ele.pfRelIso03_all < (0.203+0.963/ele.pt)))
        isolation_cut = ( low_eta_iso | high_eta_iso )
        is_good = (
            isolation_cut
            & (ele.pt > 20)
            & (ele.eta < 2.5)
            & (-2.5 < ele.eta)
            & (ele.convVeto == 1)
            # & (ele.mvaFall17V2Iso_WP80==1)
            )
        return ele[is_good]

    @zero_handler
    def hps_taus_valid(self, data):
        taus = data["Tau"]
        is_good = (
            (taus.pt > 20)
            & (taus.eta < 2.5)
            & (-2.5 < taus.eta)
            & (taus.idDeepTau2017v2p1VSe >= 1)
            & (taus.idDeepTau2017v2p1VSjet >= 1)
            & (taus.idDeepTau2017v2p1VSmu >=1)
            )
        return taus[is_good]

    @zero_handler
    def has_jets(self, data):
        return ak.num(data["jets_valid"]) >= 2

    @zero_handler
    def leading_jet(self, data, order=0):
        jets = data["jets_valid"]
        # We want to get leading jet by the score of displaced tau tagger
        # idx_leading = \
        #     ak.argsort(jets.disTauTag_score1, ascending=False)[:,order:order+1]
        idx_leading = \
            ak.argsort(jets.pt, ascending=False)[:,order:order+1]
        jets = jets[idx_leading]
        # print(jets.disTauTag_score1)
        # print(ak.is_none(jets.disTauTag_score1,axis=1))
        # print(ak.any(ak.is_none(jets.disTauTag_score1)))
        return jets

    @zero_handler
    def add_j1_j2(self, data):
        return data["jet_1"].add(data["jet_2"])

    @zero_handler
    def b_tagged_jet(self, data):
        # Jet_btagDeepFlavB satisfies the Medium (>0.2783) WP:
        # https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL18
        is_b = (data["Jet"].btagDeepFlavB > 0.2783)
        return data["Jet"][is_b]

    @zero_handler
    def b_tagged_cut(self, data):
        return ak.num(data["jet_b"]) <= 1

    @zero_handler    
    def MET_cut(self, data):
        return data["MET"].pt > 100 #GeV

    @zero_handler
    def match_nearest(self, data, coll1=None, coll2=None, dR = 0.4, not_matched = False):
        obj1 = data[coll1]
        obj2 = data[coll2]
        matches, dRlist = obj1.nearest(obj2, return_metric=True, threshold=dR)
        if not_matched:
            idx_matches_jets = ak.is_none(matches, axis=1)
        else:
            idx_matches_jets = ~ak.is_none(matches, axis=1)
        # results = obj1.mask[idx_matches_jets] # this function convert False to None e.g: [[obj, obj, None, None,..],..]
        results = obj1[idx_matches_jets] # this function drops all False e.g: [[obj, obj,..],..]
        sort_idx = ak.argsort(results.pt, axis=-1, ascending=False)
        return results[sort_idx]
                      

    @zero_handler    
    def get_mt2(self, data, name_1, name_2, name_MET = "MET") :
        
        return mt2.mt2(
            data[name_1].mass, data[name_1].px, data[name_1].py,
            data[name_2].mass, data[name_2].px, data[name_2].py,
            data[name_MET].px, data[name_MET].py,
            0, 0
        )

    @zero_handler
    def distautag_double(self, data):
        return data["jet_1"].disTauTag_score1*data["jet_2"].disTauTag_score1

    @zero_handler
    def get_matched_pfCands(self, data, match_object, dR=0.4):
        pfCands = self.match_nearest(data, coll1="PFCandidate", coll2=match_object, dR=dR)
        is_good = (
            (pfCands.pt > 5)
            & (pfCands.eta < 2.5)
            & (-2.5 < pfCands.eta)
            & (pfCands.hasTrackDetails)
        )

        # Add dxySig variable
        pfCands = pfCands[is_good]
        pfCands["dxySig"] = pfCands.dxy / pfCands.dxyError
        pfCands["Lrel"] = np.sqrt(pfCands.dxy**2 + pfCands.dz**2)
        return pfCands
    
    @zero_handler
    def gen_stau(self, data):
        return data.GenPart[ (abs(data.GenPart.pdgId) == 1000015) & 
                             (data.GenPart.hasFlags(["isLastCopy"])) ]
    
    @zero_handler
    def gen_stau_tau(self, data):
        # taus = data.gen_stau.children
        taus = data.GenPart
        taus = taus[ (abs(taus.pdgId) == 15) ]
        # return taus[:,:,0]
        return taus

    