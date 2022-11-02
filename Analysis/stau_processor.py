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
        
        selector.add_cut("b_tagged_cut", self.b_tagged_cut)
        selector.add_cut("MET_cut", self.MET_cut)

        # from jet1/2 we select only the objects that match / not match to tau
        selector.set_column("jet_1_gtau", partial(self.match_nearest, coll1="jet_1", coll2="GenVisTau", dR=0.4))
        selector.set_column("jet_2_gtau", partial(self.match_nearest, coll1="jet_2", coll2="GenVisTau", dR=0.4))

        selector.set_column("jet_1_!gtau", partial(self.match_nearest, coll1="jet_1", coll2="GenVisTau", dR=0.4, not_matched=True))
        selector.set_column("jet_2_!gtau", partial(self.match_nearest, coll1="jet_2", coll2="GenVisTau", dR=0.4, not_matched=True))

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
        selector.set_column("pfCands_jet1", partial(self.match_nearest, coll1="PFCandidate", coll2="jet_1", dR=0.4))
        selector.set_column("pfCands_jet2", partial(self.match_nearest, coll1="PFCandidate", coll2="jet_2", dR=0.4))
        selector.set_column("pfCands_jet1_DxySig", partial(self.pfCandDxySig, name="pfCands_jet1"))
        selector.set_column("pfCands_jet2_DxySig", partial(self.pfCandDxySig, name="pfCands_jet2"))

        selector.set_column("pfCands_jet1_gtau", partial(self.match_nearest, coll1="PFCandidate", coll2="jet_1_gtau", dR=0.4))
        selector.set_column("pfCands_jet2_gtau", partial(self.match_nearest, coll1="PFCandidate", coll2="jet_2_gtau", dR=0.4))
        selector.set_column("pfCands_jet1_gtau_DxySig", partial(self.pfCandDxySig, name="pfCands_jet1_gtau"))
        selector.set_column("pfCands_jet2_gtau_DxySig", partial(self.pfCandDxySig, name="pfCands_jet2_gtau"))

        selector.set_column("pfCands_jet1_!gtau", partial(self.match_nearest, coll1="PFCandidate", coll2="jet_1_!gtau", dR=0.4))
        selector.set_column("pfCands_jet2_!gtau", partial(self.match_nearest, coll1="PFCandidate", coll2="jet_2_!gtau", dR=0.4))
        selector.set_column("pfCands_jet1_!gtau_DxySig", partial(self.pfCandDxySig, name="pfCands_jet1_!gtau"))
        selector.set_column("pfCands_jet2_!gtau_DxySig", partial(self.pfCandDxySig, name="pfCands_jet2_!gtau"))

        selector.set_column("mt2_j1_j2_MET", partial(self.get_mt2, name_1 = "jet_1", name_2 = "jet_2"))

        selector.set_column("distautag_double", self.distautag_double)

    def jets_valid(self, data):
        jets = data["Jet"]

        # study kinem region
        jets = jets[(
            (self.config["jet_eta_min"] < jets.eta)
            & (jets.eta < self.config["jet_eta_max"])
            & (self.config["jet_pt_min"] < jets.pt)
            )]

        # muon veto
        jets_vetoed = jets.nearest(data["muons_valid"], return_metric=True, threshold=0.4)[0]
        idx_vetoed_jets = ak.is_none(jets_vetoed, axis=-1)
        jets = jets[idx_vetoed_jets]

        # electron veto
        jets_vetoed = jets.nearest(data["electrons_valid"], return_metric=True, threshold=0.4)[0]
        idx_vetoed_jets = ak.is_none(jets_vetoed, axis=-1)
        jets = jets[idx_vetoed_jets]

        return jets

    def muons_valid(self, data):
        muons = data["Muon"]
        is_good = (
            (muons.pt > 20)
            & (muons.eta < 2.5)
            & (-2.5 < muons.eta)
            & (muons.tightId == 1)
            & (muons.pfRelIso04_all<0.2)
            )
        return muons[is_good]

    def electrons_valid(self, data):
        ele = data["Electron"]
        is_good = (
            (ele.pt > 20)
            & (ele.eta < 2.5)
            & (-2.5 < ele.eta)
            & (ele.convVeto == 1)
            & (ele.mvaFall17V2Iso_WP80==1)
            )
        return ele[is_good]

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

    def has_jets(self, data):
        return ak.num(data["jets_valid"]) >= 2

    def leading_jet(self, data, order=0):
        jets = data["jets_valid"]
        # We want to get leading jet by the score of displaced tau tagger
        idx_leading = \
            ak.argsort(jets.disTauTag_score1, ascending=False)[:,order:order+1]
        return jets[idx_leading]

    def add_j1_j2(self, data):
        return data["jet_1"].add(data["jet_2"])

    def b_tagged_jet(self, data):
        # Jet_btagDeepFlavB satisfies the Medium (>0.2783) WP:
        # https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL18
        is_b = (data["Jet"].btagDeepFlavB > 0.2783)
        return data["Jet"][is_b]

    def b_tagged_cut(self, data):
        return ak.num(data["jet_b"]) <= 1
    
    def MET_cut(self, data):
        return data["MET"].pt > 100 #GeV

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
                      

    # def match_obj(self, data, coll1=None, coll2=None, dR = 0.4,
    #               return_coll1 = False, not_match = False ):
    #     """Function takes names of two collections
    #     and returns matched objects

    #     Parameters:
    #     coll1 -- a string, the name of first collection
    #     coll2 -- a string, the name of second collection
    #     dR -- matching threhold
    #     return_coll1 -- bool, controls return behaviour
    #     not_match -- bool, returned matched to the object or not matched
        
    #     Returns: if return_coll1 == True - return coll1 if any of coll2 match
    #     if return_coll1 == False - return coll2 elements that match to coll1
    #     """
    #     obj1 = ak.firsts(data[coll1])
    #     assert ak.count(obj1, axis=None) == ak.count(data[coll1], axis=None)
    #     obj2 = data[coll2]
    #     _dR = obj1.delta_r(obj2)

    #     if return_coll1:
    #         idx = np.expand_dims(ak.any((_dR<dR), axis=-1), axis=-1)
    #         return data[coll1].mask[~idx] if not_match else data[coll1].mask[idx]
    #     else:
    #         return obj2[(_dR>dR)] if not_match else obj2[(_dR<dR)]

    
    def get_mt2(self, data, name_1, name_2, name_MET = "MET") :
        
        return mt2.mt2(
            data[name_1].mass, data[name_1].px, data[name_1].py,
            data[name_2].mass, data[name_2].px, data[name_2].py,
            data[name_MET].px, data[name_MET].py,
            0, 0
        )

    def distautag_double(self, data):
        return data["jet_1"].disTauTag_score1*data["jet_2"].disTauTag_score1

    def pfCandDxySig(self, data, name):
        pfCand = data[name] 
        return pfCand.dxy/pfCand.dxyError