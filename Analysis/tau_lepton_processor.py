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
        selector.set_column("muons_v1", self.muons_v1)
        selector.set_column("electrons_v1", self.electrons_v1)

        # Select valid objects
        selector.set_column("muons_v2", self.muons_v2)
        selector.set_column("electrons_v2", self.electrons_v2)

        # Select valid objects
        selector.set_column("muons_v3", self.muons_v3)
        selector.set_column("electrons_v3", self.electrons_v3)

        selector.set_column("jets_valid", self.jets_valid)
        # selector.set_column("hps_taus_valid", self.hps_taus_valid)

        selector.set_column("tau_muons", self.tau_muon)
        selector.set_column("tau_elecs", self.tau_elec)

        selector.set_column("tau_muons_jet", partial(self.match_nearest, coll1="tau_muons", coll2="jets_valid", dR=0.4))
        selector.set_column("tau_elecs_jet", partial(self.match_nearest, coll1="tau_elecs", coll2="jets_valid", dR=0.4))

        selector.set_column("tau_muons_v1", partial(self.match_nearest, coll1="tau_muons", coll2="muons_v1", dR=0.4))
        selector.set_column("tau_muons_v2", partial(self.match_nearest, coll1="tau_muons", coll2="muons_v2", dR=0.4))
        selector.set_column("tau_muons_v3", partial(self.match_nearest, coll1="tau_muons", coll2="muons_v3", dR=0.4))

        selector.set_column("tau_elecs_v1", partial(self.match_nearest, coll1="tau_elecs", coll2="electrons_v1", dR=0.4))
        selector.set_column("tau_elecs_v2", partial(self.match_nearest, coll1="tau_elecs", coll2="electrons_v2", dR=0.4))
        selector.set_column("tau_elecs_v3", partial(self.match_nearest, coll1="tau_elecs", coll2="electrons_v3", dR=0.4))

        # Add basics cuts (at least 2 good jets)

        # Pick up leading jet (by score)
        # selector.set_column("jet_1", partial(self.leading_jet, order=0)) #1-st jet
        # selector.set_column("jet_2", partial(self.leading_jet, order=1)) #2-nd jet

        selector.add_cut("is_two_valid_jets", self.has_jets)
        selector.set_column("jet_b", self.b_tagged_jet)
        selector.add_cut("b_tagged_cut", self.b_tagged_cut)
        selector.add_cut("MET_cut", self.MET_cut)


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
        # print(data.metadata["dataset"])
        # print(data.metadata['filename'])
        # print(data["muons_valid"])
        # print(data["electrons_valid"])
        # print("mu:", ak.sum(ak.num(data["muons_valid"],axis=-1)), "per", ak.num(data["muons_valid"],axis=0),
        #       "e:", ak.sum(ak.num(data["electrons_valid"],axis=-1)),  "per", ak.num(data["electrons_valid"], axis=0))

        # muon veto
        # jets_vetoed = jets.nearest(data["muons_v3"], return_metric=True, threshold=0.4)[0]
        # # print("1st veto", jets_vetoed)
        # idx_vetoed_jets = ak.is_none(jets_vetoed, axis=-1)
        # # print(idx_vetoed_jets)
        # jets = jets[idx_vetoed_jets]

        # # electron veto
        # jets_vetoed = jets.nearest(data["electrons_v3"], return_metric=True, threshold=0.4)[0]
        # # print("2nd veto", jets_vetoed)
        # idx_vetoed_jets = ak.is_none(jets_vetoed, axis=-1)
        # # print(idx_vetoed_jets)
        # jets = jets[idx_vetoed_jets]

        return jets

    @zero_handler
    def muons_v1(self, data):
        muons = data["Muon"]
        is_good = (
            (muons.pt > 20)
            & (muons.eta < 2.5)
            & (-2.5 < muons.eta)
            )
        return muons[is_good]
    
    @zero_handler
    def muons_v2(self, data):
        muons = data["Muon"]
        is_good = (
            (muons.pt > 20)
            & (muons.eta < 2.5)
            & (-2.5 < muons.eta)
            & (muons.pfRelIso04_all<0.2)
            )
        return muons[is_good]
        
    @zero_handler
    def muons_v3(self, data):
        muons = data["Muon"]
        is_good = (
            (muons.pt > 20)
            & (muons.eta < 2.5)
            & (-2.5 < muons.eta)
            & (muons.pfRelIso04_all<0.2)
            & (muons.tightId == 1)
            )
        return muons[is_good]
    
    @zero_handler
    def electrons_v1(self, data):
        ele = data["Electron"]
        is_good = (
            (ele.pt > 20)
            & (ele.eta < 2.5)
            & (-2.5 < ele.eta)
            )
        return ele[is_good]

    @zero_handler
    def electrons_v2(self, data):
        ele = data["Electron"]
        is_good = (
            (ele.pt > 20)
            & (ele.eta < 2.5)
            & (-2.5 < ele.eta)
            & (ele.convVeto == 1)
            )
        return ele[is_good]

    @zero_handler
    def electrons_v3(self, data):
        ele = data["Electron"]
        low_eta_iso = ((np.abs(ele.eta) < 1.479) & (ele.pfRelIso03_all < (0.198+0.506/ele.pt)))
        high_eta_iso = ((np.abs(ele.eta) > 1.479) & (np.abs(ele.eta) < 2.5) & (ele.pfRelIso03_all < (0.203+0.963/ele.pt)))
        isolation_cut = ( low_eta_iso | high_eta_iso )
        is_good = (
            (ele.pt > 20)
            & (ele.eta < 2.5)
            & (-2.5 < ele.eta)
            & (ele.convVeto == 1)
            & isolation_cut
            )
        return ele[is_good]

    # @zero_handler
    # def electrons_valid(self, data):
    #     ele = data["Electron"]
    #     low_eta_iso = ((np.abs(ele.eta) < 1.479) & (ele.pfRelIso03_all < (0.198+0.506/ele.pt)))
    #     high_eta_iso = ((np.abs(ele.eta) > 1.479) & (np.abs(ele.eta) < 2.5) & (ele.pfRelIso03_all < (0.203+0.963/ele.pt)))
    #     isolation_cut = ( low_eta_iso | high_eta_iso )
    #     is_good = (
    #         (ele.pt > 20)
    #         & (ele.eta < 2.5)
    #         & (-2.5 < ele.eta)
    #         & (ele.convVeto == 1)
    #         & isolation_cut
    #         # & (ele.mvaFall17V2Iso_WP80==1)
    #         )
    #     return ele[is_good]

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
    def tau_muon(self, data):
        taus = data.GenPart[ (abs(data.GenPart.pdgId) == 15) & (data.GenPart.hasFlags(["isLastCopy"])) ]
        muons = taus.children
        muons = muons[ ( abs(muons.pdgId) == 13 ) & (muons.pt > 30) & (muons.eta < 2.4) & (-2.4 < muons.eta) ]
        muons = ak.firsts(muons,axis=2)
        muons = muons[ (~ak.is_none(muons, axis=1)) ]
        muons['Lxy'] = np.sqrt( muons.vertexX*muons.vertexX + muons.vertexY*muons.vertexY )
        # print(muons)
        return muons
    
    @zero_handler
    def tau_elec(self, data):
        taus = data.GenPart[ (abs(data.GenPart.pdgId) == 15) & (data.GenPart.hasFlags(["isLastCopy"])) ]
        elec = taus.children
        elec = elec[ ( abs(elec.pdgId) == 11 ) & (elec.pt > 30) & (elec.eta < 2.4) & (-2.4 < elec.eta)  ] # every tau may have only one lepton as child
        elec = ak.firsts(elec,axis=2)
        elec = elec[ (~ak.is_none(elec, axis=1)) ]
        elec['Lxy'] = np.sqrt( elec.vertexX*elec.vertexX + elec.vertexY*elec.vertexY )
        # print(elec)
        return elec