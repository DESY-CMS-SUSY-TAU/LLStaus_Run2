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
import pepper
import awkward as ak
from functools import partial


# All processors should inherit from pepper.ProcessorBasicPhysics
class Processor(pepper.ProcessorBasicPhysics):
    # We use the ConfigTTbarLL instead of its base Config, to use some of its
    # predefined extras
    config_class = pepper.ConfigTTbarLL

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

        # Select jets that pass selection cuts 
        selector.set_column("jets_valid", self.jets_valid)
        # Add basics cuts (at least 2 good jets)
        selector.add_cut("At least 2 valid jets", self.has_jets)

        # Pick up leading jet (by score)
        selector.set_column("jet_1", partial(self.leading_jet, order=0)) #1-st jet
        selector.set_column("jet_2", partial(self.leading_jet, order=1)) #2-nd jet
        selector.set_column("jet_b", self.b_tagged_jet)
        selector.set_column("sum_2jets", self.add_j1_j2)

        # Pick up matched SVs
        selector.set_column("SV_1", partial(self.match_obj, name_1="jet_1", name_2="SV", dR=0.4))
        selector.set_column("SV_2", partial(self.match_obj, name_1="jet_2", name_2="SV", dR=0.4))
        
        selector.set_column("mt2_j1_j2_MET", partial(self.get_mt2, name_1 = "jet_1", name_2 = "jet_2"))

    def jets_valid(self, data):
        jets = data["Jet"]
        # Good jets only
        is_good = (
            (self.config["jet_eta_min"] < jets.eta)
            & (jets.eta < self.config["jet_eta_max"])
            & (self.config["jet_pt_min"] < jets.pt))
        return jets[is_good]

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

    def match_obj(self, data, name_1=None, name_2=None, dR = 0.4):
        obj1 = data[name_1][:,0]
        obj2 = data[name_2]
        _dR = obj1.delta_r(obj2)
        return obj2[(_dR<dR)]
    
    def get_mt2(self, data, name_1, name_2, name_MET = "MET") :
        
        return mt2.mt2(
            data[name_1].mass, data[name_1].px, data[name_1].py,
            data[name_2].mass, data[name_2].px, data[name_2].py,
            data[name_MET].px, data[name_MET].py,
            0, 0
        )
