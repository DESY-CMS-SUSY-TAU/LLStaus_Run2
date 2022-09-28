# This file illustrates how to implement a processor, realizing the selection
# steps and outputting histograms and a cutflow with efficiencies.
# Here we create a very simplified version of the ttbar-to-dilep processor.
# One can run this processor using
# 'python3 -m pepper.runproc --debug example_processor.py example_config.json'
# Above command probably will need a little bit of time before all cuts are
# applied once. This is because a chunk of events are processed simultaneously.
# You change adjust the number of events in a chunk and thereby the memory
# usage by using the --chunksize parameter (the default value is 500000).

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

        # Pick up leading jet (by score)
        selector.set_column("jet_1st", partial(self.leading_jet, order=1)) #1st jet
        selector.set_column("jet_2nd", partial(self.leading_jet, order=2)) #2nd jet


    def leading_jet(self, data, order=0):
        jets = data["Jet"]

        # Good jets only
        is_good = (
            (self.config["jet_eta_min"] < jets.eta)
            & (jets.eta < self.config["jet_eta_max"])
            & (self.config["good_jet_pt_min"] < jets.pt))
        
        jets = jets[is_good]
        # We want to get leading jet by the score of displaced tau tagger
        idx_leading = \
            ak.argsort(jets.disTauTag_score1, ascending=False)[:,:order]

        return jets[idx_leading]