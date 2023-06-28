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
import utils.processor_stau as processor_stau
# np.set_printoptions(threshold=np.inf)


class Processor(processor_stau.ProcessorSTau):
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
        selector.set_column("pfcand_valid", self.pfcand_valid)

        selector.add_cut("MET_cut", self.MET_cut)
        selector.add_cut("muon_veto", self.muon_veto)
        selector.add_cut("elec_veto", self.elec_veto)
        
        selector.add_cut("is_two_valid_jets", self.has_jets)
        
        selector.set_column("n_distautag_loose", partial(self.n_dis_tagged, threshold = 0.7785))
        selector.set_column("n_distautag_medium", partial(self.n_dis_tagged, threshold = 0.9687))
        selector.set_column("n_distautag_tight", partial(self.n_dis_tagged, threshold = 0.9995))
        
        # GenMatch jets
        # if is_mc:
        #     selector.set_column("match_valid_jets", partial(self.gen_match_jets, jet_name="jets_valid"))
        #     selector.set_column("gentau_jets", self.gentau_jets)
            
        # Pick up leading jet (by score)
        selector.set_column("jet_1", partial(self.leading_jet, order=0)) #1-st jet
        selector.set_column("jet_2", partial(self.leading_jet, order=1)) #2-nd jet
        selector.set_column("jet_b", self.b_tagged_jet)
        selector.set_column("sum_2jets", self.add_j1_j2)
        selector.set_column("distautag_double", self.distautag_double)
        
        selector.set_column("mt2_j1_j2_MET", partial(self.get_mt2, name_1 = "jet_1", name_2 = "jet_2"))

        selector.set_column("HT_valid", self.jet_momentum_sum_valid)
        selector.set_column("HT_miss_valid", self.jet_momentum_sum_miss_valid)
        selector.set_column("HT_miss_pt_miss_valid", self.HT_miss_pt_miss_valid)

        selector.set_column("HT", self.jet_momentum_sum)
        selector.set_column("HT_miss", self.jet_momentum_sum_miss)
        selector.set_column("HT_miss_pt_miss", self.HT_miss_pt_miss)

        selector.set_column("mt_jet1", partial(self.mt, name="jet_1"))
        selector.set_column("mt_jet2", partial(self.mt, name="jet_2"))
        selector.set_column("mt_sum", self.mt_sum)
        selector.set_column("dphi_jet1_jet2", self.dphi_jet1_jet2)
        
        selector.set_column("hpstau_1", partial(self.match_nearest, coll1="hps_taus_valid", coll2="jet_1", dR=0.4))
        selector.set_column("hpstau_2", partial(self.match_nearest, coll1="hps_taus_valid", coll2="jet_2", dR=0.4))

        # Pick up matched SVs
        selector.set_column("SV_jet1", partial(self.match_nearest, coll1="SV", coll2="jet_1", dR=0.4))
        selector.set_column("SV_jet2", partial(self.match_nearest, coll1="SV", coll2="jet_2", dR=0.4))

        # Pick up matched PfCand
        self.pfCand_Sequence(selector, input_jet_n = 1) # jet_1
        self.pfCand_Sequence(selector, input_jet_n = 2) # jet_2
        
        # selector.set_cat("control_region", {"M1", "M2", "ALL"})
        # selector.set_multiple_columns(partial(self.control_region))
        
        selector.add_cut("b_tagged_1_cut", self.b_tagged_cut)   
        selector.add_cut("b_tagged_0_cut", self.b_tagged_tight_cut)
    