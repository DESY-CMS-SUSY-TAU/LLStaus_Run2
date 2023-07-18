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

# @nb.njit
# def pad_array_byindex(array_flat, array, idx_of_mother, idx_to_mother):
#     '''
#     The function pads an array by the indexes to coorespond to mother array
#     Input example:
#         array: [[Tau1, Tau2, Tau3], [Tau1, Tau2], [Tau3]] type: ak.Array
#         idx_of_mother: [[0, 1, 2, 4], [0, 1, 2, 3], [0, 1, 2]] type: ak.Array
#         idx_to_mother: [[0, -1, 2], [-1, 1], [0]]  type: ak.Array
#     Output:
#         [[Tau1, None, Tau3, None], [None, Tau2, None, None], [Tau3, None, None]]
#     '''
    
#     fill_idx = 0
#     for event_i in range(len(idx_of_mother)):
#         for i, idx_jet in enumerate(idx_of_mother[event_i]):
#             for j, idx_to_jet in enumerate(idx_to_mother[event_i]):
#                 if idx_to_jet != -1 and idx_to_jet == idx_jet:
#                     array_flat[fill_idx] = array[event_i][j]
#                 fill_idx += 1
#     return array_flat

# def wrap_padder(jets, taus):
#     array_flat = np.full(ak.count(jets["initial_idx"], axis=None), -999)   
#     idDeepTau2017v2p1VSe = pad_array_byindex(array_flat,
#                                              taus["idDeepTau2017v2p1VSe"],
#                                              jets["initial_idx"], 
#                                              taus.jetIdx)
#     idDeepTau2017v2p1VSe = ak.unflatten(idDeepTau2017v2p1VSe, ak.num(jets["initial_idx"]))
#     return idDeepTau2017v2p1VSe

# All processors should inherit from pepper.ProcessorBasicPhysics
class Processor(pepper.ProcessorSTau):
    # We use the ConfigTTbarLL instead of its base Config, to use some of its
    # predefined extras
    config_class = pepper.ConfigBasicPhysics

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
        
        # if not is_mc:
        #     selector.add_cut("Lumi", partial(
        #         self.good_lumimask, is_mc, dsname))

        # Due to the bug in the current DisTauTag production compain
        # the whole file should be dropped if Nan values are detected in the
        
        # DisTauTag score (better to analyse the whole file and drop it)
        # selector.add_cut("NanDrop", self.nan_drop)

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
        selector.set_column("hps_taus_valid", self.hps_taus_valid)
        selector.set_column("jets_valid", self.jets_valid)
        selector.set_column("pfcand_valid", self.pfcand_valid)

        selector.add_cut("MET_cut", self.MET_cut)
        selector.add_cut("is_two_valid_jets", self.has_jets)
        
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
        selector.set_column("HT", self.jet_momentum_sum)
        selector.set_column("HT_miss", self.jet_momentum_sum_miss)
        selector.set_column("HT_miss_pt_miss", self.HT_miss_pt_miss)
        selector.set_column("mt_jet1", partial(self.mt, name="jet_1"))
        selector.set_column("mt_jet2", partial(self.mt, name="jet_2"))
        selector.set_column("mt_sum", self.mt_sum)
        selector.set_column("dphi_jet1_jet2", self.dphi_jet1_jet2)

        selector.add_cut("b_tagged_1_cut", self.b_tagged_cut)   
        
        selector.set_cat("control_region", {"M1", "M2", "M3", "ALL"})
        selector.set_multiple_columns(partial(self.control_region))

        # from jet1/2 we select only the objects that match / not match to tau
        # selector.set_column("jet_1_gtau", partial(self.match_nearest, coll1="jet_1", coll2="GenVisTau", dR=0.4))
        # selector.set_column("jet_2_gtau", partial(self.match_nearest, coll1="jet_2", coll2="GenVisTau", dR=0.4))

        # selector.set_column("jet_1_!gtau", partial(self.match_nearest, coll1="jet_1", coll2="GenVisTau", dR=0.7, not_matched=True))
        # selector.set_column("jet_2_!gtau", partial(self.match_nearest, coll1="jet_2", coll2="GenVisTau", dR=0.7, not_matched=True))

        # ??????????????
        # selector.set_column("jet_1_hpstau", partial(self.match_nearest, coll1="jet_1", coll2="hps_taus_valid", dR=0.4))
        # selector.set_column("jet_2_hpstau", partial(self.match_nearest, coll1="jet_2", coll2="hps_taus_valid", dR=0.4))
        # selector.set_column("hpstau_1", partial(self.match_nearest, coll1="hps_taus_valid", coll2="jet_1_hpstau", dR=0.4))
        # selector.set_column("hpstau_2", partial(self.match_nearest, coll1="hps_taus_valid", coll2="jet_2_hpstau", dR=0.4))
        # ?????????????? -> alternative

        selector.set_column("hpstau_1", partial(self.match_nearest, coll1="hps_taus_valid", coll2="jet_1", dR=0.4))
        selector.set_column("hpstau_2", partial(self.match_nearest, coll1="hps_taus_valid", coll2="jet_2", dR=0.4))

        # selector.set_column("jet_1_!hpstau", partial(self.match_nearest, coll1="jet_1", coll2="hps_taus_valid", dR=0.7, not_matched=True))
        # selector.set_column("jet_2_!hpstau", partial(self.match_nearest, coll1="jet_2", coll2="hps_taus_valid", dR=0.7, not_matched=True))

        # selector.set_column("jet_1_gtau_hpstau", partial(self.match_nearest, coll1="jet_1_gtau", coll2="hps_taus_valid", dR=0.4))
        # selector.set_column("jet_2_gtau_hpstau", partial(self.match_nearest, coll1="jet_2_gtau", coll2="hps_taus_valid", dR=0.4))
        
        # selector.set_column("jet_1_!gtau_hpstau", partial(self.match_nearest, coll1="jet_1_!gtau", coll2="hps_taus_valid", dR=0.4))
        # selector.set_column("jet_2_!gtau_hpstau", partial(self.match_nearest, coll1="jet_2_!gtau", coll2="hps_taus_valid", dR=0.4))

        # selector.set_column("jet_1_gtau_!hpstau", partial(self.match_nearest, coll1="jet_1_gtau", coll2="hps_taus_valid", dR=0.4, not_matched=True))
        # selector.set_column("jet_2_gtau_!hpstau", partial(self.match_nearest, coll1="jet_2_gtau", coll2="hps_taus_valid", dR=0.4, not_matched=True))
        
        # selector.set_column("jet_1_!gtau_!hpstau", partial(self.match_nearest, coll1="jet_1_!gtau", coll2="hps_taus_valid", dR=0.4, not_matched=True))
        # selector.set_column("jet_2_!gtau_!hpstau", partial(self.match_nearest, coll1="jet_2_!gtau", coll2="hps_taus_valid", dR=0.4, not_matched=True))

        # Pick up matched SVs
        selector.set_column("SV_jet1", partial(self.match_nearest, coll1="SV", coll2="jet_1", dR=0.4))
        selector.set_column("SV_jet2", partial(self.match_nearest, coll1="SV", coll2="jet_2", dR=0.4))

        # Pick up matched PfCand
        self.pfCand_Sequence(selector, input_jet_n = 1) # jet_1
        self.pfCand_Sequence(selector, input_jet_n = 2) # jet_2
        
        # Fake regions
        # G1G2 - both jets are matched to gen-taus
        # G1F2 - jet1 is matched to gen-tau, jet2 is not matched
        # F1G2 - jet1 is not matched to gen-tau, jet2 is matched
        # F1F2 - both jets are not matched to gen-tau
        selector.set_cat("cat_fake", {"G1G2", "G1F2", "F1G2", "F1F2", "NOMATCH"})
        selector.set_multiple_columns(partial(self.fake_masks, is_mc=is_mc))

        # study regions:
        # PP - 2 OS prompt tau_h
        # PD - 1 prompt tau_h, 1 displaced tau_h
        # DD - 2 displaced tau_h
        # D_INCL - Inclusive for all displacement
        selector.set_cat("cat_disp", {"PP", "PD", "DD", "PDtag", "DDtag", "PP_inv", "PDtag_inv", "DDtag_inv", "D_INCL"})
        selector.set_multiple_columns(partial(self.category_masks))

        # regions of charge:
        # OS - opposite charge
        # SS - same charge
        # S_INCL - Inclusive for all charge
        # selector.set_cat("cat_charge", {"OS", "SS", "S_INCL"})
        # selector.set_multiple_columns(partial(self.charge_masks))

        # divide into the regions related to the gen-matching
        # if is_mc:
            # selector.set_column("jet1_gen", partial(self.jet_gen, jet_name="jet_1"))
            # selector.set_column("jet2_gen", partial(self.jet_gen, jet_name="jet_2"))

        selector.add_cut("b_tagged_0_cut", self.b_tagged_tight_cut)
        
        # selector.set_column("signal_bins", self.signal_bins)