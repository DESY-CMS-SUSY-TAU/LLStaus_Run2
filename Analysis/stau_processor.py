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

    @zero_handler
    def nan_drop(self, data):
        '''
        Function to drop events with Nan values in the DisTauTag score,
        alternatively one should analyse the whole file and drop it using
        bad_file_paths in the config
        '''
        if np.any(np.logical_not(np.isfinite(data["Jet"].disTauTag_score1))):
            return ak.full_like(data["event"], False)
        return ak.full_like(data["event"], True)

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


    def pfCand_Sequence(self, selector, input_jet_n=None):
        if not input_jet_n:
            raise Exception("Error: input_jet_name should not be None")
        
        jet_n = str(input_jet_n)
        selector.set_column("pfCands_jet"+jet_n, partial(self.get_matched_pfCands, match_object="jet_"+jet_n, dR=0.4))
        # selector.set_column("pfCands_jet"+jet_n+"_gtau", partial(self.get_matched_pfCands, match_object="jet_"+jet_n+"_gtau", dR=0.4))
        # selector.set_column("pfCands_jet"+jet_n+"_!gtau", partial(self.get_matched_pfCands, match_object="jet_"+jet_n+"_!gtau", dR=0.4))

        # selector.set_column("gen_stau",    self.gen_stau)
        # selector.set_column("gen_stau_tau", self.gen_stau_tau)

        # selector.set_column("jet_"+jet_n+"_!gtau_!allgtau", partial(self.match_nearest, coll1="jet_"+jet_n+"_!gtau", coll2="gen_stau_tau", dR=0.7, not_matched=True))
        # selector.set_column("jet_"+jet_n+"_!gtau_!allgtau_!stau", partial(self.match_nearest, coll1="jet_"+jet_n+"_!gtau_!allgtau", coll2="gen_stau", dR=0.7, not_matched=True))

        # # Choose PF Cands that are daughters of isolated jets
        # selector.set_column("pfCands_jet"+jet_n+"_!gtau_!allgtau", partial(self.get_matched_pfCands, match_object="jet_"+jet_n+"_!gtau_!allgtau", dR=0.4))
        # selector.set_column("pfCands_jet"+jet_n+"_!gtau_!allgtau_!stau", partial(self.get_matched_pfCands, match_object="jet_"+jet_n+"_!gtau_!allgtau_!stau", dR=0.4))

        # # Choose PF Cands that match to HPS or not match to HPS
        # selector.set_column("pfCands_jet"+jet_n+"_hpstau", partial(self.get_matched_pfCands, match_object="jet_"+jet_n+"_hpstau", dR=0.4))
        # selector.set_column("pfCands_jet"+jet_n+"_!hpstau", partial(self.get_matched_pfCands, match_object="jet_"+jet_n+"_!hpstau", dR=0.4))

        # selector.set_column("pfCands_jet"+jet_n+"_gtau_hpstau", partial(self.get_matched_pfCands, match_object="jet_"+jet_n+"_gtau_hpstau", dR=0.4))
        # selector.set_column("pfCands_jet"+jet_n+"_gtau_!hpstau", partial(self.get_matched_pfCands, match_object="jet_"+jet_n+"_gtau_!hpstau", dR=0.4))

    @zero_handler
    def gentau_jets(self, data):
        code = data["match_valid_jets"]
        jets = data["jets_valid"][(code == 1) | (code == 5)]
        return jets

    @zero_handler
    def gen_match_jets(self, data, jet_name='', dR=0.4):
        
        #     # -1 not matched to any gen-tau object
        #     # 0 matches to genJet
        #     # 1 hadronic tau
        #     # 2 matched to tau muon
        #     # 3 matched to tau electron
        #     # 4 mutched to >2 objects
        #     # 5 matched to >2 objects including hadronic tau
    
        jet = data[jet_name]

        gen_elec = data.GenPart[ (abs(data.GenPart.pdgId) == 11) &
                                    (data.GenPart.hasFlags(["isDirectTauDecayProduct"])) ]

        gen_muon = data.GenPart[ (abs(data.GenPart.pdgId) == 13) &
                                    (data.GenPart.hasFlags(["isDirectTauDecayProduct"])) ]
        gen_tauvis = data.GenVisTau
        
        matches_genjet, dRlist = jet.nearest(data.GenJet, return_metric=True, threshold=dR)
        hasJet = (~ak.is_none(matches_genjet, axis=1))
        data[jet_name, "match_genJet"] = hasJet

        matches_h, dRlist = jet.nearest(gen_tauvis, return_metric=True, threshold=dR)
        hasHadronic = (~ak.is_none(matches_h, axis=1))
        data[jet_name, "match_genTau"] = hasHadronic

        matches_muon, dRlist = jet.nearest(gen_muon, return_metric=True, threshold=dR)
        hasMuon = (~ak.is_none(matches_muon, axis=1))
        data[jet_name, "match_genMuon"] = hasMuon

        matches_elec, dRlist = jet.nearest(gen_elec, return_metric=True, threshold=dR)
        hasElec = (~ak.is_none(matches_elec, axis=1))
        data[jet_name, "match_genElec"] = hasElec
        
        code = ak.where(data[jet_name].match_genJet, 0, -1)
        code = ak.where(data[jet_name].match_genTau, 1, code)
        code = ak.where(data[jet_name].match_genMuon, 2, code)
        code = ak.where(data[jet_name].match_genElec, 3, code)

        multy_match = ak.values_astype(data[jet_name].match_genJet, "int64") + \
                      ak.values_astype(data[jet_name].match_genTau, "int64") + \
                      ak.values_astype(data[jet_name].match_genMuon, "int64") + \
                      ak.values_astype(data[jet_name].match_genElec, "int64")

        code = ak.where(multy_match > 2, 4, code)
        code = ak.where(((multy_match > 2) & hasHadronic), 5, code)

        return code

    def control_region(self, data):
        if len(data) == 0:
            return { "M1" : ak.Array([]),
                     "M2" : ak.Array([]),
                     "M3" : ak.Array([]),
                     "ALL" : ak.Array([]) }
        mask = {}
        
        # mask["SR"] = ak.firsts(data["sum_2jets"].mass) > self.config["control_region"]
        # mask["CR"] = ak.firsts(data["sum_2jets"].mass) < self.config["control_region"]
        mask["M1"] = data["mt_sum"] < self.config["region1"]
        mask["M2"] = (( data["mt_sum"] > self.config["region1"] ) & ( data["mt_sum"] < self.config["region2"] ))
        mask["M3"] = data["mt_sum"] > self.config["region2"]
        mask["ALL"] = ak.full_like(mask["M1"], True)
        
        return mask
    
    def fake_masks(self, data, is_mc):
        if(len(data) == 0):
            return { "G1G2" : ak.Array([]),
                     "G1F2" : ak.Array([]),
                     "F1G2" : ak.Array([]),
                     "F1F2" : ak.Array([]),
                     "NOMATCH" : ak.Array([]) }
        mask = {}
        Inclusive = np.full(len(data["hps_taus_valid"]), True)
        if is_mc:
            tau_vis = data.GenVisTau[ ((data.GenVisTau.pt > 30) & (abs(data.GenVisTau.eta) < 2.4) &
                                      (data.GenVisTau.parent.hasFlags(["fromHardProcess"])))
                                    ] 
            matches_h, dRlist = data["jet_1"].nearest(data.GenVisTau, return_metric=True, threshold=0.4)
            has_hps1 = ak.firsts(~ak.is_none(matches_h, axis=1))
            matches_h, dRlist = data["jet_2"].nearest(data.GenVisTau, return_metric=True, threshold=0.4)
            has_hps2 = ak.firsts(~ak.is_none(matches_h, axis=1))
            mask["G1G2"] = ( has_hps1 & has_hps2 )
            mask["G1F2"] = ( has_hps1 & (~has_hps2) )
            mask["F1G2"] = ( (~has_hps1) & has_hps2 )
            mask["F1F2"] = ( (~has_hps1) & (~has_hps2) )
        else:
            mask["G1G2"] = Inclusive
            mask["G1F2"] = Inclusive
            mask["F1G2"] = Inclusive
            mask["F1F2"] = Inclusive
        mask["NOMATCH"] = Inclusive
        return mask

    def category_masks(self, data):
        if len(data) == 0:
            return { "PP" : ak.Array([]),
                     "PD" : ak.Array([]),
                     "PDtag" : ak.Array([]),
                     "DD" : ak.Array([]),
                     "DDtag" : ak.Array([]),
                     "PP_inv" : ak.Array([]),
                     "PDtag_inv" : ak.Array([]),
                     "DDtag_inv" : ak.Array([]),
                     "D_INCL" : ak.Array([]) }

        masks = {}


        # To have jet 1 associated with hps taus
        sel_1 = ( (ak.num(data['hpstau_1'])>=1) & (ak.firsts(data['hpstau_1']).idDeepTau2017v2p1VSjet >= 6) )
        sel_1 = ak.fill_none(sel_1, False)
        sel_2 = ( np.abs(data['pfCands_jet1'].dxy) < self.config["pf_dxy"] )
        sel_2 = ak.fill_none(sel_2, False)
        sel_3 = ( np.abs(data['pfCands_jet1'].dz) < self.config["pf_dz"] )
        sel_3 = ak.fill_none(sel_3, False)
        isJet1_Prompt = ( sel_1 & sel_2 & sel_3 )
        
        
        sel_1_inv = ( (ak.num(data['hpstau_1'])>=1) & (ak.firsts(data['hpstau_1']).idDeepTau2017v2p1VSjet < 6) & (ak.firsts(data['hpstau_1']).idDeepTau2017v2p1VSjet >= 3))
        sel_1_inv = ak.fill_none(sel_1_inv, False)
        isJet1_Prompt_inv = ( sel_1_inv & sel_2 & sel_3 )
        

        # To have jet 2 associated with hps taus
        sel_1 = ( (ak.num(data['hpstau_2'])>=1) & (ak.firsts(data['hpstau_2']).idDeepTau2017v2p1VSjet >= 6) )
        sel_1 = ak.fill_none(sel_1, False)
        sel_2 = ( np.abs(data['pfCands_jet2'].dxy) < self.config["pf_dxy"] )
        sel_2 = ak.fill_none(sel_2, False)
        sel_3 = ( np.abs(data['pfCands_jet2'].dz) < self.config["pf_dz"] )
        sel_3 = ak.fill_none(sel_3, False)
        isJet2_Prompt = ( sel_1 & sel_2 & sel_3 )

        sel_1_inv = ( (ak.num(data['hpstau_2'])>=1) & (ak.firsts(data['hpstau_2']).idDeepTau2017v2p1VSjet < 6) & (ak.firsts(data['hpstau_2']).idDeepTau2017v2p1VSjet >= 3))
        sel_1_inv = ak.fill_none(sel_1_inv, False)
        isJet2_Prompt_inv = ( sel_1_inv & sel_2 & sel_3 )

        # To have jet 1 displaced
        isJet1_Displaced = ( np.abs(data['pfCands_jet1'].dxy) > self.config["pf_dxy"] )
        isJet1_Displaced = ak.fill_none(isJet1_Displaced, False)
        isJet1_Tagged = ( ak.firsts(data['jet_1'].disTauTag_score1) > 0.90 )
        isJet1_Tagged = ak.fill_none(isJet1_Tagged, False)
        isJet1_Dtag = (isJet1_Displaced & isJet1_Tagged)

        isJet1_Tagged_inv = (( ak.firsts(data['jet_1'].disTauTag_score1) < 0.90 ) & ( ak.firsts(data['jet_1'].disTauTag_score1) > 0.50 ) )
        isJet1_Tagged_inv = ak.fill_none(isJet1_Tagged_inv, False)
        isJet1_Dtag_inv = (isJet1_Displaced & isJet1_Tagged_inv)

        # To have jet 2 displaced
        isJet2_Displaced = ( np.abs(data['pfCands_jet2'].dxy) > self.config["pf_dxy"] )
        isJet2_Displaced = ak.fill_none(isJet2_Displaced, False)
        isJet2_Tagged = ( ak.firsts(data['jet_2'].disTauTag_score1) > 0.90 )
        isJet2_Tagged = ak.fill_none(isJet2_Tagged, False)
        isJet2_Dtag = (isJet2_Displaced & isJet2_Tagged)

        isJet2_Tagged_inv = (( ak.firsts(data['jet_2'].disTauTag_score1) < 0.90 ) & ( ak.firsts(data['jet_2'].disTauTag_score1) > 0.50 ) )
        isJet2_Tagged_inv = ak.fill_none(isJet2_Tagged_inv, False)
        isJet2_Dtag_inv = (isJet2_Displaced & isJet2_Tagged_inv)

        # PP category
        PP = (isJet1_Prompt & isJet2_Prompt)
        masks["PP"] = PP
        
        # PP category with inverted deepTau ID
        PP_inv = ((isJet1_Prompt_inv & isJet2_Prompt) | (isJet1_Prompt_inv & isJet2_Prompt_inv) | (isJet1_Prompt & isJet2_Prompt_inv))
        masks["PP_inv"] = PP_inv

        # PD - 1 prompt tau_h, 1 displaced tau_h
        PD = ( (isJet1_Prompt & isJet2_Displaced) | (isJet2_Prompt & isJet1_Displaced) )
        masks["PD"] = PD

        # DD - 2 displaced tau_h
        DD = (isJet1_Displaced & isJet2_Displaced)
        masks["DD"] = DD

        PDtag = ( (isJet1_Prompt & isJet2_Dtag) | (isJet2_Prompt & isJet1_Dtag) )
        masks["PDtag"] = PDtag
        
        PDtag_inv = (((isJet1_Prompt_inv & isJet2_Dtag) | (isJet1_Prompt_inv & isJet2_Dtag_inv) | (isJet1_Prompt & isJet2_Dtag_inv)) |
                     ((isJet2_Prompt_inv & isJet1_Dtag) | (isJet2_Prompt_inv & isJet1_Dtag_inv) | (isJet2_Prompt & isJet1_Dtag_inv)))
        masks["PDtag_inv"] = PDtag_inv

        DDtag = (isJet2_Dtag & isJet1_Dtag)
        masks["DDtag"] = DDtag
        
        DDtag_inv = ((isJet1_Dtag & isJet2_Dtag_inv) | (isJet1_Dtag_inv & isJet2_Dtag_inv) | (isJet1_Dtag_inv & isJet2_Dtag))  
        masks["DDtag_inv"] = DDtag_inv

        # COMB - all regions together
        masks["D_INCL"] = ak.full_like(DD, True)

        return masks

    def charge_masks(self, data):
        if len(data) == 0:
            return { "OS" : ak.Array([]),
                     "SS" : ak.Array([]),
                     "S_INCL" : ak.Array([]) }
        # SS = (ak.firsts(data['pfCands_jet1'],axis=-1).charge
        #     * ak.firsts(data['pfCands_jet2'],axis=-1).charge) > 0
        # OS = (ak.firsts(data['pfCands_jet1'],axis=-1).charge
        #     * ak.firsts(data['pfCands_jet2'],axis=-1).charge) < 0
        SS = (data['pfCands_jet1'].charge
            * data['pfCands_jet2'].charge) > 0
        OS = (data['pfCands_jet1'].charge
            * data['pfCands_jet2'].charge) < 0
        SS = ak.fill_none(SS, False)
        OS = ak.fill_none(OS, False)
        ALL_CHARGES = ak.full_like(SS, True)
        return {"OS": OS, "SS": SS, "S_INCL": ALL_CHARGES}

    @zero_handler
    def signal_bins(self, data):
        

        bins = np.full((len(data["D_INCL"])), np.nan)
        
        mt2 = data["mt2_j1_j2_MET"][:,0]
        mt_sum = data["mt_sum"]
        
        # DD region bins:
        DDMTTWO1 = (mt2 < 100)
        DDMTTWO1SMT1 = data["DDtag"] & DDMTTWO1 & ( mt_sum > 0 )   & ( mt_sum < 200)
        DDMTTWO1SMT2 = data["DDtag"] & DDMTTWO1 & ( mt_sum > 200 ) & ( mt_sum < 400) 
        DDMTTWO1SMT3 = data["DDtag"] & DDMTTWO1 & ( mt_sum > 400 ) & ( mt_sum < 800) 
        DDMTTWO1SMT4 = data["DDtag"] & DDMTTWO1 & ( mt_sum > 800 )
        DDMTTWO2SMT1 = data["DDtag"] & (mt2 > 100)
        
        # PD region bins:
        PDMTTWO1 =   (mt2 < 100) 
        PDMTTWO2 = ( (mt2 > 100) & (mt2 < 200) )
        PDMTTWO1SMT1 = data["PDtag"] & PDMTTWO1 & ( mt_sum > 0 )   & ( mt_sum < 200) 
        PDMTTWO1SMT2 = data["PDtag"] & PDMTTWO1 & ( mt_sum > 200 ) & ( mt_sum < 400)  
        PDMTTWO1SMT3 = data["PDtag"] & PDMTTWO1 & ( mt_sum > 400 ) & ( mt_sum < 600)
        PDMTTWO1SMT4 = data["PDtag"] & PDMTTWO1 & ( mt_sum > 600 ) & ( mt_sum < 800)
        PDMTTWO1SMT5 = data["PDtag"] & PDMTTWO1 & ( mt_sum > 800 )
        PDMTTWO2SMT1 = data["PDtag"] & PDMTTWO2 & ( mt_sum < 800) 
        PDMTTWO2SMT2 = data["PDtag"] & PDMTTWO2 & ( mt_sum > 800 )
        PDMTTWO3SMT1 = data["PDtag"] & (mt2 > 200)
        
        # PP region binning
        PPMTTWO1SMT1 = data["PP"] & ( (mt_sum > 200 ) & (mt_sum < 250) ) & ( (mt2 > 25) & (mt2 < 50) ) 
        PPMTTWO2SMT1 = data["PP"] & ( (mt_sum > 200 ) & (mt_sum < 250) ) & (mt2 > 50)
        PPMTTWO1SMT2 = data["PP"] & ( (mt_sum > 250 ) & (mt_sum < 300) ) & ( (mt2 > 25) & (mt2 < 50) )
        PPMTTWO2SMT2 = data["PP"] & ( (mt_sum > 250 ) & (mt_sum < 300) ) & ( (mt2 > 50) & (mt2 < 75) )
        PPMTTWO3SMT2 = data["PP"] & ( (mt_sum > 250 ) & (mt_sum < 300) ) & (mt2 > 75)
        PPMTTWO1SMT3 = data["PP"] & ( (mt_sum > 300 ) & (mt_sum < 350) ) & ( (mt2 > 25) & (mt2 < 50) )
        PPMTTWO2SMT3 = data["PP"] & ( (mt_sum > 300 ) & (mt_sum < 350) ) & ( (mt2 > 50) & (mt2 < 75) )
        PPMTTWO3SMT3 = data["PP"] & ( (mt_sum > 300 ) & (mt_sum < 350) ) & (mt2 > 75)
        PPMTTWO1SMT4 = data["PP"] & (mt_sum > 350 ) & ( (mt2 > 25) & (mt2 < 75) )
        PPMTTWO2SMT4 = data["PP"] & (mt_sum > 350 ) & ( (mt2 > 75) & (mt2 < 100) )
        PPMTTWO3SMT4 = data["PP"] & (mt_sum > 350 ) & (mt2 > 100)
        
        bins[PPMTTWO1SMT1] = 1
        bins[PPMTTWO2SMT1] = 2
        bins[PPMTTWO1SMT2] = 3
        bins[PPMTTWO2SMT2] = 4
        bins[PPMTTWO3SMT2] = 5
        bins[PPMTTWO1SMT3] = 6
        bins[PPMTTWO2SMT3] = 7
        bins[PPMTTWO3SMT3] = 8
        bins[PPMTTWO1SMT4] = 9
        bins[PPMTTWO2SMT4] = 10
        bins[PPMTTWO3SMT4] = 11
        
        bins[PDMTTWO1SMT1] = 12
        bins[PDMTTWO1SMT2] = 13
        bins[PDMTTWO1SMT3] = 14
        bins[PDMTTWO1SMT4] = 15
        bins[PDMTTWO1SMT5] = 16
        bins[PDMTTWO2SMT1] = 17
        bins[PDMTTWO2SMT2] = 18
        bins[PDMTTWO3SMT1] = 19
        
        bins[DDMTTWO1SMT1] = 20
        bins[DDMTTWO1SMT2] = 21
        bins[DDMTTWO1SMT3] = 22
        bins[DDMTTWO1SMT4] = 23
        bins[DDMTTWO2SMT1] = 24

        return bins 
        
    @zero_handler
    def jets_valid(self, data):
        jets = data["Jet"]
        jets["initial_idx"] = ak.local_index(jets, axis=1)

        # study kinem region
        jets = jets[(
            (self.config["jet_eta_min"] < jets.eta)
            & (jets.eta < self.config["jet_eta_max"])
            & (self.config["jet_pt_min"] < jets.pt)
            )]

        # muon veto
        jets_vetoed = jets.nearest(data["muons_valid"], return_metric=True, 
                                   threshold=self.config["muon_veto_dR"])[0]
        idx_vetoed_jets = ak.is_none(jets_vetoed, axis=-1)
        jets = jets[idx_vetoed_jets]

        # electron veto
        jets_vetoed = jets.nearest(data["electrons_valid"], return_metric=True,
                                   threshold=self.config["ele_veto_dR"])[0]
        idx_vetoed_jets = ak.is_none(jets_vetoed, axis=-1)
        jets = jets[idx_vetoed_jets]

        # sort by pt
        sort_idx = ak.argsort(jets.pt, axis=-1, ascending=False)
        jets = jets[sort_idx]
        
        # to index the valid jets
        
        # Add DeepTau ID to every jet that matches to HPS tau
        # if there is no matching jet, set the ID to -1
        # taus = data["hps_taus_valid"]
        # nearest_tau = jets.nearest(taus, return_metric=True, threshold=0.4)[0]
        # jets["idDeepTau2017v2p1VSjet"] = ak.where(ak.is_none(nearest_tau, axis=-1), -1, nearest_tau.idDeepTau2017v2p1VSjet)
        # jets["idDeepTau2017v2p1VSe"] = ak.where(ak.is_none(nearest_tau, axis=-1), -1, nearest_tau.idDeepTau2017v2p1VSe)
        # jets["idDeepTau2017v2p1VSmu"] = ak.where(ak.is_none(nearest_tau, axis=-1), -1, nearest_tau.idDeepTau2017v2p1VSmu)
        # jets["idDeepTau2018v2p5VSjet"] = ak.where(ak.is_none(nearest_tau, axis=-1), -1, nearest_tau.idDeepTau2018v2p5VSjet)
        # jets["idDeepTau2018v2p5VSe"] = ak.where(ak.is_none(nearest_tau, axis=-1), -1, nearest_tau.idDeepTau2018v2p5VSe)
        # jets["idDeepTau2018v2p5VSmu"] = ak.where(ak.is_none(nearest_tau, axis=-1), -1, nearest_tau.idDeepTau2018v2p5VSmu)
        # print(jets["initial_idx"])
        # print(data["hps_taus_valid"].jetIdx)
        # ids = wrap_padder(jets, taus)
        # exit()
        
        jets["valid_idx"] = ak.local_index(jets, axis=1)
        
        return jets

    @zero_handler
    def muons_valid(self, data):
        muons = data["Muon"]
        is_good = (
              (muons.pt > self.config["muon_pt"])
            & (muons.eta < self.config["muon_eta_max"])
            & (muons.eta > self.config["muon_eta_min"])
            # & (muons[ self.config["muonID"] == 1)
            & (muons.pfRelIso04_all < self.config["muon_pfRelIso04_all"])
            )
        return muons[is_good]

    @zero_handler
    def electrons_valid(self, data):
        ele = data["Electron"]
        ele_low_eta_iso = eval(self.config["ele_low_eta_iso"])
        ele_high_eta_iso = eval(self.config["ele_high_eta_iso"])
        isolation_cut = ( ele_low_eta_iso | ele_high_eta_iso )
        is_good = (
            isolation_cut
            & (ele.pt > self.config["muon_pt"])
            & (ele.eta < self.config["ele_eta_max"])
            & (ele.eta > self.config["ele_eta_min"])
            & (ele[self.config["eleVeto"]] == 1)
            # & (ele[self.config["eleID"]]==1)
            )
        return ele[is_good]

    @zero_handler
    def hps_taus_valid(self, data):
        taus = data["Tau"]
        is_good = (
            (taus.pt > self.config["tau_pt"])
            & (taus.eta < self.config["tau_eta_max"])
            & (taus.eta > self.config["tau_eta_min"])
            # & (taus.idDeepTau2017v2p1VSjet >= self.config["tau_idDeepTau_vsjet"])
            # & (taus.idDeepTau2017v2p1VSmu >= self.config["tau_idDeepTau_vsmu"])
            # & (taus.idDeepTau2017v2p1VSe >= self.config["tau_idDeepTau_vsele"])
            # & (taus.idDeepTau2018v2p5VSe >= self.config["tau_idDeepTau_vsele"])
            # & (taus.idDeepTau2018v2p5VSjet >= self.config["tau_idDeepTau_vsjet"])
            # & (taus.idDeepTau2018v2p5VSmu >= self.config["tau_idDeepTau_vsmu"])
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
        # print("score outside [0..1]: ", ak.any((jets.disTauTag_score1>1.0 | jets.disTauTag_score1<0.0 )))
        return jets

    @zero_handler
    def add_j1_j2(self, data):
        return data["jet_1"].add(data["jet_2"])

    @zero_handler
    def b_tagged_jet(self, data):
        # To leading jets are excluded! 
        jet_not_signal = data["jets_valid"][:,2:]
        # Jet_btagDeepFlavB satisfies the Medium (>0.2783) WP:
        # https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL18
        b_tagged_idx = (jet_not_signal.btagDeepFlavB > 0.2783)
        return jet_not_signal[b_tagged_idx]

    @zero_handler
    def b_tagged_cut(self, data):
        return ak.num(data["jet_b"]) <= 1
    
    @zero_handler
    def b_tagged_tight_cut(self, data):
        return ak.num(data["jet_b"]) == 0

    @zero_handler    
    def MET_cut(self, data):
        return data["MET"].pt > self.config["MET_pt"] #GeV

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
        # invisible_parts = data.GenPart[ (abs(data.GenPart.pdgId) == 1000022) & 
        #                                 (data.GenPart.hasFlags(["isLastCopy"])) ]
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
    def pfcand_valid(self, data):
        pfCands = data["PFCandidate"]
        is_good = (
            (pfCands.pt > self.config["pfcand_pt"])
            & (pfCands.eta < self.config["pfcand_eta_max"])
            & (pfCands.eta > self.config["pfcand_eta_min"])
            & (pfCands[self.config["track"]])
        )
        return pfCands[is_good]

    @zero_handler
    def get_matched_pfCands(self, data, match_object, dR=0.4):
        pfCands = self.match_nearest(data, coll1="pfcand_valid", coll2=match_object, dR=dR)
        # print(pfCands)
        pfCands = ak.firsts(pfCands, axis=-1) # take only leading pfCand (optional)
        # print(pfCands)
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

    @zero_handler
    def jet_momentum_sum(self, data):
        jets = data["jets_valid"]
        HT = ak.sum(jets.pt, axis=-1)
        return HT
    
    @zero_handler
    def jet_momentum_sum_miss(self, data):
        jets = data["jets_valid"]
        px = ak.sum(jets.px, axis=-1)
        py = ak.sum(jets.py, axis=-1)
        return np.sqrt(px*px + py*py)
    
    def delta_phi(self, phi1_ak, phi2_ak):
        phi1 = np.array(phi1_ak)
        phi2 = np.array(phi2_ak)
        assert phi1.shape == phi2.shape
        d = phi1 - phi2
        indx_pos = d>np.pi
        d[indx_pos] -= np.pi*2
        indx_neg = d<=-np.pi
        d[indx_neg] += np.pi*2
        return d

    @zero_handler
    def mt(self, data, name):
        visible = ak.firsts(data[name])
        MET = data["MET"]
        one_min_cs = 1.0 - np.cos(self.delta_phi(visible.phi, MET.phi))
        prod = 2*visible.pt*MET.pt
        return np.sqrt( prod * one_min_cs) 
    
    @zero_handler
    def HT_miss_pt_miss(self, data):
        return data["HT_miss"]/data["MET"].pt
    
    @zero_handler
    def mt_sum(self, data):
        return data["mt_jet1"]+data["mt_jet2"]
    
    @zero_handler
    def dphi_jet1_jet2(self, data):
        return self.delta_phi(ak.firsts(data["jet_1"].phi),
                              ak.firsts(data["jet_2"].phi))
