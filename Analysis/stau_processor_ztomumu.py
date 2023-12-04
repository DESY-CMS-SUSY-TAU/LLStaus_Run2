# This file illustrates how to implement a processor, realizing the selection
# steps and outputting histograms and a cutflow with efficiencies.
# Here we create a very simplified version of the ttbar-to-dilep processor.
# One can run this processor using
# 'python3 -m pepper.runproc --debug example_processor.py example_config.json'
# Above command probably will need a little bit of time before all cuts are
# applied once. This is because a chunk of events are processed simultaneously.
# You change adjust the number of events in a chunk and thereby the memory
from nis import match
import pepper
import awkward as ak
from functools import partial
import numpy as np
import mt2
import numba as nb
import coffea
import uproot
import os

from coffea.nanoevents import NanoAODSchema

class Processor(pepper.ProcessorBasicPhysics):
    # We use the ConfigTTbarLL instead of its base Config, to use some of its
    # predefined extras
    config_class = pepper.ConfigSTau
    
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

        if "pileup_reweighting" not in config:
            logger.warning("No pileup reweigthing specified")

        if "DY_lo_sfs" not in config:
            logger.warning("No DY k-factor specified")

        if "DY_jet_reweight" not in config or len(config["DY_jet_reweight"]) == 0:
            logger.warning("No DY jet-binned sample stitching is specified")

        if "muon_sf" not in config or len(config["muon_sf"]) == 0:
            logger.warning("No muon scale factors specified")

        if "single_mu_trigger_sfs" not in config or\
            len(config["single_mu_trigger_sfs"]) == 0:
            logger.warning("No single muon trigger scale factors specified")
        
        if "jet_fake_rate" not in config and len(config["jet_fake_rate"]) == 0:
            self.predict_jet_fakes = False
            logger.warning("No jet fake rate specified")
        else:
            self.predict_jet_fakes = config["predict_yield"]

    def process_selection(self, selector, dsname, is_mc, filler):
        
        # Triggers
        pos_triggers, neg_triggers = pepper.misc.get_trigger_paths_for(
            dsname,
            is_mc,
            self.config["dataset_trigger_map"],
            self.config["dataset_trigger_order"])
        selector.add_cut("Trigger", partial(
            self.passing_trigger, pos_triggers, neg_triggers))
        
        if is_mc and ( dsname.startswith("DYJetsToLL_M-50") or \
                       dsname.startswith("DY1JetsToLL_M-50") or \
                       dsname.startswith("DY2JetsToLL_M-50") or \
                       dsname.startswith("DY3JetsToLL_M-50") or \
                       dsname.startswith("DY4JetsToLL_M-50") ):
            selector.add_cut("DY jet reweighting",
                partial(self.do_dy_jet_reweighting))

        if is_mc and "pileup_reweighting" in self.config:
            selector.add_cut("Pileup reweighting", partial(
                self.do_pileup_reweighting, dsname))

        # MET cut
        selector.add_cut("MET", self.MET_cut_max)
        selector.add_cut("MET filters", partial(self.met_filters, is_mc))
        
        # # veto all loose muons
        # selector.set_column("muon_veto", self.muons_veto)
        # selector.add_cut("loose_muon_veto", self.loose_muon_veto_cut)
        
        # # veto all loose electrons
        # selector.set_column("electron_veto", self.electron_veto)
        # selector.add_cut("loose_electron_veto", self.loose_electron_veto_cut)
        
        # 2 muons
        selector.set_column("Muon", self.select_muons) 
        selector.add_cut("two_muons", self.two_muons_cut)
        selector.set_column("sum_mumu", self.sum_mumu)
        if is_mc:
            selector.set_column("sum_ll_gen", self.sum_ll_gen)

        selector.add_cut("dy_gen_sfs", partial(self.get_dy_gen_sfs, is_mc=is_mc, dsname=dsname))
        selector.add_cut("muon_sfs", partial(self.get_muon_sfs, is_mc=is_mc))

        selector.add_cut("mass_window", self.mass_window)
        selector.add_cut("charge", self.charge)
        
        selector.set_column("Jet_select", self.jet_selection)
        selector.add_cut("has_more_two_jets", self.has_more_two_jets)
        selector.set_column("Jet_obj1", partial(self.leading_jet, order=0))
        selector.set_column("Jet_obj2", partial(self.leading_jet, order=1))
        selector.set_column("Jet_select", self.getloose_jets)
 
        selector.set_column("PfCands", self.pfcand_valid)
        selector.set_column("Jet_lead_pfcand", partial(self.get_matched_pfCands, match_object="Jet_select", dR=0.4))
        selector.set_column("Jet_select", self.set_jet_dxy)
        selector.add_cut("charge2", self.charge)
        
        selector.add_cut("two_loose_jets", self.has_two_jets)

        # Variables related to the two jets:
        selector.set_column("sum_jj", self.sum_jj)
        selector.set_column("mt2_j1_j2_MET", self.get_mt2)
        selector.set_multiple_columns(self.missing_energy)
        selector.set_multiple_columns(self.mt_jets)
        selector.set_column("dphi_jet1_jet2", self.dphi_jet1_jet2)
        selector.add_cut("dphi_min_cut", self.dphi_min_cut)

        # selector.set_column("logger", self.skim_jets)
        
        # Tagger part for calculating scale factors
        # Scale factors should be calculated -
        # before cuts on the number of the jets
        selector.set_multiple_columns(self.set_njets_pass)
        if self.config["predict_yield"]:
            selector.set_multiple_columns(partial(self.predict_yield, weight=selector.systematics["weight"]))
        
        selector.set_column("Jet_obj1", partial(self.leading_jet, order=0))
        selector.set_column("Jet_obj2", partial(self.leading_jet, order=1))
        
        selector.add_cut("two_loose_jets_final", self.has_two_jets)
        
        selector.set_column("Jet_select", self.gettight_jets)
        selector.add_cut("two_tight_jets", self.has_two_jets)
    
    # @zero_handler
    # def muons_veto(self, data):
    #     muons = data["Muon"]
        
    #     is_good = (
    #           (muons.pt > self.config["muon_veto_pt_min"])
    #         & (muons.eta < self.config["muon_veto_eta_max"])
    #         & (muons.eta > self.config["muon_veto_eta_min"])
    #         & (muons[self.config["muon_veto_ID"]] == 1)
    #         )
        
    #     is_tag =(
    #           (muons.pt > self.config["muon_pt_min"])
    #         & (muons.eta < self.config["muon_eta_max"])
    #         & (muons.eta > self.config["muon_eta_min"])
    #         & (muons[self.config["muon_ID"]] == 1)
    #         & (muons.pfIsoId >= self.config["muon_pfIsoId"])
    #         & (abs(muons.dxy) <= self.config["muon_absdxy"])
    #         & (abs(muons.dz) <= self.config["muon_absdz"])
    #         )
        
    #     return muons[(is_good & (~is_tag))]
        
    # @zero_handler
    # def loose_muon_veto_cut(self, data):
    #     return ak.num(data["muon_veto"])==0
    
    # @zero_handler
    # def loose_electron_veto_cut(self, data):
    #     return ak.num(data["electron_veto"])==0
    
    # @zero_handler
    # def electron_veto(self, data):
    #     ele = data["Electron"]
    #     # ele_low_eta_iso =  ((np.abs(ele.eta) < 1.479) & (ele.pfRelIso03_all < (0.198+0.506/ele.pt)))
    #     # ele_high_eta_iso = ((np.abs(ele.eta) > 1.479) & (np.abs(ele.eta) < 2.5) & (ele.pfRelIso03_all < (0.203+0.963/ele.pt)))
    #     # isolation_cut = ( ele_low_eta_iso | ele_high_eta_iso )
    #     is_good = (
    #         # isolation_cut
    #         & (ele.pt > self.config["elec_veto_pt"])
    #         & (ele.eta < self.config["elec_veto_eta_min"])
    #         & (ele.eta > self.config["elec_veto_eta_max"])
    #         & (ele[self.config["elec_Veto"]] == 1)
    #         )
    #     return ele[is_good]
    
    @zero_handler
    def sum_jj(self, data):
        return data['Jet_select'][:,0].add(data['Jet_select'][:,1])
    
    @zero_handler
    def missing_energy(self, data):
        jets = data["Jet_select"]
        HT_valid = ak.sum(jets.pt, axis=-1)
    
        px = ak.sum(jets.px, axis=-1)
        py = ak.sum(jets.py, axis=-1)
        HT_miss_valid = np.sqrt(px*px + py*py)

        jets = data["Jet"]
        HT = ak.sum(jets.pt, axis=-1)

        px = ak.sum(jets.px, axis=-1)
        py = ak.sum(jets.py, axis=-1)
        HT_miss= np.sqrt(px*px + py*py)
        
        return {
            "HT_valid" : HT_valid, "HT_miss_valid" : HT_miss_valid,
            "HT" : HT, "HT_miss" : HT_miss
        }
    
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
    def mt_jets(self, data):
        jet1 = data["Jet_select"][:,0]
        jet2 = data["Jet_select"][:,1]
        
        MET = data["MET"]
        one_min_cs = 1.0 - np.cos(self.delta_phi(jet1.phi, MET.phi))
        prod = 2*jet1.pt*MET.pt
        mt_j1 = np.sqrt( prod * one_min_cs)
        
        one_min_cs = 1.0 - np.cos(self.delta_phi(jet2.phi, MET.phi))
        prod = 2*jet2.pt*MET.pt
        mt_j2 = np.sqrt( prod * one_min_cs) 
    
        return {
            "mt_jet1" : mt_j1,
            "mt_jet2" : mt_j2,
            "mt_sum" : mt_j1 + mt_j2
        }
        
    @zero_handler
    def dphi_jet1_jet2(self, data):
        return self.delta_phi(data["Jet_select"][:,0].phi,
                              data["Jet_select"][:,1].phi)
    
    
    @zero_handler    
    def get_mt2(self, data):
        jet1 = data["Jet_select"][:,0]
        jet2 = data["Jet_select"][:,1]
        met = data["MET"]
        return mt2.mt2(
            jet1.mass, jet1.px, jet1.py,
            jet2.mass, jet2.px, jet2.py,
            met.px, met.py,
            0, 0
        )
    
    # @zero_handler
    # def skim_jets(self, data):
    #     jets = data["Jet_select"]
    #     file_name = os.path.basename(data.metadata["filename"])
    #     dataset_name = data.metadata["dataset"]
    #     prefix = f"evn_{data.metadata['entrystart']}_{data.metadata['entrystop']}_"
    #     path = f"{self.config["skim_path"]}/{dataset_name}"
    #     if not os.path.exists(path):
    #         os.makedirs(path)
    #     ak.to_parquet(jets, f"{path}/"+prefix+file_name+".parquet")
    #     return np.ones(len(data))
    
    @zero_handler
    def has_two_jets(self, data):
        jets = data["Jet_select"]
        return ak.num(jets) == 2
    
    @zero_handler
    def has_more_two_jets(self, data):
        jets = data["Jet_select"]
        return ak.num(jets) >= 2

    @zero_handler
    def getloose_jets(self, data):
        jets = data["Jet_select"]
        jets = jets[(jets.disTauTag_score1 >= self.config["loose_thr"])]
        return jets

    @zero_handler
    def gettight_jets(self, data):
        jets = data["Jet_select"]
        jets = jets[(jets.disTauTag_score1 >= self.config["tight_thr"])]
        return jets

    @zero_handler
    def do_dy_jet_reweighting(self, data):
        njet = data["LHE"]["Njets"]
        weights = self.config["DY_jet_reweight"][njet]
        return weights

    @zero_handler
    def MET_cut_max(self, data):
        return data["MET"].pt < self.config["MET_cut_max"]
    
    @zero_handler
    def select_muons(self, data):
        muons = data["Muon"]
        is_good = (
              (muons.pt > self.config["muon_pt_min"])
            & (muons.eta < self.config["muon_eta_max"])
            & (muons.eta > self.config["muon_eta_min"])
            & (muons[self.config["muon_ID"]] == 1)
            & (muons.pfIsoId >= self.config["muon_pfIsoId"])
            & (abs(muons.dxy) <= self.config["muon_absdxy"])
            & (abs(muons.dz) <= self.config["muon_absdz"])
            )
        muons = muons[is_good]
        
        trig_muons = data["TrigObj"]
        trig_muons = trig_muons[
              (abs(trig_muons.id) == 13)
            & (trig_muons.filterBits >= 8)
            ]
        
        matches, dRlist = muons.nearest(trig_muons, return_metric=True, threshold=0.4)
        idx_matches_muon = ~ak.is_none(matches, axis=1)
        results = muons[idx_matches_muon]
        
        return results
    
    @zero_handler
    def two_muons_cut(self, data):
        return ak.num(data["Muon"]) == 2
    
    @zero_handler
    def sum_mumu(self, data):
        return data['Muon'][:,0].add(data['Muon'][:,1])
    
    @zero_handler
    def sum_ll_gen(self, data):
        part = data["GenPart"]
        part = part[ part.hasFlags("isLastCopy")
                & (part.hasFlags("fromHardProcess")
                & ((abs(part["pdgId"]) == 11) 
                | (abs(part["pdgId"]) == 13)
                | (abs(part["pdgId"]) == 12)
                | (abs(part["pdgId"]) == 14)
                | (abs(part["pdgId"]) == 15)
                | (abs(part["pdgId"]) == 16)))
                | (part.hasFlags("isDirectHardProcessTauDecayProduct"))
        ]
        sum_p4 = part.sum(axis=1) # sum over all particles in event
        return sum_p4

    @zero_handler
    def get_dy_gen_sfs(self, data, is_mc, dsname):
        weight = np.ones(len(data))
        if is_mc and dsname.startswith("DY"):
            z_boson = data["sum_ll_gen"]
            dy_gen_sfs = self.config["DY_lo_sfs"](mass=z_boson.mass, pt=z_boson.pt)
            weight *= ak.to_numpy(dy_gen_sfs)
            return weight
        else:
            return weight

    @zero_handler
    def mass_window(self, data):
        is_good = (
              (data["sum_mumu"].mass > self.config["mass_ll_lower"])
            & (data["sum_mumu"].mass < self.config["mass_ll_upper"])
            )
        return is_good
    
    @zero_handler
    def charge(self, data):
        return data["sum_mumu"].charge == 0
    
    @zero_handler
    def get_muon_sfs(self, data, is_mc):
        weight = np.ones(len(data))
        if is_mc:
            id_iso_sfs, systematics = self.muon_id_iso_sfs(data)
            weight *= ak.to_numpy(id_iso_sfs)
            mu_trigger_sfs = self.apply_mu_trigger_sfs(data) # systematics are not implemented yet
            weight *= ak.to_numpy(mu_trigger_sfs)
            return weight, systematics
        else:
            return weight

    def apply_mu_trigger_sfs(self, data):
        weight = np.ones(len(data))
        # Single muon trigger efficiency applied for leading muon
        muon = ak.firsts(data["Muon"])
        sfs = self.config["single_mu_trigger_sfs"][0](pt=muon.pt, abseta=abs(muon.eta))
        return weight * sfs

    def muon_id_iso_sfs(self, data):
        """Compute identification and isolation scale factors for
           leptons (electrons and muons)."""
        muons = data["Muon"]
        weight = np.ones(len(data))
        systematics = {}
        # Muon identification and isolation efficiency
        for i, sffunc in enumerate(self.config["muon_sf"]):
            params = {}
            for dimlabel in sffunc.dimlabels:
                if dimlabel == "abseta":
                    params["abseta"] = abs(muons.eta)
                else:
                    params[dimlabel] = getattr(muons, dimlabel)
            central = ak.prod(sffunc(**params), axis=1)
            key = f"muonsf{i}"
            if self.config["compute_systematics"]:
                if ("split_muon_uncertainty" not in self.config
                        or not self.config["split_muon_uncertainty"]):
                    unctypes = ("",)
                else:
                    unctypes = ("stat ", "syst ")
                for unctype in unctypes:
                    up = ak.prod(sffunc(
                        **params, variation=f"{unctype}up"), axis=1)
                    down = ak.prod(sffunc(
                        **params, variation=f"{unctype}down"), axis=1)
                    systematics[key + unctype.replace(" ", "")] = (
                        up / central, down / central)
            weight = weight * central
        return weight, systematics

    @zero_handler
    def jet_selection(self, data):
        jets = data["Jet"]
        jets = jets[(
            (self.config["jet_eta_min"] < jets.eta)
            & (jets.eta < self.config["jet_eta_max"])
            & (self.config["jet_pt_min"] < jets.pt)
            & (jets.jetId >= self.config["jet_jetId"] )
            )]
        matches_h, dRlist = jets.nearest(data["Muon"], return_metric=True, threshold=0.4)
        isoJets = jets[ak.is_none(matches_h, axis=-1)]
        return isoJets
    
    @zero_handler
    def pfcand_valid(self, data):
        pfCands = data["PFCandidate"]
        is_good = (
            (pfCands.pt > self.config["pfcand_pt"])
            & (pfCands.eta < self.config["pfcand_eta_max"])
            & (pfCands.eta > self.config["pfcand_eta_min"])
            & (pfCands[self.config["track"]])
        )
        pfCands_selected = pfCands[is_good]
        sort_idx = ak.argsort(pfCands_selected.pt, axis=-1, ascending=False)
        return pfCands_selected[sort_idx]
    
    @zero_handler
    def match_jet_to_pfcand(self, data, jet_name = None, pf_name = None, dR = 0.4):
        '''
        This function match all particle-flow candidates to every jet in the [jet_name] collection
        (function is alternative to match_nearest, but works for multydim case and expected to be much slower)
        '''
        jets = data[jet_name]
        pfcands = data[pf_name]
        # here return_combinations=True is needed to return _pfcands_unzipped
        # which broadcasted in the way to duplicate pdcands list per every jet in event
        (dr, (_, _pfcands_unzipped)) = jets.metric_table(pfcands, metric=coffea.nanoevents.methods.vector.LorentzVector.delta_r, return_combinations=True)
        pfcands_matched = _pfcands_unzipped[(dr < dR)]
        return pfcands_matched
    
    @zero_handler
    def get_matched_pfCands(self, data, match_object, dR=0.4):
        pfCands = self.match_jet_to_pfcand(data, jet_name=match_object, pf_name="PfCands", dR=dR)
        pfCands_lead = ak.firsts(pfCands, axis=-1)
        pfCands_lead["dxysig"] = pfCands_lead.dxy / pfCands_lead.dxyError
        pfCands_lead["Lrel"] = np.sqrt(pfCands_lead.dxy**2 + pfCands_lead.dz**2)
        pfCands_lead["dxy_weight"] = ak.mean(pfCands.dxy, weight=pfCands.pt, axis=-1)
        pfCands_lead["dxysig_weight"] = ak.mean(pfCands.dxy / pfCands.dxyError, weight=pfCands.pt, axis=-1)
        return pfCands_lead
    
    @zero_handler
    def set_jet_dxy(self, data):
        jets = data["Jet_select"]
        # Mask jets with dxy nan (no selected pfcands matching)
        bad_jets = ak.is_none(data["Jet_lead_pfcand"].dxy, axis=-1)
        jets = ak.mask(jets, ~bad_jets) # mask bad jets to keep coorect shape
        jets["dz"] = np.abs(data["Jet_lead_pfcand"].dz)
        jets["dxy"] = np.abs(data["Jet_lead_pfcand"].dxy)
        jets["dxy_weight"] = np.abs(data["Jet_lead_pfcand"].dxy_weight)
        jets["dxysig"] = np.abs(data["Jet_lead_pfcand"].dxysig)
        jets["dxysig_weight"] = np.abs(data["Jet_lead_pfcand"].dxysig_weight)
        jets = jets[~bad_jets] # remove bad jets
        jets = jets[jets.dxy >= self.config["jet_dxy_min"]]
        return jets
    
    @zero_handler
    def leading_jet(self, data, order=0):
        jets = data["Jet_select"]
        idx_leading = \
            ak.argsort(jets.pt, ascending=False)[:,order:order+1]
        jets = jets[idx_leading]
        return jets
    
    @zero_handler
    def jets_available(self, data, n_available=1):
        jets = data["Jet_select"][~ak.is_none(data["Jet_select"].pt, axis=-1)]
        return ak.num(jets) >= n_available
    
    @zero_handler
    def dphi_min_cut(self, data):
        return abs(data["dphi_jet1_jet2"]) > self.config["dphi_j1j2"]

    @zero_handler
    def set_njets_pass(self, data):
        jets_score = data["Jet_select"].disTauTag_score1
        n_pass = []
        for score in self.config["score_pass"]:
            jets_pass = jets_score[(jets_score>=score)]
            passed = ak.num(jets_pass, axis=1)
            n_pass.append( passed )
        n_pass = ak.from_regular(
            np.stack(n_pass, axis=1), axis=-1)
        tight_wp = self.config["score_pass"].index(self.config["tight_thr"])
        tight_pass = n_pass[:,tight_wp]
        tight_bin0 = (n_pass[:,tight_wp] == 0)
        tight_bin1 = (n_pass[:,tight_wp] == 1)
        tight_bin2 = (n_pass[:,tight_wp] == 2)
        return {
            "n_pass" : n_pass,
            "n_pass_score_bin" : ak.local_index(n_pass, axis=1),
            "tight_pass" : tight_pass,
            "tight_bin0" : tight_bin0,
            "tight_bin1" : tight_bin1,
            "tight_bin2" : tight_bin2
        }
        
    @zero_handler
    def predict_yield(self, data, weight=None):
        jets = data["Jet_select"]
        
        weights_bin0to1 = []
        weights_bin0to2 = []
        weights_bin1to2 = []
        
        for score in self.config["score_pass"]:
            
            fake =  self.config["jet_fake_rate"](jet_pt=jets.pt, jet_dxy=jets.dxy)
            
            # from bin 0 to bin 1 and 2
            events_0tag = (ak.num(jets[(jets.disTauTag_score1 >= score)]) == 0)
            masked_jets = ak.broadcast_arrays(jets.disTauTag_score1, True)[1]
            masked_jets = events_0tag * masked_jets
            fake0 = ak.mask(fake, masked_jets)
            f_1, f_2 = fake0[:,0], fake0[:,1]
            from0to1 = ( f_1*(1-f_2) + f_2*(1-f_1) ) / ((1-f_2)*(1-f_1))
            from0to2 = ( f_1*f_2 ) / ((1-f_2)*(1-f_1))
            from0to1 = ak.fill_none(from0to1, 0.0)
            from0to2 = ak.fill_none(from0to2, 0.0)
            weights_bin0to1.append(from0to1)
            weights_bin0to2.append(from0to2)
            
            # from bin 1 to bin 2
            events_1tag = (ak.num(jets[(jets.disTauTag_score1 >= score)]) == 1)
            masked_jets = ak.broadcast_arrays(jets.disTauTag_score1, True)[1]
            masked_jets = events_1tag * masked_jets
            fake1 = ak.mask(fake, masked_jets)
            f_1, f_2 = fake1[:,0], fake1[:,1]
            from1to2 = ( f_1*f_2 ) / (f_1*(1-f_2) + f_2*(1-f_1))
            from1to2 = ak.fill_none(from1to2, 0.0)
            weights_bin1to2.append(from1to2)
            
        yield_bin0to1 = ak.from_regular(np.stack(weights_bin0to1, axis=1), axis=-1)
        yield_bin0to2 = ak.from_regular(np.stack(weights_bin0to2, axis=1), axis=-1)
        yield_bin1to2 = ak.from_regular(np.stack(weights_bin1to2, axis=1), axis=-1)
            
        # now we need to each predicted yield assign cooresponding score bin
        score_bin = ak.local_index(yield_bin0to1, axis=1) + 1 # +1 because we skip first bin
        
        # One of the WP (called tight) is used for the final analysis
        # the weight for this will be saved saparetly to have one weight per event
        tight_wp = self.config["score_pass"].index(self.config["tight_thr"])
        
        return {"yield_bin0to1" : weight*yield_bin0to1,
                "yield_bin1to2" : weight*yield_bin1to2,
                "yield_bin0to2" : weight*yield_bin0to2,
                "tight_yield_bin0to1" : weight*yield_bin0to1[:,tight_wp],
                "tight_yield_bin1to2" : weight*yield_bin1to2[:,tight_wp],
                "tight_yield_bin0to2" : weight*yield_bin0to2[:,tight_wp],
                "score_bin"     : score_bin}