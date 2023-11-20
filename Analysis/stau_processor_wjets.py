from nis import match
from unittest import result
import awkward as ak
from functools import partial
import numpy as np
import mt2
import numba as nb
import coffea
import uproot
import pepper

from coffea.nanoevents import NanoAODSchema
# np.set_printoptions(threshold=np.inf)


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
            
        if "muon_sf" not in config or len(config["muon_sf"]) == 0:
            logger.warning("No muon scale factors specified")

        if "single_mu_trigger_sfs" not in config or\
            len(config["single_mu_trigger_sfs"]) == 0:
            logger.warning("No single muon trigger scale factors specified")

        # It is not recommended to put anything as member variable into a
        # a Processor because the Processor instance is sent as raw bytes
        # between nodes when running on HTCondor.

    def process_selection(self, selector, dsname, is_mc, filler):
        
        
        if self.config["flavour_fake_rate_study"]:
            # Alternative processor which aims
            # calculating jet fake rate basing
            # on the gen level flavour matching
            # used only with WJet dataset at the moment
            selector.systematics["weight"] = \
                ak.full_like(selector.systematics["weight"], 1.0)
            self.process_flav_study(selector, dsname, is_mc, filler)
            return
        
        # Triggers
        pos_triggers, neg_triggers = pepper.misc.get_trigger_paths_for(
            dsname, is_mc, self.config["dataset_trigger_map"],
            self.config["dataset_trigger_order"])
        selector.add_cut("Trigger", partial(
            self.passing_trigger, pos_triggers, neg_triggers))
        
        if is_mc and "pileup_reweighting" in self.config:
            selector.add_cut("Pileup reweighting", partial(
                self.do_pileup_reweighting, dsname))
        
        if is_mc and ( dsname.startswith("WJetsToLNu") or \
                    dsname.startswith("W1JetsToLNu") or \
                    dsname.startswith("W2JetsToLNu") or \
                    dsname.startswith("W3JetsToLNu") or \
                    dsname.startswith("W4JetsToLNu") ):
            selector.systematics["weight"] = \
                ak.full_like(selector.systematics["weight"], 1.0)
            selector.add_cut("W jet reweighting",
                partial(self.do_w_jet_reweighting))

        # MET cut
        selector.add_cut("MET", self.MET_cut)

        selector.add_cut("MET filters", partial(self.met_filters, is_mc))
        
        # one tight tag muon
        selector.set_column("muon_tag", self.muons_tag) 
        selector.add_cut("single_muon_tag", self.single_muon_tag_cut)
        
        selector.add_cut("muon_sfs", partial(self.get_muon_sfs, is_mc=is_mc))
        
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
        selector.set_column("Jet_select", self.jet_selection)
        selector.set_column("Jet_select", self.getloose_jets)

        selector.set_column("PfCands", self.pfcand_valid)
        selector.set_column("Jet_lead_pfcand", partial(self.get_matched_pfCands, match_object="Jet_select", dR=0.4))
        selector.set_column("Jet_select", self.set_jet_dxy)
        selector.add_cut("mt_muon2", self.mt_muon_cut)

        selector.add_cut("two_loose_jets", self.has_two_jets)
        
        # Variables related to the two jets:
        selector.set_column("sum_jj", self.sum_jj)
        selector.set_multiple_columns(self.mt_jets)
        selector.set_column("dphi_jet1_jet2", self.dphi_jet1_jet2)
        selector.add_cut("dphi_min_cut", self.dphi_min_cut)
        selector.set_column("mt2_j1_j2_MET", self.get_mt2)
        
        # Tagger part for calculating scale factors
        # Scale factors should be calculated -
        # before cuts on the number of the jets
        selector.set_multiple_columns(self.set_njets_pass)
        if self.config["predict_yield"]:
            selector.set_multiple_columns(partial(self.predict_yield, weight=selector.systematics["weight"]))

        selector.add_cut("two_loose_jets_final", self.has_two_jets)
        
        # selector.set_column("Jet_select", self.gettight_jets)
        # selector.add_cut("two_tight_jets", self.has_two_jets)        
    
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
    def do_w_jet_reweighting(self, data):
        njet = data["LHE"]["Njets"]
        weights = self.config["W_jet_reweight"][njet]
        return weights

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
            & (muons.pfIsoId >= self.config["muon_pfIsoId"])
            & (abs(muons.dxy) <= self.config["muon_absdxy"])
            & (abs(muons.dz) <= self.config["muon_absdz"])
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
            & (muons.pfIsoId >= self.config["muon_pfIsoId"])
            & (abs(muons.dxy) <= self.config["muon_absdxy"])
            & (abs(muons.dz) <= self.config["muon_absdz"])
            )
        
        return muons[(is_good & (~is_tag))]
        
    
    @zero_handler
    def electron_veto(self, data):
        ele = data["Electron"]
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
            & (jets.jetId >= self.config["jet_jetId"] )
            )]
        matches_h, dRlist = jets.nearest(data["muon_tag"], return_metric=True, threshold=self.config["tag_muon_veto_dR"])
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
        pfCands_lead["ip3d"] = np.sqrt(pfCands_lead.dxy**2 + pfCands_lead.dz**2)
        pfCands_lead["dxy_weight"] = ak.mean(pfCands.dxy, weight=pfCands.pt, axis=-1)
        pfCands_lead["dxysig_weight"] = ak.mean(pfCands.dxy / pfCands.dxyError, weight=pfCands.pt, axis=-1)
        return pfCands_lead
    
    @zero_handler
    def set_jet_dxy(self, data):
        jets = data["Jet_select"]
        # Mask jets with dxy nan (no selected pfcands matching)
        bad_jets = ak.is_none(data["Jet_lead_pfcand"].dxy, axis=-1)
        jets = ak.mask(jets, ~bad_jets) # mask bad jets to keep correct shape
        jets["dz"] = np.abs(data["Jet_lead_pfcand"].dz)
        jets["dxy"] = np.abs(data["Jet_lead_pfcand"].dxy)
        jets["dxy_weight"] = np.abs(data["Jet_lead_pfcand"].dxy_weight)
        jets["dxysig"] = np.abs(data["Jet_lead_pfcand"].dxysig)
        jets["dxysig_weight"] = np.abs(data["Jet_lead_pfcand"].dxysig_weight)
        jets["ip3d"] = data["Jet_lead_pfcand"].ip3d
        jets = jets[~bad_jets] # remove bad jets
        # add dxy cut
        jets = jets[jets.dxy >= self.config["jet_dxy_min"]]
        return jets
    
    @zero_handler
    def n_jets(self, data):
        return ak.num(data["valid_jets"])
   
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
        
        # # from bin 0 to bin 1
        # weights_bin0to1 = []
        # for score in self.config["score_pass"][1:]: # skip first bin because it is just 1
        #     events_0tag = (ak.num(jets[(jets.disTauTag_score1 >= score)]) == 0) # events with 0 tag
        #     jets_notag = (jets.disTauTag_score1 < score) # to have a per jet mask
        #     jets_counted = jets[events_0tag * jets_notag] # to select only jets in events with 0 tag
        #     # fake_sf =  self.config["jet_fake_rate"](jet_dxy=jets_counted.dxy, jet_pt=jets_counted.pt, jet_score=score)
        #     fake_sf = ak.full_like(jets_counted.pt, self.config["jet_fake_rate"](jet_score=[score])[0])
        #     weight_sfs = ak.sum(fake_sf, axis=1)
        #     weights_bin0to1.append(weight_sfs)
        # yield_bin0to1 = ak.from_regular(np.stack(weights_bin0to1, axis=1), axis=-1)

        # # from bin 1 to bin 2
        # weights_bin1to2 = []
        # for score in self.config["score_pass"][1:]: # skip first bin because it is just 1
        #     events_1tag = (ak.num(jets[(jets.disTauTag_score1 >= score)]) == 1) # events with 1 tag
        #     jets_notag = (jets.disTauTag_score1 < score) # to have a per jet mask and not to count the tagged jet
        #     jets_counted = jets[events_1tag * jets_notag]  # to select only jets in events with 1 tag
        #     fake_sf = ak.full_like(jets_counted.pt, self.config["jet_fake_rate"](jet_score=[score])[0])
        #     weight_sfs = ak.sum(fake_sf, axis=1)
        #     weights_bin1to2.append(weight_sfs)
        # yield_bin1to2 = ak.from_regular(np.stack(weights_bin1to2, axis=1), axis=-1)
        
        # # from bin 0 to bin 2
        # weights_bin0to2 = []
        # for score in self.config["score_pass"][1:]: # skip first bin because it is just 1
        #     events_0tag = (ak.num(jets[(jets.disTauTag_score1 >= score)]) == 0) # events with 0 tag
        #     jets_notag = (jets.disTauTag_score1 < score) # to have a per jet mask
        #     jets_counted = jets[events_0tag * jets_notag] # to select only jets in events with 1 tag
        #     fake_sf = ak.full_like(jets_counted.pt, self.config["jet_fake_rate"](jet_score=[score])[0])
        #     combinations = ak.combinations(fake_sf, 2, axis=1) # to have all possible combinations of 2 jets
        #     combinations_unzipped = ak.unzip(combinations)
        #     products = combinations_unzipped[0] * combinations_unzipped[1]
        #     weight_sfs = ak.sum(products, axis=1)
        #     weights_bin0to2.append(weight_sfs)
        # yield_bin0to2 = ak.from_regular(np.stack(weights_bin0to2, axis=1), axis=-1)

        
        weights_bin0to1 = []
        weights_bin0to2 = []
        weights_bin1to2 = []
        
        for score in self.config["score_pass"]:
            
            fake =  self.config["jet_fake_rate"](jet_pt=jets.pt, jet_dxy=jets.dxy)
            # fake =  self.config["jet_fake_rate"](jet_pt=jets.pt)
            
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
        
    @zero_handler
    def mt(self, data, name):
        visible = ak.firsts(data[name])
        MET = data["MET"]
        one_min_cs = 1.0 - np.cos(self.delta_phi(visible.phi, MET.phi))
        prod = 2*visible.pt*MET.pt
        return np.sqrt( prod * one_min_cs)
    
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
        
    @zero_handler
    def has_two_jets(self, data):
        jets = data["Jet_select"]
        return ak.num(jets) == 2

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
    def dphi_min_cut(self, data):
        return abs(data["dphi_jet1_jet2"]) > self.config["dphi_j1j2"]

    # Gen Study -------------------------------------------------------------------------
    
    def process_flav_study(self, selector, dsname, is_mc, filler):
        
        # if not (dsname.startswith("W") or \
        #         dsname.startswith("DY")) :
        #     raise NameError('Error: Gen study of jet flavour available only for WJet.')
        # set event weight to one because of the wronggly esigned
    
        selector.set_multiple_columns(self.gen_lep)    
        # add cuts and selections on the jets
        selector.set_column("Jet_select", self.jet_selection_genveto)
        selector.set_column("Jet_select", self.getloose_jets)
        selector.set_column("PfCands", self.pfcand_valid)
        selector.set_column("Jet_lead_pfcand", partial(self.get_matched_pfCands, match_object="Jet_select", dR=0.4))
        selector.set_column("Jet_select", self.set_jet_dxy)
        selector.add_cut("two_loose_jets", self.has_two_jets)
        # Variables related to the two jets:
        selector.set_column("sum_jj", self.sum_jj)
        selector.set_multiple_columns(self.mt_jets)
        selector.set_column("dphi_jet1_jet2", self.dphi_jet1_jet2)
        selector.set_column("mt2_j1_j2_MET", self.get_mt2)
        selector.add_cut("two_loose_jets_final", self.has_two_jets)
                
    @zero_handler
    def gen_lep(self, data):
        gen_tau = data.GenPart[
            (abs(data.GenPart.pdgId) == 15)
            & data.GenPart.hasFlags(["isHardProcess"])
            & data.GenPart.hasFlags(["isFirstCopy"])
        ]
        gen_mu = data.GenPart[
            (abs(data.GenPart.pdgId) == 13)
            & data.GenPart.hasFlags(["isHardProcess"])
            & data.GenPart.hasFlags(["isFirstCopy"])
        ]
        gen_ele = data.GenPart[
            (abs(data.GenPart.pdgId) == 11)
            & data.GenPart.hasFlags(["isHardProcess"])
            & data.GenPart.hasFlags(["isFirstCopy"])
        ]
        return {
            "gen_tau" : gen_tau,
            "gen_mu"  : gen_mu,
            "gen_ele" : gen_ele
        }
        
    @zero_handler
    def jet_selection_genveto(self, data):
        
        jets = data["Jet"]
        jets = jets[(
            (self.config["jet_eta_min"] < jets.eta)
            & (jets.eta < self.config["jet_eta_max"])
            & (self.config["jet_pt_min"] < jets.pt)
            & (jets.jetId >= self.config["jet_jetId"] )
            )]
        
        # add hadronic tau flavour:
        tau_vis = data.GenVisTau[ ((data.GenVisTau.pt > 30) & (abs(data.GenVisTau.eta) < 2.4) &
                                  (data.GenVisTau.parent.hasFlags(["fromHardProcess"])))
                                ]
        matches_tauhad, _ = jets.nearest(tau_vis, return_metric=True, threshold=self.config["tag_muon_veto_dR"])
        updFlavour = ak.where(~ak.is_none(matches_tauhad, axis=-1), 15, jets.partonFlavour)
        jets["partonFlavour_upd"] = updFlavour
        
        # mask the jets if they match to gen level but not hadron tau
        matches_tau, _ = jets.nearest(data["gen_tau"], return_metric=True, threshold=self.config["tag_muon_veto_dR"])
        matches_mu, _  = jets.nearest(data["gen_mu"], return_metric=True, threshold=self.config["tag_muon_veto_dR"])
        matches_ele, _ = jets.nearest(data["gen_ele"], return_metric=True, threshold=self.config["tag_muon_veto_dR"])
        
        jets = jets[
            ak.is_none(matches_ele, axis=1) &
            ak.is_none(matches_mu, axis=1) &
            ( ak.is_none(matches_tau, axis=1) | (~ak.is_none(matches_tau, axis=1) & ~ak.is_none(matches_tauhad, axis=1)))
        ]
        
        return jets