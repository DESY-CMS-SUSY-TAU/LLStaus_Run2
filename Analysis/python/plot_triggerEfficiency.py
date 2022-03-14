import aghast
import awkward
import boost_histogram
import coffea
import coffea.hist
import coffea.processor
import dataclasses
import logging
import numba
import numpy
import os
import sortedcontainers
import uproot3

import ROOT
ROOT.gROOT.SetBatch(True)

#import CMS_lumi, tdrstyle
import utils

from coffea.nanoevents import NanoEventsFactory, NanoAODSchema


@numba.njit
def find_objpair(events_obj1, events_obj2, areSameObjs, builder) :
    
    for objs1, objs2 in zip(events_obj1, events_obj2) :
        
        builder.begin_list()
        
        nObj1 = len(objs1)
        nObj2 = len(objs2)
        
        #print(nObj1, nObj2)
        
        obj1_idx = -1
        obj2_idx = -1
        
        sum_pt_max = 0
        
        for iObj1 in range(nObj1):
            
            obj1 = objs1[iObj1]
            
            if (not (obj1.pt > 20 and abs(obj1.eta) < 2.5)) :
                
                continue
            
            for iObj2 in range(nObj2) :
                
                obj2 = objs2[iObj2]
                
                if (areSameObjs and iObj2 <= iObj1) :
                    
                    continue
                
                if (not (obj2.pt > 20 and abs(obj2.eta) < 2.5)) :
                    
                    continue
                
                chg1 = obj1.charge
                chg2 = obj2.charge
                
                if (chg1 * chg2 > 0) :
                    
                    continue
                
                sum_pt = obj1.pt + obj2.pt
                
                if (
                    (obj1_idx < 0 and obj2_idx < 0) or
                    sum_pt > sum_pt_max
                ) :
                    
                    obj1_idx = iObj1
                    obj2_idx = iObj2
                    
                    sum_pt_max = sum_pt
        
        if (obj1_idx >= 0 and obj2_idx >= 0) :
            
            # Sort only if they're the same objects
            if (areSameObjs and objs2[obj2_idx].pt > objs1[obj1_idx].pt) :
                
                obj1_idx, obj2_idx = obj2_idx, obj1_idx
            
            #print(obj1_idx, obj2_idx)
            
            builder.integer(obj1_idx)
            builder.integer(obj2_idx)
        
        
        builder.end_list()
    
    return builder


@numba.njit
def find_taupair(events_taus, builder):
    
    for taus in events_taus :
        
        builder.begin_list()
        
        nTau = len(taus)
        
        tau1_idx = -1
        tau2_idx = -1
        
        sum_pt_max = 0
        
        for iTau1 in range(nTau):
            
            if (not (taus[iTau1].pt > 20 and abs(taus[iTau1].eta) < 2.5)) :
                
                continue
            
            for iTau2 in range(iTau1+1, nTau) :
                
                if (not (taus[iTau2].pt > 20 and abs(taus[iTau2].eta) < 2.5)) :
                    
                    continue
                
                if (taus[iTau1].charge * taus[iTau2].charge > 0) :
                    
                    continue
                
                sum_pt = taus[iTau1].pt + taus[iTau2].pt
                
                if (
                    (tau1_idx < 0 and tau2_idx < 0) or
                    sum_pt > sum_pt_max
                ) :
                    
                    tau1_idx = iTau1
                    tau2_idx = iTau2
                    
                    sum_pt_max = sum_pt
        
        if (tau1_idx >= 0 and tau2_idx >= 0) :
            
            if (taus[tau2_idx].pt > taus[tau1_idx].pt) :
                
                tau1_idx, tau2_idx = tau2_idx, tau1_idx
            
            #print(tau1_idx, tau2_idx)
            
            builder.integer(tau1_idx)
            builder.integer(tau2_idx)
        
        
        builder.end_list()
    
    return builder


@dataclasses.dataclass
class MyProcessor(coffea.processor.ProcessorABC) :
    
    genObj1opt          : str
    genObj2opt          : str
    trigger_str         : str       = None
    
    
    def __post_init__(self) :
        
        self.dataset_axis = coffea.hist.Cat("dataset", "Dataset")
        
        self._accumulator = coffea.processor.dict_accumulator({
            "MET_pt": coffea.hist.Hist(
                "Counts",
                self.dataset_axis,
                coffea.hist.Bin("MET_pt", "MET_pt", 100, 0, 1000),
            ),
            
            "GenMET_pt": coffea.hist.Hist(
                "Counts",
                self.dataset_axis,
                coffea.hist.Bin("MET_pt", "MET_pt", 100, 0, 1000),
            ),
        })
    
    
    @property
    def accumulator(self) :
        
        return self._accumulator
    
    
    # we will receive a NanoEvents instead of a coffea DataFrame
    def process(self, events) :
        
        output = self.accumulator.identity()
        
        # Will only work for leptons
        events["GenPart", "charge"] = -events.GenPart.pdgId
        
        d_events_obj = {}
        
        for objOpt in [self.genObj1opt, self.genObj2opt] :
            
            if (objOpt == "mu") :
                
                genMuon = events.GenPart[
                    (abs(events.GenPart.pdgId) == 13)
                    & events.GenPart.hasFlags(["isLastCopy"])
                    & events.GenPart.hasFlags(["isTauDecayProduct"])
                ]
                
                d_events_obj[objOpt] = genMuon
            
            elif (objOpt == "tauh") :
                
                d_events_obj[objOpt] = events.GenVisTau
            
            else :
                
                logging.error("%s not a valid option." %(objOpt))
            
        print(self.genObj1opt, self.genObj2opt)
        print(d_events_obj)
        
        objPair_idx = find_objpair(
            events_obj1 = d_events_obj[self.genObj1opt],
            events_obj2 = d_events_obj[self.genObj2opt],
            areSameObjs = (self.genObj1opt == self.genObj2opt),
            builder = awkward.ArrayBuilder(),
        ).snapshot()
        
        # Skip processing as it is an EmptyArray
        if awkward.all(awkward.num(objPair_idx) == 0) :
            
            return output
        
        events = events[awkward.num(objPair_idx, axis = 1) >= 2]
        
        if (self.trigger_str is not None and len(self.trigger_str)) :
            
            #events = events[events.HLT.PFMETNoMu110_PFMHTNoMu110_IDTight | events.HLT.PFMETNoMu120_PFMHTNoMu120_IDTight | ]
            events = events[eval(self.trigger_str)]
        
        output["MET_pt"].fill(
            dataset = events.metadata["dataset"],
            MET_pt = events.MET.pt,
            weight = numpy.ones(len(events)),
        )
        
        output["GenMET_pt"].fill(
            dataset = events.metadata["dataset"],
            MET_pt = events.GenMET.pt,  
            weight = numpy.ones(len(events)),
        )
        
        #print(len(events))
        
        return output
    
    
    def postprocess(self, accumulator):
        
        return accumulator





def main() :
    
    datasets = sortedcontainers.SortedDict({
        "stau100_lsp1_ctau1000mm": ["tmp/nanoaod_stau100_lsp1_ctau1000mm.root"],
        "stau250_lsp1_ctau1000mm": ["tmp/nanoaod_stau250_lsp1_ctau1000mm.root"],
        "stau400_lsp1_ctau1000mm": ["tmp/nanoaod_stau400_lsp1_ctau1000mm.root"],
    })
    
    #trigger_str = ""
    
    trigger_str = (
        "HLT_PFMET120_PFMHT120_IDTight"
        " | HLT_PFMET130_PFMHT130_IDTight "
        " | HLT_PFMET140_PFMHT140_IDTight "
        " | HLT_PFMETNoMu110_PFMHTNoMu110_IDTight "
        " | HLT_PFMETNoMu120_PFMHTNoMu120_IDTight "
        " | HLT_PFMETNoMu130_PFMHTNoMu130_IDTight "
        " | HLT_PFMETNoMu140_PFMHTNoMu140_IDTight "
        
        " | HLT_IsoMu24 "
        #" | HLT_TkMu100 "
        
        #" | HLT_IsoMu27_MET90 "
        
        #" | HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg "
        #" | HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg "
        #" | HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg "
        #" | HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg "
        
        " | HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100 "
        " | HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET110 "
        " | HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET120 "
        " | HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET130 "
        " | HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET140 "
        
        " | HLT_MET105_IsoTrk50 "
        " | HLT_MET120_IsoTrk50 "
        
        #" | HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1 "
        #" | HLT_IsoMu24_eta2p1_MediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_CrossL1 "
        
        #" | HLT_Mu43NoFiltersNoVtx_Photon43_CaloIdL"
        
        #" | HLT_IsoTrackHB "
    )
    
    trigger_str = trigger_str.replace("HLT_", "events.HLT.")
    
    genObj1opt, genObj2opt = "mu", "tauh"
    #genObj1opt, genObj2opt = "tauh", "tauh"
    
    #extraTag = ""
    #extraTag = "excl-HLT-MET-IsoTrk-PFTau"
    extraTag = "incl-HLT-MET-IsoTrk-PFTau"
    
    output_den = coffea.processor.run_uproot_job(
        datasets,
        "Events",
        MyProcessor(
            genObj1opt = genObj1opt,
            genObj2opt = genObj2opt,
        ),
        coffea.processor.iterative_executor,
        {"schema": NanoAODSchema},
    )
    
    output_num = coffea.processor.run_uproot_job(
        datasets,
        "Events",
        MyProcessor(
            trigger_str = trigger_str,
            genObj1opt = genObj1opt,
            genObj2opt = genObj2opt,
        ),
        coffea.processor.iterative_executor,
        {"schema": NanoAODSchema},
    )
    
    
    d_hist = {
        "MET_pt": {
            "xtitle": "p^{miss}_{T} [GeV]"
        },
        
        "GenMET_pt": {
            "xtitle": "p^{miss}_{T, gen} [GeV]"
        },
    }
    
    
    for dataset in datasets:
        
        tag = "{genObj1opt}-{genObj2opt}{extraTag}".format(
            genObj1opt = genObj1opt,
            genObj2opt = genObj2opt,
            extraTag = "_%s" %(extraTag) if len(extraTag) else ""
        )
        
        outdir = "output/plots/%s/%s" %(dataset, tag)
        os.system("mkdir -p %s" %(outdir))
        
        for histName in d_hist.keys() :
            
            print("*"*10, dataset, histName, "*"*10)
            
            l_hist = []
            
            hist_num = output_num[histName].integrate("dataset", dataset)
            hist_den = output_den[histName].integrate("dataset", dataset)
            
            h1_num = aghast.to_root(aghast.from_numpy(hist_num.to_boost().to_numpy()), "%s_%s_num" %(dataset, histName))
            h1_den = aghast.to_root(aghast.from_numpy(hist_den.to_boost().to_numpy()), "%s_%s_den" %(dataset, histName))
            
            count_num = h1_num.Integral()
            count_den = h1_den.Integral()
            
            print("h1_num: integral %f" %(h1_num.Integral()))
            print("h1_den: integral %f" %(h1_den.Integral()))
            
            h1_num.SetLineColor(2)
            h1_num.SetLineWidth(3)
            h1_num.SetMarkerSize(0)
            h1_num.SetTitle("After HLT (frac = %0.2f)" %(count_num/count_den))
            
            h1_den.SetLineColor(4)
            h1_den.SetLineWidth(3)
            h1_den.SetMarkerSize(0)
            h1_den.SetTitle("Before HLT")
            
            l_hist.extend([h1_den, h1_num])
            
            outfile = "%s/%s.pdf" %(outdir, histName)
            
            utils.root_plot1D(
                l_hist = l_hist,
                ratio_num_den_pairs = [(h1_num, h1_den)],
                outfile = outfile,
                xrange = (0, 500),
                yrange = (0, max([_h.GetMaximum() for _h in l_hist])),
                logx = False, logy = False,
                ytitle = "Counts",
                xtitle_ratio = d_hist[histName]["xtitle"], ytitle_ratio = "Efficiency",
                centertitlex = True, centertitley = True,
                centerlabelx = False, centerlabely = False,
                gridx = False, gridy = False,
                ndivisionsx = None,
                stackdrawopt = "nostack",
                legendpos = "UR",
                legendtitle = "#splitline{%s}{(%s)}" %(dataset, tag),
                legendncol = 1,
                #legendtextsize = 0.04,
                legendwidthscale = 1.3,
                legendheightscale = 1.5,
                lumiText = "2018 (13 TeV)"
            )
            #
            #
            #
            #h1_den.Draw("hist same")
            #h1_num.Draw("hist")
            #
            #CMS_lumi.CMS_lumi(canvas, iPeriod = 0, iPosX = 0)
            #CMS_lumi.CMS_lumi(pad = canvas, iPeriod = 0, iPosX = 0, CMSextraText = "Simulation Preliminary", lumiText = "2018 (13 TeV)")
            #
            #
            #canvas.SaveAs(outfile)
        
        
        #outdir_root = "output/root_files/%s" %(dataset)
        #os.system("mkdir -p %s" %(outdir))
        #outfile_root = "%s/output.root" %(outdir)
        #
        #with uproot3.recreate(outfile) as root_file:
        #    
        #    for histName in output_den.keys() :
        #        
        #        hist_num = output_num[histName].integrate("dataset", dataset)
        #        hist_den = output_den[histName].integrate("dataset", dataset)
        #        
        #        h1_num = aghast.to_root(aghast.from_numpy(hist_num.to_boost().to_numpy()), "%s_num" %(histName))
        #        h1_den = aghast.to_root(aghast.from_numpy(hist_den.to_boost().to_numpy()), "%s_den" %(histName))
        #        
        #        #if not isinstance(hist, coffea.hist.Hist) :
        #        #    
        #        #    continue
        #        
        #        #h1_tmp = hist.integrate("dataset", dataset)
        #        #
        #        #print("*"*10, dataset, histName, "*"*10)
        #        #
        #        #root_file[histName] = coffea.hist.export1d(h1_tmp)
        #        
        #    #h1_tmp = hist.integrate("dataset", dataset)
        #    #
        #    #h1_root = aghast.to_root(aghast.from_numpy(h1_tmp.to_boost().to_numpy()), "root_hist")
        #    #h1_root.Draw("hist")
        #    #
        #    #h1_root.SetLineColor(2)
        #    #h1_root.SetLineWidth(3)
        #    #
        #    #print(h1_root.Integral())
        #    #print(h1_root.GetEntries())
        
        #h1_tmp = hist[dataset].copy()
        #
        #h1_root = aghast.to_root(aghast.from_numpy(h1_tmp.to_boost().to_numpy()), "root_hist")
        #h1_root.Draw("hist")
        #
        #h1_root.SetLineColor(2)
        #h1_root.SetLineWidth(3)
        #
        #print(h1_root.Integral())
        #print(h1_root.GetEntries())
    
    
    
    
    
    
    
    #CMS_lumi.CMS_lumi(canvas, iPeriod = 0, iPosX = 0)
    #CMS_lumi.CMS_lumi(pad = canvas, iPeriod = 0, iPosX = 0, CMSextraText = "Simulation Preliminary", lumiText = "2018 (13 TeV)")
    #
    #canvas.SaveAs("plots/test.pdf")
    
    return 0


if __name__ == "__main__" :
    
    main()
