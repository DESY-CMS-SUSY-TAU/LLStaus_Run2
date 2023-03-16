#!/usr/bin/env python3

import os

import ROOT
ROOT.gROOT.SetBatch(True)

import utils.utils


def main() :
    
    infilename = "output_stau_kinematics_study.root"
    infile = ROOT.TFile.Open(infilename, "READ")
    
    outdir = "plots/stau_kinematics_study"
    os.system(f"mkdir -p {outdir}")
    
    l_stautype = [
        "stau_LH",
        "stau_RH",
        "stau_MM",
    ]
    
    l_histname = [
        "GenVisTauh_pt",
        "GenVisTaul_pt",
        "GenVisTaul_e_by_stau_mass",
    ]
    
    d_xrange = {}
    d_xrange["GenVisTauh_pt"] = (0, 500)
    d_xrange["GenVisTaul_pt"] = (0, 500)
    d_xrange["GenVisTaul_e_by_stau_mass"] = (0, 0.5)
    
    d_yrange = {}
    d_yrange["GenVisTauh_pt"] = (0, 0.1)
    d_yrange["GenVisTaul_pt"] = (0, 0.2)
    d_yrange["GenVisTaul_e_by_stau_mass"] = (0, 0.05)
    
    d_hist = {}
    
    for stautype in l_stautype :
        
        d_hist[stautype] = {}
        
        for histname in l_histname :
            
            hist = infile.Get(f"{stautype}/{histname}").Clone()
            hist.SetDirectory(0)
            
            hist.Scale(1.0 / hist.GetEntries())
            
            d_hist[stautype][histname] = hist
    
    infile.Close()
    
    for histname in l_histname :
        
        l_hist = []
        outfilename = f"{outdir}/{histname}.pdf"
        
        xrange = d_xrange[histname]
        yrange = d_yrange[histname]
        
        for istautype, stautype in enumerate(l_stautype) :
            
            hist = d_hist[stautype][histname].Clone()
            hist.SetTitle(stautype)
            hist.SetLineColor(istautype+1)
            hist.SetLineWidth(3)
            hist.SetMarkerSize(0)
            
            l_hist.append(hist)
        
        utils.utils.root_plot1D(
            l_hist = l_hist,
            #ratio_num_den_pairs = [(h1_num, h1_den)],
            outfile = outfilename,
            xrange = xrange,
            yrange = yrange,
            logx = False, logy = False,
            xtitle = histname,
            ytitle = "a.u.",
            #xtitle_ratio = d_hist[histName]["xtitle"], ytitle_ratio = "Efficiency",
            centertitlex = True, centertitley = True,
            centerlabelx = False, centerlabely = False,
            gridx = True, gridy = False,
            ndivisionsx = None,
            stackdrawopt = "nostack",
            legendpos = "UR",
            legendtitle = "[stau(250), lsp(1)]",
            legendncol = 1,
            #legendtextsize = 0.04,
            legendwidthscale = 0.8,
            #legendheightscale = 1.5,
            lumiText = "2018 (13 TeV)"
        )


if __name__ == "__main__" :
    
    main()