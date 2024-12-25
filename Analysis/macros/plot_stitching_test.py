import os
import json
from argparse import ArgumentParser
import re
import ROOT
import itertools

from pepper import Config

import sys
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))
sys.path.append(parent_dir)
from utils.plotter import plot1D, plot2D
from utils.utils import ColorIterator, root_plot1D, root_plot2D

######################################################################################

# ext = "pdf"
# output_dir = "./output_stritching/"
# name_file = "DY_stitching_2017"
# limi = 41.4

# Inclusive_path = "/afs/desy.de/user/m/mykytaua/nfscms/softLLSTAU/LLStaus_Run2/Analysis/output_iteration_4/2017/output_zmumu/zmumu_v3_stitching_v2_done/hists/Cut_000_BeforeCuts_LHE_Vpt.root"
# Inclusive_hist = "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8"
# Inclusive_scale = 6233550 * limi / 103344974

# Combined_path = "/afs/desy.de/user/m/mykytaua/nfscms/softLLSTAU/LLStaus_Run2/Analysis/output_iteration_4/2017/output_zmumu/zmumu_v3_stitching_v2_done/hists/Cut_001_DY jet reweighting_LHE_Vpt.root"
# Combined_hists = [
#     "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8",
#     "DY1JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8",
#     "DY2JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8",
#     "DY3JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8",
#     "DY4JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8"
#     ]
# Combined_scales = limi
# legend = "Z(#mu#mu)+jets"
# year = "2017"

######################################################################################

# ext = "png"
# output_dir = "./output_stritching/"
# name_file = "Wjet_stitching_2017"
# limi = 41.4

# Inclusive_path = "/afs/desy.de/user/m/mykytaua/nfscms/softLLSTAU/LLStaus_Run2/Analysis/output_iteration_5/2017/output_wjet/wjet_v1_stitch/hists/Cut_001_reset_weight_LHE_Vpt.root"
# Inclusive_hist = "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8"
# Inclusive_scale = 61526700 * limi / 78981243

# Combined_path = "/afs/desy.de/user/m/mykytaua/nfscms/softLLSTAU/LLStaus_Run2/Analysis/output_iteration_5/2017/output_wjet/wjet_v1_stitch/hists/Cut_002_W jet reweighting_LHE_Vpt.root"
# Combined_hists = [
#     "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
#     "W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
#     "W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
#     "W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
#     "W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8"
#     ]
# Combined_scales = limi
# legend = "W(#mu#nu)+jets"
# year = "2017"

######################################################################################

# ext = "png"
# name_file = "Wjet_stitching_2016pre"
# output_dir = "./output_stritching/"
# limi = 19.5

# Inclusive_path = "/afs/desy.de/user/m/mykytaua/nfscms/softLLSTAU/LLStaus_Run2/Analysis/output_iteration_5/2016pre/output_wjet_2016preVFP/wjet_v1_stitch/hists/Cut_002_reset_weight_LHE_Vpt.root"
# Inclusive_hist = "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8"
# Inclusive_scale = 61526700 * limi / 74676454.0

# Combined_path = "/afs/desy.de/user/m/mykytaua/nfscms/softLLSTAU/LLStaus_Run2/Analysis/output_iteration_5/2016pre/output_wjet_2016preVFP/wjet_v1_stitch/hists/Cut_003_W jet reweighting_LHE_Vpt.root"
# Combined_hists = [
#     "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
#     "W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
#     "W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
#     "W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
#     "W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8"
#     ]
# Combined_scales = limi
# legend = "W(#mu#nu)+jets"
# year = "2016pre"

######################################################################################

# ext = "png"
# name_file = "DY_stitching_2016pre"
# output_dir = "./output_stritching/"
# limi = 19.5

# Inclusive_path = "/afs/desy.de/user/m/mykytaua/nfscms/softLLSTAU/LLStaus_Run2/Analysis/output_iteration_5/2016pre/output_wjet_2016preVFP/wjet_v1_stitch/hists/Cut_000_BeforeCuts_LHE_Vpt.root"
# Inclusive_hist = "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8"
# Inclusive_scale = 6233550 * limi / 95170542

# Combined_path = "/afs/desy.de/user/m/mykytaua/nfscms/softLLSTAU/LLStaus_Run2/Analysis/output_iteration_5/2016pre/output_wjet_2016preVFP/wjet_v1_stitch/hists/Cut_001_DY jet reweighting_LHE_Vpt.root"
# Combined_hists = [
#     "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8",
#     "DY1JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8",
#     "DY2JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8",
#     "DY3JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8",
#     "DY4JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8"
#     ]
# Combined_scales = limi
# legend = "Z(#mu#mu)+jets"
# year = "2016pre"

######################################################################################

# ext = "png"
# name_file = "Wjet_stitching_2016post"
# output_dir = "./output_stritching/"
# limi = 19.5
# 
# Inclusive_path = "/nfs/dust/cms/user/mykytaua/softLLSTAU/LLStaus_Run2/Analysis/output_iteration_5/2016post/output_wjet/wjet_v1_stitch/hists/Cut_002_reset_weight_LHE_Vpt.root"
# Inclusive_hist = "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8"
# Inclusive_scale = 61526700 * limi / 82572775.0
# 
# Combined_path = "/nfs/dust/cms/user/mykytaua/softLLSTAU/LLStaus_Run2/Analysis/output_iteration_5/2016post/output_wjet/wjet_v1_stitch/hists/Cut_003_W jet reweighting_LHE_Vpt.root"
# Combined_hists = [
#     "WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
#     "W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
#     "W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
#     "W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8",
#     "W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8"
#     ]
# Combined_scales = limi
# legend = "W(#mu#nu)+jets"
# year = "2016post"

######################################################################################

ext = "png"
name_file = "DY_stitching_2016post"
output_dir = "./output_stritching/"
limi = 19.5

Inclusive_path = "/nfs/dust/cms/user/mykytaua/softLLSTAU/LLStaus_Run2/Analysis/output_iteration_5/2016post/output_wjet/wjet_v1_stitch/hists/Cut_000_BeforeCuts_LHE_Vpt.root"
Inclusive_hist = "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8"
Inclusive_scale = 6233550 * limi / 82448537

Combined_path = "/nfs/dust/cms/user/mykytaua/softLLSTAU/LLStaus_Run2/Analysis/output_iteration_5/2016post/output_wjet/wjet_v1_stitch/hists/Cut_001_DY jet reweighting_LHE_Vpt.root"
Combined_hists = [
    "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8",
    "DY1JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8",
    "DY2JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8",
    "DY3JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8",
    "DY4JetsToLL_M-50_MatchEWPDG20_TuneCP5_13TeV-madgraphMLM-pythia8"
    ]
Combined_scales = limi
legend = "Z(#mu#mu)+jets"
year = "2016post"

if __name__ == "__main__":
    
    try: os.stat(output_dir)
    except: os.mkdir(output_dir)

    file = ROOT.TFile.Open(Inclusive_path, 'read')
    hist_data = file.Get(Inclusive_hist+"/hist")
    hist_data.Scale(Inclusive_scale)
    hist_data.SetDirectory(0)
    hist_data.SetMarkerStyle(22)
    hist_data.SetMarkerSize(2)
    hist_data.SetMarkerColor(4)
    hist_data.SetLineColor(4)
    hist_data.SetFillColor(0)
    hist_data.SetLineWidth(2)
    hist_data.SetTitle(legend+" inclusive")
    # hist_data.Print("all")

    file_data = ROOT.TFile.Open(Combined_path, 'read')
    sum_hist = None
    for hist_name in Combined_hists:
        hist_prediction = file_data.Get(hist_name+"/hist")
        hist_prediction.Scale(Combined_scales)
        hist_prediction.SetDirectory(0)
        print(hist_prediction.GetEntries())
        if sum_hist == None: sum_hist = hist_prediction.Clone()
        else: sum_hist.Add(hist_prediction)
    sum_hist.SetMarkerStyle(21)
    sum_hist.SetMarkerSize(2)
    sum_hist.SetMarkerColor(2)
    sum_hist.SetLineColor(2)
    sum_hist.SetFillColor(0)
    sum_hist.SetLineWidth(2)
    sum_hist.SetTitle(legend+" combined")
    # sum_hist.Print("all")

    root_plot1D(
        l_hist = [sum_hist],
        l_hist_overlay = [hist_data],
        outfile = f"{output_dir}/{name_file}.{ext}",
        xrange = (10, 200),
        yrange = (0, sum_hist.GetBinContent(7)*1.4),
        # yrange = (0.0, 0.2),
        logx = False, logy = False,
        logx_ratio = False, logy_ratio = False,
        include_overflow = False,
        # xtitle = "Jet pt (GeV)",
        xtitle = "LHE Vp_{t} [GeV]",
        ytitle = "Events",
        # xtitle_ratio = "Jet pt (GeV)",
        xtitle_ratio = "LHE V p_{T} [GeV]",
        ytitle_ratio = "incl./stitched",
        centertitlex = True, centertitley = True,
        centerlabelx = False, centerlabely = False,
        gridx = True, gridy = True,
        ndivisionsx = None,
        stackdrawopt = "nostack",
        # normilize = True,
        normilize_overlay = False,
        legendpos = "UR",
        legendtitle = f"",
        legendncol = 1,
        legendtextsize = 0.05,
        legendwidthscale = 2.0,
        legendheightscale = 2.3,
        lumiText = f"{year} (13 TeV)",
        CMSextraText = "Internal",
        signal_to_background_ratio = True,
        ratio_mode = "B",
        yrange_ratio = (0.95, 1.05),
        ndivisionsy_ratio = (5, 5, 0),
        draw_errors = False
    )