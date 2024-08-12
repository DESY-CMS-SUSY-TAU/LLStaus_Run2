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

ext = "png"

output_dir = "./output_fake_compare_run2/"

period_names = [
    "2016preVFP",
    "2016postVFP",
    "2017",
    "2018",
]

dy_fake_rates = [
    "/nfs/dust/cms/user/mykytaua/softLLSTAU/LLStaus_Run2/Analysis/output_iteration_4/2016/output_zmumu_2016pre/zmumu_v1_fake_yield/fake_rate_ext/",
    "/nfs/dust/cms/user/mykytaua/softLLSTAU/LLStaus_Run2/Analysis/output_iteration_4/2016/output_zmumu_2016post/zmumu_v1_fake_yield/fake_rate_ext/",
    "/nfs/dust/cms/user/mykytaua/softLLSTAU/LLStaus_Run2/Analysis/output_iteration_4/2017/output_zmumu/zmumu_v1/fake_rate_ext/",
    "/nfs/dust/cms/user/mykytaua/softLLSTAU/LLStaus_Run2/Analysis/output_iteration_4/2018/output_zmumu/zmumu_v1/fake_rate_ext/",
]

wjet_fake_rates = [
    "/nfs/dust/cms/user/mykytaua/softLLSTAU/LLStaus_Run2/Analysis/output_iteration_4/2016/output_wjet_2016preVFP/wjet_v2_fake/fake_rate_ext/",
    "/nfs/dust/cms/user/mykytaua/softLLSTAU/LLStaus_Run2/Analysis/output_iteration_4/2016/output_wjet_2016postVFP/wjet_v2_fake/fake_rate_ext/",
    "/nfs/dust/cms/user/mykytaua/softLLSTAU/LLStaus_Run2/Analysis/output_iteration_4/2017/output_wjet/wjet_v1_fake/fake_rate_ext/",
    "/nfs/dust/cms/user/mykytaua/softLLSTAU/LLStaus_Run2/Analysis/output_iteration_4/2018/output_wjet/wjet_fake_v1/fake_rate_ext/",

]

histograms = { # hist : [tree_name, label, [xmin, xmax], [ymin, ymax]]
    "fake_rate_jet_dxy.root" : ["fake_rate_jet_dxy", "d_{xy} (cm)", [0.5, 35.0], [0, 1.2]],
    "fake_rate_jet_pt.root" : ["fake_rate_jet_pt", "p_{T} (GeV)", [30.0, 1000], [0, 0.5]],
}



# plot 1D histograms
if __name__ == "__main__":
    
    try: os.stat(output_dir)
    except: os.mkdir(output_dir)
    
    for period, (dy_path, wjet_path) in enumerate( zip(dy_fake_rates, wjet_fake_rates) ):
        for hist_name in list(histograms.keys()):
        
            tree_name = histograms[hist_name][0]
            label = histograms[hist_name][1]    
            xmin = histograms[hist_name][2][0]
            xmax = histograms[hist_name][2][1]
            ymin = histograms[hist_name][3][0]
            ymax = histograms[hist_name][3][1]
        
            file = ROOT.TFile.Open(dy_path + hist_name, 'read')
            hist_data = file.Get(tree_name)
            hist_data.SetDirectory(0)
            hist_data.SetMarkerStyle(22)
            hist_data.SetMarkerSize(2)
            hist_data.SetMarkerColor(4)
            hist_data.SetLineColor(4)
            hist_data.SetFillColor(0)
            hist_data.SetLineWidth(2)
            hist_data.SetTitle("Z(#mu#mu)+jets fake rate")
            hist_data.Print("all")

            file_data = ROOT.TFile.Open(wjet_path + hist_name, 'read')
            hist_prediction = file_data.Get(tree_name)
            hist_prediction.SetDirectory(0)
            hist_prediction.SetMarkerStyle(21)
            hist_prediction.SetMarkerSize(2)
            hist_prediction.SetMarkerColor(2)
            hist_prediction.SetLineColor(2)
            hist_prediction.SetFillColor(0)
            hist_prediction.SetLineWidth(2)
            hist_prediction.SetTitle("W(#mu#nu)+jets fake rate")
            hist_prediction.Print("all")

            root_plot1D(
                l_hist = [hist_prediction, hist_data],
                l_hist_overlay = [],
                outfile = f"{output_dir}/{tree_name}_{period_names[period]}.{ext}",
                xrange = [xmin, xmax],
                yrange = (ymin, ymax),
                # yrange = (0.0, 0.2),
                logx = True, logy = False,
                logx_ratio = False, logy_ratio = False,
                include_overflow = False,
                # xtitle = "Jet pt (GeV)",
                xtitle = label,
                ytitle = "FR.",
                # xtitle_ratio = "Jet pt (GeV)",
                xtitle_ratio = label,
                ytitle_ratio = "ratio",
                centertitlex = True, centertitley = True,
                centerlabelx = False, centerlabely = False,
                gridx = True, gridy = True,
                ndivisionsx = None,
                stackdrawopt = "nostack",
                # normilize = True,
                normilize_overlay = False,
                legendpos = "UL",
                legendtitle = f"",
                legendncol = 1,
                legendtextsize = 0.05,
                legendwidthscale = 2.0,
                legendheightscale = 2.5,
                lumiText = f"{period_names[period]} (13 TeV)",
                CMSextraText = "Private work",
                signal_to_background_ratio = False,
                draw_errors = False
            )
