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

ext = "png"

output_dir = "./output_fake_2D_run2/"

period_names = [
    "2016preVFP",
    "2016postVFP",
    "2017",
    "2018",
]

wjet_fake_rates = [
    "/nfs/dust/cms/user/mykytaua/softLLSTAU/LLStaus_Run2/Analysis/output_iteration_4/2016/output_wjet_2016preVFP/wjet_v2_fake/fake_rate_ext/",
    "/nfs/dust/cms/user/mykytaua/softLLSTAU/LLStaus_Run2/Analysis/output_iteration_4/2016/output_wjet_2016postVFP/wjet_v2_fake/fake_rate_ext/",
    "/nfs/dust/cms/user/mykytaua/softLLSTAU/LLStaus_Run2/Analysis/output_iteration_4/2017/output_wjet/wjet_v1_fake/fake_rate_ext/",
    "/nfs/dust/cms/user/mykytaua/softLLSTAU/LLStaus_Run2/Analysis/output_iteration_4/2018/output_wjet/wjet_fake_v1/fake_rate_ext/",

]

histograms = { # hist : [tree_name, label_x, label_y, [xmin, xmax], [ymin, ymax]]
    "fake_rate_jet_dxy_pt.root" : ["fake_rate_jet_dxy_pt", "p_{T} [GeV]", "d_{xy} [cm]", [0.5, 30.0], [20, 1000]],
}

# plot 1D histograms
if __name__ == "__main__":
    
    try: os.stat(output_dir)
    except: os.mkdir(output_dir)
    
    for period,  wjet_path in enumerate( wjet_fake_rates ):
        for hist_name in list(histograms.keys()):
            
            tree_name = histograms[hist_name][0]
            label_x = histograms[hist_name][1]
            label_y = histograms[hist_name][2]  
            xmin = histograms[hist_name][3][0]
            xmax = histograms[hist_name][3][1]
            ymin = histograms[hist_name][4][0]
            ymax = histograms[hist_name][4][1]

            file = ROOT.TFile.Open(wjet_path + hist_name, 'read')
            hist_data = file.Get(tree_name)
            hist_data.SetDirectory(0)
            hist_data.SetTitle("W(#mu#nu)+jets fake rate")
            # hist_data.Print("all")
            
            #Set axis Titles and user ranges
            hist_data.GetXaxis().SetTitle(label_x)
            hist_data.GetYaxis().SetTitle(label_y)
            # hist_data.GetXaxis().SetRangeUser(xmin, xmax)
            # hist_data.GetYaxis().SetRangeUser(ymin, ymax)
            
            canvas = ROOT.TCanvas("canvas", "W(#mu#nu)+jets fake rate", 1200, 1000)
            
            # remove stat window from canvas
            ROOT.gStyle.SetOptStat(0)
            
            # Shift axis title closer to the plot
            hist_data.GetXaxis().SetTitleOffset(1.4)
            hist_data.GetYaxis().SetTitleOffset(0.9)

            # Set smaller font size for axis titles
            hist_data.GetXaxis().SetTitleSize(0.04)
            hist_data.GetYaxis().SetTitleSize(0.04)
            
            ROOT.gStyle.SetPaintTextFormat(".3f")

            # Adjust margins to ensure labels are visible
            canvas.SetLeftMargin(0.10)
            canvas.SetRightMargin(0.15)
            canvas.SetBottomMargin(0.12)
            
            canvas.SetLogx(1)
            canvas.SetLogy(1)
            
            
            # Draw the histogram with COLZ, TEXT45, and E options
            hist_data.Draw("COLZ TEXT45 E")
            
            # Update the canvas to reflect changes
            canvas.Update()

            # Save the canvas to a file
            canvas.SaveAs(output_dir+f"hist_name{period_names[period]}.pdf")
            