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

parser = ArgumentParser(
    description="Plot histograms from previously created scale factors to compare over different flavour id"
    "to run: ")
parser.add_argument(
    "--hist1", nargs='+', help="Histograms over the flavor: fake_rate_score_jet_score_jet_pt_jet_dxy.root")
parser.add_argument(
    "--hist2", nargs='+', help="Histogram inclusive fake_rate_score_jet_score_jet_pt_jet_dxy.root")
parser.add_argument(
    "--outdir", help="Output directory. If not given, output to the directory "
    "where histfile is located", default='hist_output')
parser.add_argument(
    "--ext", choices=["pdf", "svg", "png"], help="Output file format",
    default="png")


args = parser.parse_args()

try: os.stat(args.outdir)
except: os.mkdir(args.outdir)

if True:
    # output_name = "fake_jet_tau_"
    # tree_name = "fake_rate_score_jet_score_jet_pt_jet_dxy"
    # # axis: ["jet_score", "jet_pt", "jet_dxy"]
    # jet_score_bins = [7, 8, 9, 10] 
    # # [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 2.0, 4.0, 10.0, 16.0, 20.0, 30.0, 50.0],
    # # jet_dxy_bin = []
    # # [20, 30, 40, 50, 70, 120, 1000],
    # jet_pt_bins = [1, 2, 3, 4, 5]
    
    # names_score = [0.1716, 0.4612, 0.786, 0.846, 0.9872, 0.9900, 0.9972, 0.9996]
    # names_pt = [20, 30, 40, 50, 70, 120, 1000]
    # names_dxy = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 2.0, 4.0, 10.0, 16.0, 20.0, 30.0, 50.0]
    # names_flav = [1, 4, 5, 6, 15, 16, 21, 22, 25]
    names_flav = ["u/d/s", "c/b", 6, "tau", 16, "g", 22]
    #"jet_flav" : [1, 4, 6, 15, 16, 21, 22],

    output_name = "fake_WJets_gen_vs_data_dxy_"
    # tree_name = "fake_rate_jet_prob_dxy"
    # jet_score_bins = [2,3,4,5,6,7,8] 
    # jet_pt_bins = [1,2,3,4,5,6,7]
    # jet_flav_bins = [1, 2, 3, 5, 7]
    jet_flav_bins = [1, 2, 4, 6]

    # tree_name = "fake_rate_jet_flav_jet_pt_jet_dxy"
    # tree_name2  = "fake_rate_eta_jet_pt_jet_dxy"
    tree_name = "fake_rate_jet_flav_dxy"
    # tree_name2 = "fake_rate_jet_eta_pt"
    tree_name2 = "fake_rate_jet_dxy"

    # print(list(itertools.product(jet_score_bins, jet_pt_bins)))

    # for jet_score_bin, jet_pt_bin in list(itertools.product(jet_score_bins, jet_pt_bins)):
    # for jet_score_bin in jet_score_bins:
    hist_data = []
    for i, jet_flav_bin in enumerate(jet_flav_bins):
        # for jet_pt_bin in jet_pt_bins:
        file = ROOT.TFile.Open(str(args.hist1[1]), 'read')
        # hist_data = file.Get(data).ProjectionY("hist1",jet_score_bin, jet_score_bin,
        #                                        jet_dxy_bin, jet_dxy_bin)
        # hist_data = file.Get(tree_name).ProjectionZ("hist1",jet_flav_bin, jet_flav_bin,
        #                                                     jet_pt_bin, jet_pt_bin)
        # hist_data = file.Get(tree_name).ProjectionX("hist1",jet_flav_bin, jet_flav_bin)
        # hist_data.SetDirectory(0)
        # hist_data.SetMarkerStyle(22)
        # hist_data.SetMarkerSize(1)
        # hist_data.SetMarkerColor(46)
        # hist_data.SetLineColor(46)
        # hist_data.SetFillColor(0)
        # hist_data.SetLineWidth(2)
        # hist_data.SetTitle(args.hist1[0])
        # print(hist_data)
        hist_data.append(file.Get(tree_name).ProjectionX("hist1",jet_flav_bin, jet_flav_bin))
        hist_data[-1].SetDirectory(0)
        hist_data[-1].SetMarkerStyle(22)
        hist_data[-1].SetMarkerSize(1)
        hist_data[-1].SetMarkerColor(ColorIterator(i, 10))
        hist_data[-1].SetLineColor(ColorIterator(i, 10))
        hist_data[-1].SetFillColor(0)
        hist_data[-1].SetLineWidth(2)
        # hist_data[-1].SetTitle(args.hist1[0])
        flav_bin_name = str(names_flav[jet_flav_bin-1]).replace(".", "p")
        hist_data[-1].SetTitle(f"flav={flav_bin_name}")

    file_data = ROOT.TFile.Open(str(args.hist2[1]), 'read')
    # hist_prediction = file_data.Get(tree_name2).ProjectionZ("hist2",1, 1,
                                                            # jet_pt_bin, jet_pt_bin)
    # hist_prediction = file_data.Get(tree_name2).ProjectionY("hist2",jet_score_bin, jet_score_bin)
    # hist_prediction = file_data.Get(tree_name2).ProjectionX("hist2",1, 1)
    hist_prediction = file_data.Get(tree_name2)
    hist_prediction.SetDirectory(0)
    hist_prediction.SetMarkerStyle(21)
    hist_prediction.SetMarkerSize(1)
    hist_prediction.SetMarkerColor(30)
    hist_prediction.SetLineColor(30)
    hist_prediction.SetFillColor(0)
    hist_prediction.SetLineWidth(2)
    hist_prediction.SetTitle(args.hist2[0])
        # print(hist_prediction.Print("all"))
        # print(hist_prediction)

        # bin_names = ["", "" ,"", ">0.1716", ">0.4612", ">0.6631", ">0.786", ">0.846", ">0.9529", ">0.9889", ">0.9972", ">0.9996"]
        # for hist in [hist_data, hist_prediction]:
        #     for bin_i, label in enumerate(bin_names):
        #         hist.GetXaxis().SetBinLabel(bin_i, label)
        #         hist.GetXaxis().SetTitle("")

        # score_bin_name = str(names_score[jet_score_bin-1]).replace(".", "p")
    # flav_bin_name = str(names_flav[jet_flav_bin-1]).replace(".", "p")
    # pt_bin_name = str(names_pt[jet_pt_bin-1]).replace(".", "p")
    flav_bin_name = "all"
    pt_bin_name = "incl"
    root_plot1D(
        l_hist = [hist_prediction],
        l_hist_overlay = hist_data,
        outfile = args.outdir + "/" + output_name + f'flav{flav_bin_name}_pt{pt_bin_name}.' + args.ext,
        xrange = [hist_prediction.GetXaxis().GetXmin(), hist_prediction.GetXaxis().GetXmax()],
        yrange = (0.001,  10*max(hist_prediction.GetMaximum(), hist_prediction.GetMaximum())),
        # yrange = (0.0, 0.2),
        logx = True, logy = True,
        logx_ratio = True, logy_ratio = False,
        include_overflow = False,
        # xtitle = "Jet pt (GeV)",
        xtitle = "Jet dxy (cm)",
        ytitle = "Fake rate",
        # xtitle_ratio = "Jet pt (GeV)",
        xtitle_ratio = "Jet dxy (cm)",
        ytitle_ratio = "ratio",
        centertitlex = True, centertitley = True,
        centerlabelx = False, centerlabely = False,
        gridx = True, gridy = True,
        ndivisionsx = None,
        stackdrawopt = "nostack",
        # normilize = True,
        normilize_overlay = False,
        legendpos = "UR",
        legendtitle = f"",
        legendncol = 2,
        legendtextsize = 0.055,
        legendwidthscale = 1.0,
        legendheightscale = 3.0,
        lumiText = "2018 (13 TeV)",
        signal_to_background_ratio = True,
        ratio_mode = "DATA",
        yrange_ratio = (0.0, 2),
        ndivisionsy_ratio = (4, 5, 0),
        # yrange_ratio = (0.8, 1.2),
        draw_errors = True
    )