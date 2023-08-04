import os
import json
from argparse import ArgumentParser
import re
import ROOT

from pepper import Config
from utils.plotter import plot1D, plot2D
from utils.utils import ColorIterator, root_plot1D, root_plot2D

parser = ArgumentParser(
    description="Plot histograms from previously created histograms")
parser.add_argument(
    "plot_config", help="Path to a configuration file for plotting")
parser.add_argument(
    "histfile", nargs="+", help="Coffea file with a single histogram or a "
    "JSON file containing histogram info. See output of select_events.py")
parser.add_argument(
    "--outdir", help="Output directory. If not given, output to the directory "
    "where histfile is located", default='hist_output')
parser.add_argument(
    "--cutflow", help="cutflow file", required=True)
parser.add_argument(
    "--ext", choices=["pdf", "svg", "png"], help="Output file format",
    default="svg")


args = parser.parse_args()

try: os.stat(args.outdir)
except: os.mkdir(args.outdir)

config = Config(args.plot_config)

try:
    with open(config["crosssections"]) as f:
        crosssections = json.load(f)
except:
    raise ValueError('Error reading/open file with crosssections')

try:
    with open(args.cutflow) as f:
        cutflow = json.load(f)
except:
    raise ValueError('Error reading/open file with cutflow')

if True:

    # ---------------------------
    # data_TH2 = os.path.dirname(args.histfile[0]) + "/Cut-007_charge_n_jet_pass.root"
    # data_TH2_bin = 4
    # predicted_TH1 = os.path.dirname(args.histfile[0]) + "/Cut-007_charge_yield_bin0to2.root"
    # data = "SingleMuon"
    # output_name = "data_preditcion_bin0to2.png"
    #---------------------------
    data_TH2 = os.path.dirname(args.histfile[0]) + "/Cut-007_charge_n_jet_pass.root"
    data_TH2_bin = 3
    predicted_TH1 = os.path.dirname(args.histfile[0]) + "/Cut-007_charge_yield_bin0to1.root"
    data = "SingleMuon"
    output_name = "data_preditcion_bin0to1.png"
    #---------------------------
    # data_TH2 = os.path.dirname(args.histfile[0]) + "/Cut-007_charge_n_jet_pass.root"
    # data_TH2_bin = 4
    # predicted_TH1 = os.path.dirname(args.histfile[0]) + "/Cut-007_charge_yield_bin1to2.root"
    # data = "SingleMuon"
    # output_name = "data_preditcion_bin1to2.png"
    #---------------------------

    file = ROOT.TFile.Open(str(data_TH2), 'read')
    hist_data = file.Get(data).ProjectionX("SingleMuon pred.",data_TH2_bin,data_TH2_bin).Clone("SingleMuon pred.")
    hist_data.SetDirectory(0)
    # print(hist_prediction.Print("all"))
    hist_data.SetDirectory(0)
    hist_data.SetMarkerStyle(8)
    hist_data.SetMarkerSize(1)
    hist_data.SetMarkerColor(1)
    hist_data.SetLineWidth(2)
    hist_data.SetTitle("Data")
    print(hist_data)

    file_data = ROOT.TFile.Open(str(predicted_TH1), 'read')
    hist_prediction = file_data.Get(data).Clone("SingleMuon")
    hist_prediction.SetMarkerStyle(21)
    hist_prediction.SetMarkerColor(30)
    hist_prediction.SetLineColor(30)
    hist_prediction.SetFillColor(0)
    hist_prediction.SetLineWidth(2)
    hist_prediction.SetTitle("Pred.")
    # print(hist_prediction.Print("all"))
    print(hist_prediction)

    bin_names = ["", "" ,"", ">0.1716", ">0.4612", ">0.6631", ">0.786", ">0.846", ">0.9529", ">0.9889", ">0.9972", ">0.9996"]
    for hist in [hist_data, hist_prediction]:
        for bin_i, label in enumerate(bin_names):
            hist.GetXaxis().SetBinLabel(bin_i, label)
            hist.GetXaxis().SetTitle("")

    root_plot1D(
        l_hist = [hist_prediction],
        l_hist_overlay = [hist_data],
        outfile = args.outdir + "/" + output_name,
        xrange = [1, 10],
        yrange = (1.0,  10*hist_data.GetMaximum()),
        logx = False, logy = True,
        logx_ratio = False, logy_ratio = False,
        include_overflow = False,
        xtitle = "tag bins",
        ytitle = "events",
        xtitle_ratio = "tag bins",
        ytitle_ratio = "DATA/MC",
        centertitlex = True, centertitley = True,
        centerlabelx = False, centerlabely = False,
        gridx = True, gridy = True,
        ndivisionsx = None,
        stackdrawopt = "",
        # normilize = True,
        normilize_overlay = False,
        legendpos = "UR",
        legendtitle = f"",
        legendncol = 2,
        legendtextsize = 0.055,
        legendwidthscale = 1.0,
        legendheightscale = 1.0,
        lumiText = "2018 (13 TeV)",
        signal_to_background_ratio = True,
        ratio_mode = "DATA",
        yrange_ratio = (0.0, 2.5),
        # yrange_ratio = (0.8, 1.2),
        draw_errors = True
    )