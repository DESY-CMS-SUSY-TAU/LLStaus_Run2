
from turtle import color
import ROOT
ROOT.gROOT.SetBatch(True)

import os
import json
import numpy as np
from argparse import ArgumentParser

from pepper import Config

from utils.utils import ColorIterator, root_plot1D

parser = ArgumentParser(
    description="Plot histogram ratio from previously created histograms.")
parser.add_argument(
    "plot_config", help="Path to a configuration file for plotting")
parser.add_argument(
    "--outdir", help="Output directory. If not given, output to the directory "
    "where histfile is located")
parser.add_argument(
    "--cutflow", help="cutflow file")

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
        basedir = os.path.dirname(args.cutflow)+"/hists/"
except:
    raise ValueError('Error reading/open file with cutflow')

# histfiles = []
# histnames = []
# for histfile in args.histfile:
#     if histfile.endswith(".json"):
#         dirname = os.path.dirname(histfile)
#         with open(histfile) as f:
#             f = json.load(f)
#             for keys, histfile in zip(*f):
#                 if len(keys) != 2:
#                     continue
#                 histfiles.append(os.path.join(dirname, histfile))
#                 histnames.append(keys[1])
#     else:
#         raise ValueError('Json should be provided')

# print(basedir)
# print(config["Histograms_ratio"])
# exit()

# for _histfile, _histname in zip(histfiles,histnames):
for _histogram_setup in config["Histograms_ratio"]:
    for _datagroup in config["Labels"]:
        
        print(_histogram_setup, _datagroup)
        
        _dataset = config["Labels"][_datagroup]
        _histograms = []
        _base_histogram = []
        for _idx, _filename in enumerate(config["Histograms_ratio"][_histogram_setup]):
            files_setup = config["Histograms_ratio"][_histogram_setup][_filename]
            file = ROOT.TFile.Open(basedir + str(_filename), 'read')

            hist = file.Get(_dataset)
            hist.SetDirectory(0)
            hist.Rebin(2)
            hist_put = hist.Clone()
            hist_put.SetDirectory(0)
            
            if files_setup[0] == "base":
                _base_histogram.append(hist)
            else:
                if not _base_histogram:
                    raise Exception("Error: Base histogram should be first.")

            n = hist_put.Integral()
            color = ColorIterator(_idx, 2)
            hist_put.Divide(_base_histogram[-1])
            hist_put.SetTitle(files_setup[1]+ f'(evnt:{n})')
            hist_put.SetLineColor(color)
            hist_put.SetLineWidth(2)
            hist_put.SetMarkerSize(0)
            _histograms.append(hist_put)
            file.Close()

        xrange_min = _base_histogram[0].GetXaxis().GetXmin()
        xrange_max = _base_histogram[0].GetXaxis().GetXmax()

        root_plot1D(
            l_hist = _histograms,
            outfile = args.outdir + "/" + _histogram_setup + "_" + _dataset + ".png",
            xrange = [0.0, 50.0],
            yrange = (0.0, 1.6),
            # yrange = (1,  1.1*y_max),
            logx = False, logy = False,
            include_overflow = False,
            # ytitle = _histograms["signal"][0].GetYaxis().GetTitle(),
            xtitle = _base_histogram[0].GetXaxis().GetTitle(),
            ytitle = "eff.",
            centertitlex = True, centertitley = True,
            centerlabelx = False, centerlabely = False,
            gridx = True, gridy = True,
            ndivisionsx = None,
            stackdrawopt = "nostack",
            # normilize = True,
            legendpos = "UR",
            legendtitle = f"",
            legendncol = 1,
            legendtextsize = 0.03,
            legendwidthscale = 1.5,
            legendheightscale = 0.8,
            lumiText = "2018 (13 TeV)"
        )


    # file = ROOT.TFile.Open(str(_histfile), 'read')

    # _histograms = {"background":[], "signal":[]}
    # for _group_idx, _group_name in enumerate(config["Labels"].keys()):

    #     isSignal = "signal" if _group_name in config["Signal_samples"] else "background"
        
    #     for _idx, _histogram_name in enumerate(config["Labels"][_group_name]):
            
    #         # rescaling:
    #         hist = file.Get(_histogram_name)
    #         N = cutflow[_histogram_name]["all"]["Before cuts"]
    #         hist.Scale( (crosssections[_histogram_name] * config["luminosity"]) / N)

    #         if _histname in config["SetupBins"]:
    #             hist.Rebin(config["SetupBins"][_histname][2])
            
    #         if _idx == 0:
    #             _histograms[isSignal].append(hist)
    #         else:
    #             _histograms[isSignal][-1].Add(hist)

    #     if isSignal == "signal":
    #         color_setup = config["Signal_samples"][_group_name]  
    #         line_color = color_setup[1]
    #         fill_color = color_setup[0]
    #         _histograms[isSignal][-1].SetLineStyle(2)
    #     else:
    #         color_setup = config["MC_bkgd"][_group_name]  
    #         line_color = color_setup[1]
    #         fill_color = color_setup[0]
        
    #     _histograms[isSignal][-1].SetLineColor(line_color)
    #     _histograms[isSignal][-1].SetFillColor(fill_color)
    #     _histograms[isSignal][-1].SetLineWidth(2)
    #     _histograms[isSignal][-1].SetMarkerSize(0)
    #     _histograms[isSignal][-1].SetTitle(_group_name)
    
    # # get maximum for the y-scale
    # y_max = _histograms["background"][0].GetMaximum()
    # for _h in _histograms["background"]:
    #     y_max = max(y_max,_h.GetMaximum())

    # # sort histogram from min to max
    # _histograms_background_entries = []
    # _histograms_background_sorted = []
    # for _h in _histograms["background"]:
    #     _histograms_background_entries.append(_h.Integral())
    # _sorted_hist = np.argsort(_histograms_background_entries)
    # for _idx in _sorted_hist:
    #     _histograms_background_sorted.append(_histograms["background"][_idx])

    # # read the binning if available:

    # if _histname in config["SetupBins"]:
    #     xrange_min = config["SetupBins"][_histname][0]
    #     xrange_max = config["SetupBins"][_histname][1]
    #     overflow =  bool(config["SetupBins"][_histname][3])
    #     # print(overflow)
    # else:
    #     xrange_min = _histograms["signal"][0].GetXaxis().GetXmin()
    #     xrange_max = _histograms["signal"][0].GetXaxis().GetXmax()
    #     overflow =  True


    # root_plot1D(
    #     l_hist = _histograms_background_sorted,
    #     l_hist_overlay = _histograms["signal"],
    #     outfile = args.outdir + "/" + os.path.splitext(os.path.basename(_histfile))[0] + ".png",
    #     xrange = [xrange_min, xrange_max],
    #     yrange = (0.1,  100*y_max),
    #     # yrange = (1,  1.1*y_max),
    #     logx = False, logy = True,
    #     include_overflow = overflow,
    #     # ytitle = _histograms["signal"][0].GetYaxis().GetTitle(),
    #     xtitle = _histograms["signal"][0].GetXaxis().GetTitle(),
    #     ytitle = "arb. unit",
    #     centertitlex = True, centertitley = True,
    #     centerlabelx = False, centerlabely = False,
    #     gridx = True, gridy = True,
    #     ndivisionsx = None,
    #     stackdrawopt = "",
    #     # normilize = True,
    #     legendpos = "UR",
    #     legendtitle = f"",
    #     legendncol = 3,
    #     legendtextsize = 0.03,
    #     legendwidthscale = 1.9,
    #     legendheightscale = 0.4,
    #     lumiText = "2018 (13 TeV)"
    # )