
import ROOT
ROOT.gROOT.SetBatch(True)

import os
import json
import numpy as np
from argparse import ArgumentParser

from pepper import Config

from utils.utils import ColorIterator, root_plot1D

parser = ArgumentParser(
    description="Plot histograms from previously created histograms")
parser.add_argument(
    "plot_config", help="Path to a configuration file for plotting")
parser.add_argument(
    "histfile", nargs="+", help="Coffea file with a single histogram or a "
    "JSON file containing histogram info. See output of select_events.py")
parser.add_argument(
    "--outdir", help="Output directory. If not given, output to the directory "
    "where histfile is located")
parser.add_argument(
    "--cutflow", help="cutflow file")
parser.add_argument(
    "--ext", choices=["pdf", "svg", "png"], help="Output file format",
    default="svg")
parser.add_argument(
    "-s", "--signals", nargs='*', default=["None"], help="Set of signal "
    "points to plot. Can be All, None (default) or a name (or series of names)"
    " of a specific set defined in the config")

args = parser.parse_args()

try: os.stat(args.outdir)
except: os.mkdir(args.outdir)

config = Config(args.plot_config)

try:
    with open(config["crosssections"]) as f:
        crosssections = json.load(f)
except:
    raise ValueError('Can not open file with crosssections')

try:
    with open(args.cutflow) as f:
        cutflow = json.load(f)
except:
    raise ValueError('Can not open file with cutflow')

histfiles = []
for histfile in args.histfile:
    if histfile.endswith(".json"):
        dirname = os.path.dirname(histfile)
        with open(histfile) as f:
            f = json.load(f)
            for keys, histfile in zip(*f):
                if len(keys) != 2:
                    continue
                histfiles.append(os.path.join(dirname, histfile))
    else:
        histfiles.append(histfile)


for _histfile in histfiles:

    file = ROOT.TFile.Open(str(_histfile), 'read')

    _histograms = {"background":[], "signal":[]}
    for _group_idx, _group_name in enumerate(config["Labels"].keys()):

        isSignal = "signal" if _group_name in config["Signal_samples"] else "background"
        
        for _idx, _histogram_name in enumerate(config["Labels"][_group_name]):
            
            # rescaling:
            hist = file.Get(_histogram_name)
            N = cutflow[_histogram_name]["all"]["Before cuts"]
            hist.Scale( (crosssections[_histogram_name] * config["luminosity"]) / N)
            
            if _idx == 0:
                _histograms[isSignal].append(hist)
            else:
                _histograms[isSignal][-1].Add(hist)


        _histograms[isSignal][-1].SetLineColor(ColorIterator(_group_idx))
        _histograms[isSignal][-1].SetFillColor(ColorIterator(_group_idx))
        _histograms[isSignal][-1].SetLineWidth(2)
        _histograms[isSignal][-1].SetMarkerSize(0)
        _histograms[isSignal][-1].SetTitle(_group_name)
    
    # get maximum for the y-scale
    y_max = _histograms["background"][0].GetMaximum()
    for _h in _histograms["background"]:
        y_max = max(y_max,_h.GetMaximum())

    # sort histogram from min to max
    _histograms_background_entries = []
    _histograms_background_sorted = []
    for _h in _histograms["background"]:
        _histograms_background_entries.append(_h.Integral())
    _sorted_hist = np.argsort(_histograms_background_entries)
    for _idx in _sorted_hist:
        _histograms_background_sorted.append(_histograms["background"][_idx])

    root_plot1D(
        l_hist = _histograms_background_sorted,
        l_hist_overlay = _histograms["signal"],
        outfile = args.outdir + "/" + os.path.splitext(os.path.basename(_histfile))[0] + ".png",
        xrange = [_histograms["signal"][0].GetXaxis().GetXmin(), 
                  _histograms["signal"][0].GetXaxis().GetXmin()],
        yrange = (0.1,  80*y_max),
        logx = False, logy = True,
        ytitle = _histograms["signal"][0].GetYaxis().GetTitle(),
        xtitle = _histograms["signal"][0].GetXaxis().GetTitle(),
        centertitlex = True, centertitley = True,
        centerlabelx = False, centerlabely = False,
        gridx = True, gridy = True,
        ndivisionsx = None,
        stackdrawopt = "",
        legendpos = "UL",
        legendtitle = f"",
        legendncol = 2,
        legendtextsize = 0.04,
        legendwidthscale = 1.5,
        legendheightscale = 0.5,
        lumiText = "2018 (13 TeV)"
    )