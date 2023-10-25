import os
import json
from argparse import ArgumentParser
import re

from pepper import Config
from utils.plotter import plot1D, plot2D, plot_predict


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
parser.add_argument(
    '-m','--mode', nargs='+', default=[], help="The set of histograms to be plotted "
    "Available options: [1D, 2D]", required=True)
parser.add_argument(
    '-d','--data', action='store_true', help="If True Data/MC comparison will be plotted"
        "otherwise signal/background will be plotted")


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

if "1D" in args.mode:
    histfiles = []
    histnames = []
    for histfile in args.histfile:
        if histfile.endswith(".json"):
            dirname = os.path.dirname(histfile)
            with open(histfile) as f:
                f = json.load(f)
                for keys, histfile in zip(*f):
                    if len(keys) != 2:
                        continue
                    # 2D histogram will be filtered out (labeled by TH2)
                    if re.search(".*_TH2", keys[1]):
                        continue
                    if keys[1] in config["mask_hists"]:
                        continue
                    histfiles.append(os.path.join(dirname, histfile))
                    histnames.append(keys[1])
        else:
            raise ValueError('Json should be provided')

    plot1D(histfiles, histnames, config, crosssections, cutflow, args.outdir, args.data)

if "2D" in args.mode:
    histfiles = []
    histnames = []
    for histfile in args.histfile:
        if histfile.endswith(".json"):
            dirname = os.path.dirname(histfile)
            with open(histfile) as f:
                f = json.load(f)
                for keys, histfile in zip(*f):
                    if len(keys) != 2:
                        continue
                    # 2D histogram will be filtered out (labeled by TH2)
                    if not re.search(".*_TH2", keys[1]):
                        continue
                    # if not re.search(".*mt2_sum_mt_TH2", keys[1]):
                    #     continue
                    histfiles.append(os.path.join(dirname, histfile))
                    histnames.append(keys[1])
        else:
            raise ValueError('Json should be provided')

    plot2D(histfiles, histnames, config, crosssections, cutflow, args.outdir)
    
if "prediction" in args.mode:
    if not args.histfile[0].endswith(".json"):
        raise ValueError('Json should be provided')
    dirname = os.path.dirname(args.histfile[0])
    plot_predict(dirname, config, crosssections, cutflow, args.outdir)
