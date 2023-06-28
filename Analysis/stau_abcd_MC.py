import os
import json
from argparse import ArgumentParser
import re
import ROOT
ROOT.gROOT.SetBatch(True)
import itertools
import numpy as np

from pepper import Config
# from utils.plotter import plot1D, plot2D
from utils.plotter  import ColorIterator, root_plot1D, root_plot2D

## Two fakes and fake genuine

# region_contruct = {  
#      "A" : [ ["M3"], ["F1F2", "F1G2","G1F2"], ["DDtag"] ],
#      "B" : [ ["M1","M2"], ["F1F2", "F1G2","G1F2"], ["DDtag"] ],
#      "C" : [ ["M3"], ["F1F2", "F1G2","G1F2"], ["DDtag_inv"] ],
#      "D" : [ ["M1","M2"], ["F1F2", "F1G2","G1F2"], ["DDtag_inv"] ],
#     }

# region_contruct = {  
#      "A" : [ ["M3"], ["F1F2", "F1G2","G1F2"], ["PDtag"] ],
#      "B" : [ ["M1","M2"], ["F1F2", "F1G2","G1F2"], ["PDtag"] ],
#      "C" : [ ["M3"], ["F1F2", "F1G2","G1F2"], ["PDtag_inv"] ],
#      "D" : [ ["M1","M2"], ["F1F2", "F1G2","G1F2"], ["PDtag_inv"] ],
#     }

# region_contruct = {  
#      "A" : [ ["M3"], ["F1F2", "F1G2","G1F2"], ["PP"] ],
#      "B" : [ ["M1","M2"], ["F1F2", "F1G2","G1F2"], ["PP"] ],
#      "C" : [ ["M3"], ["F1F2", "F1G2","G1F2"], ["PP_inv"] ],
#      "D" : [ ["M1","M2"], ["F1F2", "F1G2","G1F2"], ["PP_inv"] ],
#     }

## Two fakes

# region_contruct = {  
#      "A" : [ ["M3"], ["F1F2"], ["DDtag"] ],
#      "B" : [ ["M1","M2"], ["F1F2"], ["DDtag"] ],
#      "C" : [ ["M3"], ["F1F2"], ["DDtag_inv"] ],
#      "D" : [ ["M1","M2"], ["F1F2"], ["DDtag_inv"] ],
#     }

# region_contruct = {  
#      "A" : [ ["M3"], ["F1F2"], ["PDtag"] ],
#      "B" : [ ["M1","M2"], ["F1F2"], ["PDtag"] ],
#      "C" : [ ["M3"], ["F1F2"], ["PDtag_inv"] ],
#      "D" : [ ["M1","M2"], ["F1F2"], ["PDtag_inv"] ],
#     }

# region_contruct = {  
#      "A" : [ ["M3"], ["F1F2"], ["PP"] ],
#      "B" : [ ["M1","M2"], ["F1F2"], ["PP"] ],
#      "C" : [ ["M3"], ["F1F2"], ["PP_inv"] ],
#      "D" : [ ["M1","M2"], ["F1F2"], ["PP_inv"] ],
#     }

## Two genuine

# region_contruct = {  
#      "A" : [ ["M3"], ["G1G2"], ["DDtag"] ],
#      "B" : [ ["M1","M2"], ["G1G2"], ["DDtag"] ],
#      "C" : [ ["M3"], ["G1G2"], ["DDtag_inv"] ],
#      "D" : [ ["M1","M2"], ["G1G2"], ["DDtag_inv"] ],
#     }

region_contruct = {  
     "A" : [ ["M3"], ["G1G2"], ["PDtag"] ],
     "B" : [ ["M1","M2"], ["G1G2"], ["PDtag"] ],
     "C" : [ ["M3"], ["G1G2"], ["PDtag_inv"] ],
     "D" : [ ["M1","M2"], ["G1G2"], ["PDtag_inv"] ],
    }

region_contruct = {  
     "A" : [ ["M3"], ["G1G2"], ["PP"] ],
     "B" : [ ["M1","M2"], ["G1G2"], ["PP"] ],
     "C" : [ ["M3"], ["G1G2"], ["PP_inv"] ],
     "D" : [ ["M1","M2"], ["G1G2"], ["PP_inv"] ],
    }

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

output = args.outdir + "/" + region_contruct["A"][-1][0]

for _histfile, _histname in zip(histfiles,histnames):
        
    if not any([cut in str(_histfile) for cut in config["cuts"]]):
        continue
    
    # sum up histograms for each abcd-region     
    _histograms_predict = {"A":[], "B":[], "C":[], "D":[]}
    _histograms_stack = {"A":[], "B":[], "C":[], "D":[]}
    
    for region_name in ["A", "B", "C", "D"]:
    
        # The plots are done for specific cut step specified in config.
        file = ROOT.TFile.Open(str(_histfile), 'read')
        for _group_idx, _group_name in enumerate(config["Labels"].keys()):
            
            if _group_name in config["Signal_samples"] or _group_name in config["Data"]:
                continue
            
            # Accumulate histograms amoung categories accroding to region_contruct
            categories_list = list(itertools.product(*region_contruct[region_name]))
            categories_list = ["_".join(cat) for cat in categories_list]
            
            for _ci, _categ in enumerate(categories_list):
                
                # Accumulate the dataset for the particular data group as specified in config “Labels”.
                for _idx, _histogram_data in enumerate(config["Labels"][_group_name]):
                    
                    # Rescaling according to cross-section and luminosity
                    print("Reading data:", _histogram_data + "_" + _categ)
                    hist = file.Get(_histogram_data + "_" + _categ)
                    hist.SetDirectory(0)
                    N = cutflow[_histogram_data]["all"]["Before cuts"]
                    hist.Scale( (crosssections[_histogram_data] * config["luminosity"]) / N)

                    if _histname in config["SetupBins"]:
                        hist.Rebin(config["SetupBins"][_histname][2])
                        if config["SetupBins"][_histname][4]:
                            for bin_i, label in enumerate(config["SetupBins"][_histname][4]):
                                hist.GetXaxis().SetBinLabel(bin_i, label)
                                hist.GetXaxis().SetTitle("")

                    hist_copy = hist.Clone()
                    hist_copy.SetDirectory(0)
                    # To have histogram for every group
                    if _ci==0 and _idx == 0:
                        _histograms_stack[region_name].append(hist_copy)
                    else:
                        _histograms_stack[region_name][-1].Add(hist_copy)
                    
                    hist_copy = hist.Clone()
                    hist_copy.SetDirectory(0)
                    # To sum up histograms for all groups
                    if _group_idx==0 and _ci==0 and _idx == 0:
                        _histograms_predict[region_name].append(hist_copy)
                    else:
                        _histograms_predict[region_name][-1].Add(hist_copy)
                
            color_setup = config["MC_bkgd"][_group_name]  
            line_color = color_setup[1]
            fill_color = color_setup[0]
            
            _histograms_stack[region_name][-1].SetMarkerSize(0)
            _histograms_stack[region_name][-1].SetLineWidth(4)
            _histograms_stack[region_name][-1].SetLineColor(line_color)
            _histograms_stack[region_name][-1].SetFillColor(fill_color)
            _histograms_stack[region_name][-1].SetTitle(_group_name)
            
        _histograms_predict[region_name][-1].SetMarkerStyle(8)
        _histograms_predict[region_name][-1].SetMarkerSize(2)
        _histograms_predict[region_name][-1].SetLineWidth(0)
        _histograms_predict[region_name][-1].SetFillColor(0)
        _histograms_predict[region_name][-1].SetTitle("A=B/C*D")
                
    # get maximum for the y-scale
    y_max = _histograms_stack["A"][0].GetMaximum()
    for _h in _histograms_stack["A"]:
        y_max = max(y_max,_h.GetMaximum())

    # sort histogram from min to max
    _histograms_background_entries = []
    _histograms_background_sorted = []
    for _h in _histograms_stack["A"]:
        _histograms_background_entries.append(_h.Integral())
    _sorted_hist = np.argsort(_histograms_background_entries)
    for _idx in _sorted_hist:
        _histograms_background_sorted.append( _histograms_stack["A"][_idx])

    # read the binning if available:
    if _histname in config["SetupBins"]:
        xrange_min = config["SetupBins"][_histname][0]
        xrange_max = config["SetupBins"][_histname][1]
        overflow =  bool(config["SetupBins"][_histname][3])
    else:
        xrange_min = _histograms_stack["A"][0].GetXaxis().GetXmin()
        xrange_max = _histograms_stack["A"][0].GetXaxis().GetXmax()
        overflow =  True
    
    
    # ABCD algebra
    _histograms_predict["A"][0] = _histograms_predict["C"][0].Divide(_histograms_predict["D"][0])*_histograms_predict["B"][0]
    
    root_plot1D(
            l_hist = _histograms_background_sorted,
            l_hist_overlay = _histograms_predict["A"],
            outfile = output + "/" + os.path.splitext(os.path.basename(_histfile))[0] + ".png",
            xrange = [xrange_min, xrange_max],
            yrange = (0.001,  10000*y_max),
            logx = False, logy = True,
            logx_ratio = False, logy_ratio = False,
            include_overflow = overflow,
            xtitle = _histograms_stack["A"][0].GetXaxis().GetTitle(),
            ytitle = "events",
            xtitle_ratio = _histograms_stack["A"][0].GetXaxis().GetTitle(),
            ytitle_ratio = "MC_{pred}/MC",
            centertitlex = True, centertitley = True,
            centerlabelx = False, centerlabely = False,
            gridx = True, gridy = True,
            ndivisionsx = None,
            stackdrawopt = "",
            # normilize = True,
            normilize_overlay = False,
            legendpos = "UR",
            legendtitle = f"",
            legendncol = 3,
            legendtextsize = 0.025,
            legendwidthscale = 1.9,
            legendheightscale = 0.26,
            lumiText = "2018 (13 TeV)",
            signal_to_background_ratio = True,
            ratio_mode = "DATA",
            yrange_ratio = (0.0, 5.0),
            draw_errors = True
            )
    
    # exit()