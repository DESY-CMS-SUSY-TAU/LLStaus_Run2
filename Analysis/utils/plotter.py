import ROOT
ROOT.gROOT.SetBatch(True)
import itertools
import numpy as np
import os

from .utils import ColorIterator, root_plot1D, root_plot2D

def plot1D(histfiles, histnames, config, xsec, cutflow, output_path, isData):

    categories_list = list(itertools.product(*config["Categories"]))
    # categories_list = [f"{cat1}_{cat2}_{cat3}" for cat1,cat2,cat3 in categories_list]
    categories_list = ["_".join(cat) for cat in categories_list]
    categories_list = [""]

    for _ci, _categ in enumerate(categories_list):

        output = output_path + "/" + _categ

        for _histfile, _histname in zip(histfiles,histnames):
            
            # The plots are done for specific cut step specified in config.
            if not any([cut in str(_histfile) for cut in config["cuts"]]):
                continue

            print("OPEN:", _histfile)
            file = ROOT.TFile.Open(str(_histfile), 'read')

            _histograms = {"background":[], "signal":[], "data":[]}
            for _group_idx, _group_name in enumerate(config["Labels"].keys()):

                if _group_name in config["Signal_samples"]:
                    isSignal = "signal" 
                elif _group_name in config["Data"]:
                    isSignal = "data" 
                else:
                    isSignal = "background"
                
                # Accumulate the dataset for the particular data group as specified in config “Labels”.
                for _idx, _histogram_data in enumerate(config["Labels"][_group_name]):
                    
                    # Rescaling according to cross-section and luminosity
                    # print("Histogram: ",  _histname)
                    # print("Reading data:", _histogram_data + "_" + _categ)
                    # hist = file.Get(_histogram_data + "_" + _categ)
                    # print(_histogram_data)
                    hist = file.Get(_histogram_data)
                    
                    if isSignal != "data": # Scaling of the MC to the lumi and xsection
                        # N = cutflow[_histogram_data]["all"]["NanDrop"] #After Nan dropper
                        N = cutflow[_histogram_data]["all"]["Before cuts"]
                        hist.Scale( (xsec[_histogram_data] * config["luminosity"]) / N)

                    if _histname in config["SetupBins"]:
                        rebin_setup = config["SetupBins"][_histname][2]
                        if type(rebin_setup) == list:
                            hist = hist.Rebin(len(rebin_setup)-1, hist.GetName()+"_rebin", np.array(rebin_setup))
                        else: 
                            hist.Rebin(rebin_setup)
                            
                        if config["SetupBins"][_histname][4]:
                            for bin_i, label in enumerate(config["SetupBins"][_histname][4]):
                                hist.GetXaxis().SetBinLabel(bin_i, label)
                                hist.GetXaxis().SetTitle("")
          
                    if _idx == 0:
                        _histograms[isSignal].append(hist)
                    else:
                        _histograms[isSignal][-1].Add(hist)

                if isSignal == "signal":
                    color_setup = config["Signal_samples"][_group_name]  
                    line_color = color_setup[1]
                    fill_color = color_setup[0]
                    _histograms[isSignal][-1].SetLineStyle(2)
                    _histograms[isSignal][-1].SetMarkerSize(0)
                    _histograms[isSignal][-1].SetLineWidth(4)
                elif isSignal == "data":
                    color_setup = config["Data"][_group_name]  
                    line_color = color_setup[1]
                    fill_color = color_setup[0] 
                    _histograms[isSignal][-1].SetMarkerStyle(8)
                    _histograms[isSignal][-1].SetMarkerSize(2)
                    _histograms[isSignal][-1].SetLineWidth(1)
                else:
                    color_setup = config["MC_bkgd"][_group_name]  
                    line_color = color_setup[1]
                    fill_color = color_setup[0]
                    _histograms[isSignal][-1].SetMarkerSize(0)
                    _histograms[isSignal][-1].SetLineWidth(4)
                
                _histograms[isSignal][-1].SetLineColor(line_color)
                _histograms[isSignal][-1].SetFillColor(fill_color)
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

            # read the binning if available:
            if _histname in config["SetupBins"]:
                xrange_min = config["SetupBins"][_histname][0]
                xrange_max = config["SetupBins"][_histname][1]
                overflow =  bool(config["SetupBins"][_histname][3])
            else:
                xrange_min = _histograms["background"][0].GetXaxis().GetXmin()
                xrange_max = _histograms["background"][0].GetXaxis().GetXmax()
                overflow =  True

            if isData:
                root_plot1D(
                    l_hist = _histograms_background_sorted,
                    l_hist_overlay = _histograms["data"],
                    outfile = output + "/" + os.path.splitext(os.path.basename(_histfile))[0] + ".png",
                    xrange = [xrange_min, xrange_max],
                    yrange = (0.001,  10000*y_max),
                    logx = False, logy = True,
                    logx_ratio = False, logy_ratio = False,
                    include_overflow = overflow,
                    xtitle = _histograms["background"][0].GetXaxis().GetTitle(),
                    ytitle = "events",
                    xtitle_ratio = _histograms["background"][0].GetXaxis().GetTitle(),
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
            
            else:
                root_plot1D(
                    l_hist = _histograms_background_sorted,
                    l_hist_overlay = _histograms["signal"],
                    outfile = output + "/" + os.path.splitext(os.path.basename(_histfile))[0] + ".png",
                    xrange = [xrange_min, xrange_max],
                    # yrange = (0.0001,  1000*y_max),
                    yrange = (0.0,  1.5*y_max),
                    logx = False, logy = False,
                    include_overflow = overflow,
                    xtitle = _histograms["background"][0].GetXaxis().GetTitle(),
                    ytitle = "events",
                    xtitle_ratio = _histograms["background"][0].GetXaxis().GetTitle(),
                    ytitle_ratio = "S/#sqrt{S+B}",
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
                    ratio_mode = "SB",
                    yrange_ratio = (1E-04, 1),
                    draw_errors = True
                )

def plot2D(histfiles, histnames, config, xsec, cutflow, output_path):

    categories_list = list(itertools.product(*config["Categories"]))
    categories_list = [f"{cat1}_{cat2}" for cat1,cat2 in categories_list]

    for _ci, _categ in enumerate(categories_list):

        output = output_path + "/" + _categ

        for _histfile, _histname in zip(histfiles,histnames):

            # The plots are done for specific triggers specified in config.
            if not any([cut in str(_histfile) for cut in config["cuts"]]):
                continue

            file = ROOT.TFile.Open(str(_histfile), 'read')

            _histograms = {"background":[], "signal":[]}
            for _group_idx, _group_name in enumerate(config["Labels"].keys()):

                isSignal = "signal" if _group_name in config["Signal_samples"] else "background"
                
                # Accumulate the dataset for the particular data group as specified in config “Labels”.
                for _idx, _histogram_data in enumerate(config["Labels"][_group_name]):
                    
                    # Rescaling according to cross-section and luminosity
                    # print("Reading hist:", _histogram_data + "_" + _categ)
                    hist = file.Get(_histogram_data + "_" + _categ)
                    N = cutflow[_histogram_data]["all"]["Before cuts"] #After Nan dropper
                    hist.Scale( (xsec[_histogram_data] * config["luminosity"]) / N)

                    if _histname in config["SetupBins"]:
                        hist.Rebin2D(config["SetupBins"][_histname][4], config["SetupBins"][_histname][5])
                    
                    if _idx == 0:
                        _histograms[isSignal].append(hist)
                    else:
                        _histograms[isSignal][-1].Add(hist)

                if isSignal == "signal":
                    color_setup = config["Signal_samples"][_group_name]  
                    line_color = color_setup[1]
                    fill_color = color_setup[0]
                    _histograms[isSignal][-1].SetLineStyle(2)
                else:
                    color_setup = config["MC_bkgd"][_group_name]  
                    line_color = color_setup[1]
                    fill_color = color_setup[0]
                
                # _histograms[isSignal][-1].SetLineColor(line_color)
                # _histograms[isSignal][-1].SetFillColor(fill_color)
                # _histograms[isSignal][-1].SetLineWidth(2)
                # _histograms[isSignal][-1].SetMarkerSize(0)
                # _histograms[isSignal][-1].SetTitle(_group_name)
            
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

            # read the binning if available:
            if _histname in config["SetupBins"]:
                xrange_min = config["SetupBins"][_histname][0]
                xrange_max = config["SetupBins"][_histname][1]
                yrange_min = config["SetupBins"][_histname][2]
                yrange_max = config["SetupBins"][_histname][3]
                text_colz = config["SetupBins"][_histname][6]
                log_z = config["SetupBins"][_histname][7]
            else:
                xrange_min = _histograms["background"][0].GetXaxis().GetXmin()
                xrange_max = _histograms["background"][0].GetXaxis().GetXmax()
                yrange_min = _histograms["background"][0].GetYaxis().GetXmin()
                yrange_max = _histograms["background"][0].GetYaxis().GetXmax()
                text_colz = False
                log_z = True

            root_plot2D(
                l_hist = _histograms_background_sorted,
                l_hist_overlay = _histograms["signal"],
                outfile = output + "/" + os.path.splitext(os.path.basename(_histfile))[0] + ".png",
                xrange = [xrange_min, xrange_max],
                yrange = [yrange_min, yrange_max],
                # yrange = (1,  1.1*y_max),
                logx = False, logy = False, logz = log_z,
                # ytitle = _histograms["signal"][0].GetYaxis().GetTitle(),
                xtitle = _histograms["background"][0].GetXaxis().GetTitle(),
                ytitle = _histograms["background"][0].GetYaxis().GetTitle(),
                centertitlex = True, centertitley = True,
                centerlabelx = False, centerlabely = False,
                gridx = True, gridy = True,
                ndivisionsx = None,
                stackdrawopt = "",
                ratio_mode="SB",
                # normilize = True,
                normilize_overlay = False,
                legendpos = "UR",
                legendtitle = f"",
                legendncol = 3,
                legendtextsize = 0.03,
                legendwidthscale = 1.9,
                legendheightscale = 0.4,
                lumiText = "2018 (13 TeV)",
                signal_to_background_ratio = True,
                text_colz = text_colz,
            )