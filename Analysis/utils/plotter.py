import ROOT
ROOT.gROOT.SetBatch(True)
import itertools
import numpy as np
import os

from .utils import ColorIterator, root_plot1D, root_plot2D

def plot_predict(dirname, config, xsec, cutflow, output_path):
    
    for prediction_bin, data_bin in zip(config["prediction_hist"]["predictions"], config["prediction_hist"]["bin_data"]):
        for hist in config["prediction_hist"]["hists"]:

            print("Prediction bin:", prediction_bin, "Histogram:", hist)
            
            # Extracting histogram information
            x_min, x_max, rebin_setup, overflow, bin_labels = \
                config["prediction_hist"]["hists"][hist]
            cut = config["prediction_hist"]["cut"]
            
            path_predict = dirname+"/"+cut+"_"+hist+"_yield_"+prediction_bin+".root"
            print(path_predict)
            file_predict = ROOT.TFile.Open(path_predict, 'read')
            hist_prediction = None
            for data_group in config["Data"].keys():
                for data_name in config["Labels"][data_group]:
                    print("Extract prediction:", data_name)
                    if hist_prediction == None:
                        hist_prediction = file_predict.Get(data_name)
                    else:
                        hist_prediction.Add(file_predict.Get(data_name))
                        
            if type(rebin_setup) == list:
                hist_prediction = hist_prediction.Rebin(len(rebin_setup)-1, hist_prediction.GetName()+"_rebin", np.array(rebin_setup, dtype=np.double))
            else: 
                hist_prediction.Rebin(rebin_setup)
                
            hist_prediction.SetMarkerStyle(21)
            hist_prediction.SetMarkerColor(31)
            hist_prediction.SetLineColor(30)
            hist_prediction.SetFillColor(0)
            hist_prediction.SetLineWidth(3)
            hist_prediction.SetTitle(f"Pred. ({str(prediction_bin)})")
            # print("Prediction:")
            # print(hist_prediction.Print("all"))
            hists_main = [hist_prediction]
            
            # ---- Part to assign signal histogram
            path_data = dirname+"/"+cut+"_"+hist+"_pass.root"
            print(path_data)
            file_n_pass_sig = ROOT.TFile.Open(path_data, 'read') 
            signal_hists = []
            if config["prediction_hist"]["plot_signal"]:
                for _group_idx, _group_name in enumerate(config["Signal_samples"]):
                    # Accumulate the dataset for the particular data group as specified in config “Labels”.
                    for _dataset_idx, _histogram_data in enumerate(config["Labels"][_group_name]):
                        print("Adding signal dataset:", _histogram_data)
                        if isinstance(data_bin, str):
                            _hist = file_n_pass_sig.Get(_histogram_data+"_"+data_bin)
                        elif isinstance(data_bin, int):
                            _hist = file_n_pass_sig.Get(_histogram_data)
                            _hist = _hist.ProjectionX(_histogram_data+"_proj", data_bin, data_bin)
                        _hist.SetDirectory(0)
                        N = cutflow[_histogram_data]["all"]["Before cuts"]
                        scale =  xsec[_histogram_data] * config["luminosity"] / N
                        # _hist.Scale( (xsec[_histogram_data] * config["luminosity"]) / N)
                        for bin_i in range(0, _hist.GetNbinsX()+2):
                            _hist.SetBinContent(bin_i, _hist.GetBinContent(bin_i)*scale)
                            _hist.SetBinError(bin_i, _hist.GetBinError(bin_i)*scale)
                        if _dataset_idx == 0:
                            signal_hists.append(_hist)
                        else:
                            signal_hists[-1].Add(_hist)
                    if type(rebin_setup) == list:
                        signal_hists[-1] = signal_hists[-1].Rebin(len(rebin_setup)-1, signal_hists[-1].GetName()+"_rebin", np.array(rebin_setup, dtype=np.double))
                    else: 
                        signal_hists[-1].Rebin(rebin_setup)
                    color_setup = config["Signal_samples"][_group_name]  
                    line_color = color_setup[1]
                    fill_color = color_setup[0]
                    signal_hists[-1].SetLineStyle(2)
                    signal_hists[-1].SetMarkerSize(0)
                    signal_hists[-1].SetLineWidth(2)
                    signal_hists[-1].SetLineColor(line_color)
                    signal_hists[-1].SetFillColor(fill_color)
                    signal_hists[-1].SetTitle(_group_name)

            # print(signal_hists)
            
            # ---- Part to assign true histogram
            file_n_pass = ROOT.TFile.Open(path_data, 'read')
            if config["prediction_hist"]["plot_unblind"]:
                hist_data = None
                for data_group in config["Data"].keys():
                    for data_name in config["Labels"][data_group]:
                        # print(file_n_pass.ls())
                        # print(data_name+"_"+data_bin)
                        if isinstance(data_bin, str):
                            _hist_data = file_n_pass.Get(data_name+"_"+data_bin)
                        elif isinstance(data_bin, int):
                            _hist_data = file_n_pass.Get(data_name)
                            _hist_data = _hist_data.ProjectionX(data_name+"_proj", data_bin, data_bin)
                       
                        if hist_data == None:
                            hist_data = _hist_data
                        else: hist_data.Add(_hist_data)
                        
                if type(rebin_setup) == list:
                    hist_data = hist_data.Rebin(len(rebin_setup)-1, hist_data.GetName()+"_rebin", np.array(rebin_setup, dtype=np.double))
                else: 
                    hist_data.Rebin(rebin_setup)
                hist_data.SetDirectory(0)
                hist_data.SetMarkerStyle(8)
                hist_data.SetMarkerSize(1)
                hist_data.SetMarkerColor(1)
                hist_data.SetLineWidth(3)
                hist_data.SetTitle(f"Data ({(str(data_bin))})")
                signal_hists.append(hist_data)
                # print("Data:")
                # print(hist_data.Print("all"))            

            root_plot1D(
                l_hist = hists_main,
                l_hist_overlay = signal_hists,
                # l_hist_overlay = [hist_data] if config["prediction_hist"]["plot_unblind"] else [],
                outfile = output_path + "/" + hist + "_" + prediction_bin + ".png",
                xrange = [x_min, x_max],
                yrange = (0.01,  1000*hist_prediction.GetMaximum()),
                logx = False, logy = True,
                logx_ratio = False, logy_ratio = False,
                include_overflow = overflow,
                xtitle = signal_hists[-1].GetXaxis().GetTitle(),
                ytitle = "events",
                xtitle_ratio = signal_hists[-1].GetXaxis().GetTitle(),
                ytitle_ratio = "DATA/MC",
                centertitlex = True, centertitley = True,
                centerlabelx = False, centerlabely = False,
                gridx = True, gridy = True,
                ndivisionsx = None,
                stackdrawopt = "",
                ratio_mode = "DATA",
                normilize = False,
                normilize_overlay = False,
                legendpos = "UL",
                legendtitle = f"",
                legendncol = 2,
                legendtextsize = 0.040,
                legendwidthscale = 2.0,
                legendheightscale = 3.0,
                lumiText = "2018 (13 TeV)",
                yrange_ratio = (0.0, 3.0),
                ndivisionsy_ratio = (6, 5, 0),
                signal_to_background_ratio = True,
                draw_errors = False
            )


def plot1D(histfiles, histnames, config, xsec, cutflow, output_path, isData):

    categories_list = list(itertools.product(*config["Categories"]))
    # categories_list = [f"{cat1}_{cat2}_{cat3}" for cat1,cat2,cat3 in categories_list]
    categories_list = ["_".join(cat) for cat in categories_list]
    categories_list = [""]

    for _ci, _categ in enumerate(categories_list):

        output = output_path + "/" + _categ

        for _histfile, _histname in zip(histfiles,histnames):
            
            # print([cut in str(_histfile) for cut in config["cuts"]])
            if not any([cut in str(_histfile) for cut in config["cuts"]]):
                continue

            # print("OPEN:", _histfile)
            file = ROOT.TFile.Open(str(_histfile), 'read')

            _histograms = {"background":[], "signal":[], "data":[]}
            for _group_idx, _group_name in enumerate(config["Labels"].keys()):

                if _group_name in config["Signal_samples"]:
                    isSignal = "signal" 
                elif _group_name in config["Data"]:
                    isSignal = "data" 
                else:
                    isSignal = "background"
                
                has_group_entry = False
                # Accumulate the dataset for the particular data group as specified in config “Labels”.
                for _idx, _histogram_data in enumerate(config["Labels"][_group_name]):
                    
                    # Rescaling according to cross-section and luminosity
                    # print("Reading data:", _histogram_data + "_" + _categ)
                    # hist = file.Get(_histogram_data + "_" + _categ)
                    hist = file.Get(_histogram_data)

                    if not hist:
                        print("Warning: Histogram not found! ", end='')
                        print("Histogram->", file, _histname, _histogram_data)
                        continue
                    
                    if isSignal != "data": # Scaling of the MC to the lumi and xsection
                        if config["DY_stitching_applied"] and (
                                "DYJetsToLL_M-50" in _histogram_data or
                                "DY1JetsToLL_M-50" in _histogram_data or
                                "DY2JetsToLL_M-50" in _histogram_data or
                                "DY3JetsToLL_M-50" in _histogram_data or
                                "DY4JetsToLL_M-50" in _histogram_data ):
                            # print("Stitching:", _histogram_data)
                            hist.Scale(config["luminosity"])
                        elif config["W_stitching_applied"] and (
                                ("WJetsToLNu" in _histogram_data and (not "TTWJets" in _histogram_data)) or
                                "W1JetsToLNu" in _histogram_data or
                                "W2JetsToLNu" in _histogram_data or
                                "W3JetsToLNu" in _histogram_data or
                                "W4JetsToLNu" in _histogram_data ):
                            # print("Stitching:", _histogram_data)
                            hist.Scale(config["luminosity"])
                        else:
                            # N = cutflow[_histogram_data]["all"]["NanDrop"] #After Nan dropper
                            N = cutflow[_histogram_data]["all"]["Before cuts"]
                            hist.Scale( (xsec[_histogram_data] * config["luminosity"]) / N)

                    if _histname in config["SetupBins"]:
                        rebin_setup = config["SetupBins"][_histname][2]
                        if type(rebin_setup) == list:
                            hist = hist.Rebin(len(rebin_setup)-1, hist.GetName()+"_rebin", np.array(rebin_setup, dtype=np.double))
                        else: 
                            hist.Rebin(rebin_setup)
                            
                        if config["SetupBins"][_histname][4]:
                            for bin_i, label in enumerate(config["SetupBins"][_histname][4]):
                                hist.GetXaxis().SetBinLabel(bin_i, label)
                                hist.GetXaxis().SetTitle("")
          
                    if not has_group_entry:
                        # print("Append:", _histogram_data)
                        _histograms[isSignal].append(hist)
                        has_group_entry = True
                    else:
                        # print("Add:", _histogram_data)
                        _histograms[isSignal][-1].Add(hist)

                if has_group_entry: # if there is at least one histogram in input
                    
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
                        _histograms[isSignal][-1].SetMarkerSize(1)
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
                    print("Set title:", _group_name)
            # exit()
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
                    # yrange = (0.0,  1.5*y_max), 
                    # logx = False, logy = False,
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
                    legendtextsize = 0.035,
                    legendwidthscale = 1.9,
                    legendheightscale = 0.36,
                    lumiText = "2018 (13 TeV)",
                    signal_to_background_ratio = True,
                    ratio_mode = "DATA",
                    yrange_ratio = (0.0, 2.0),
                    draw_errors = True
                )
            
            else:
                root_plot1D(
                    l_hist = _histograms_background_sorted,
                    l_hist_overlay = _histograms["signal"],
                    outfile = output + "/" + os.path.splitext(os.path.basename(_histfile))[0] + ".png",
                    xrange = [xrange_min, xrange_max],
                    yrange = (0.0001,  1000*y_max),
                    # yrange = (0.0,  1.5*y_max),
                    logx = False, logy = True,
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
                    signal_to_background_ratio = False,
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
