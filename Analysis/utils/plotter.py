import ROOT
ROOT.gROOT.SetBatch(True)
import itertools
import numpy as np
import os

from .utils import ColorIterator, root_plot1D, root_plot2D, root_plots2D_simple
ROOT.gInterpreter.Declare('#include "utils/histogram2d.cpp"')

def plot_predict(dirname, config, xsec, cutflow, output_path):
    
    for prediction_bin, data_bin in zip(config["prediction_hist"]["predictions"], config["prediction_hist"]["bin_data"]):
        for hist in config["prediction_hist"]["hists"]:

            print("Prediction bin:", prediction_bin, "Histogram:", hist)
            
            # Extracting histogram information
            x_min, x_max, rebin_setup, overflow, bin_labels = \
                config["prediction_hist"]["hists"][hist]
            cut = config["prediction_hist"]["cut"]
            
            # ---- Part to assign prediction histogram
            path_predict = dirname+"/"+cut+"_"+hist+"_yield_"+prediction_bin+".root"
            print(path_predict)
            file_predict = ROOT.TFile.Open(path_predict, 'read')
            # print(file_predict.ls())
            hist_prediction = None
            isDATA = False
            # for data_group in config["Data"].keys():
            #     for data_name in config["Labels"][data_group]:

            for _group_idx, _group_name in enumerate(config["Labels"].keys()):
                # if (not _group_name in config["MC_bkgd"]) and (not _group_name in config["Signal_samples"]):
                #     continue
                # if not( _group_name in config["Data"]): continue
                # Accumulate the dataset for the particular data group as specified in config "Labels".
                for _idx, data_name in enumerate(config["Labels"][_group_name]):

                    print("Extract prediction:", data_name)
                    open_tag = data_name+"/hist"
                    _hist_predict = file_predict.Get(open_tag)
                    if not _hist_predict:
                        print("Warning: Histogram not found! ", end='')
                        print("Histogram->", file_predict, open_tag)
                        continue
                    _hist_predict.SetDirectory(0)

                    if _group_name in config["Data"]:
                        isDATA = True
                        pass
                    else:
                        if isDATA: raise("Can not combine data and MC")

                        if config["DY_stitching_applied"] and (
                                "DYJetsToLL_M-50" in data_name or
                                "DY1JetsToLL_M-50" in data_name or
                                "DY2JetsToLL_M-50" in data_name or
                                "DY3JetsToLL_M-50" in data_name or
                                "DY4JetsToLL_M-50" in data_name ):
                            # print("Stitching:", data_name)
                            _hist_predict.Scale(config["luminosity"])
                        elif config["W_stitching_applied"] and (
                                ("WJetsToLNu" in data_name and (not "TTWJets" in data_name)) or
                                "W1JetsToLNu" in data_name or
                                "W2JetsToLNu" in data_name or
                                "W3JetsToLNu" in data_name or
                                "W4JetsToLNu" in data_name ):
                            # print("Stitching:", data_name)
                            _hist_predict.Scale(config["luminosity"])
                        else:
                            # N = cutflow[_histogram_data]["all"]["NanDrop"] #After Nan dropper
                            N = cutflow[data_name]["all"]["BeforeCuts"]
                            _hist_predict.Scale( (xsec[data_name] * config["luminosity"]) / N)

                    if hist_prediction == None:
                        hist_prediction = _hist_predict
                    else:
                        hist_prediction.Add(_hist_predict)
                        
            if type(rebin_setup) == list:
                hist_prediction = hist_prediction.Rebin(len(rebin_setup)-1, hist_prediction.GetName()+"_rebin", np.array(rebin_setup, dtype=np.double))
            else: 
                hist_prediction.Rebin(rebin_setup)
                
            hist_prediction.SetMarkerStyle(21)
            hist_prediction.SetMarkerColor(618)
            hist_prediction.SetLineColor(618)
            hist_prediction.SetFillColor(0)
            hist_prediction.SetLineWidth(3)
            hist_prediction.SetTitle(f"Pred. ({str(prediction_bin)})")
            
            # print("Prediction:")
            # print(hist_prediction.Print("all"))
            hists_main = [hist_prediction]
            
            signal_hists = []
            # ---- Part to assign true histogram
            path_data = dirname+"/"+cut+"_"+hist+"_pass.root"
            print(path_data)
            file_n_pass_sig = ROOT.TFile.Open(path_data, 'read')
            file_n_pass = ROOT.TFile.Open(path_data, 'read')
            # print(file_n_pass_sig.ls())
            # print(file_n_pass.ls())
            if config["prediction_hist"]["plot_unblind"]:
                hist_data = None
                isDATA = False
                # for data_group in config["Data"].keys():
                #     for data_name in config["Labels"][data_group]:

                for _group_idx, _group_name in enumerate(config["Labels"].keys()):
                    # if (not _group_name in config["MC_bkgd"]) and (not _group_name in config["Signal_samples"]):
                    #     continue
                    # if not( _group_name in config["Data"]): continue
                    # Accumulate the dataset for the particular data group as specified in config "Labels".
                    for _idx, data_name in enumerate(config["Labels"][_group_name]):

                        # print(file_n_pass.ls())
                        # print(data_name+"_"+data_bin)
                        if isinstance(data_bin, str):
                            open_tag = data_name+"/"+data_bin+"/hist"
                            _hist_data = file_n_pass.Get(open_tag)
                        elif isinstance(data_bin, int):
                            _hist_data = file_n_pass.Get(data_name)
                            _hist_data = _hist_data.ProjectionX(data_name+"_proj", data_bin, data_bin)
                        if not _hist_data:
                            print("Warning: Histogram not found! ", end='')
                            print("Histogram->", file_n_pass, open_tag)
                            continue
                        _hist_data.SetDirectory(0)

                        if _group_name in config["Data"]:
                            isDATA = True
                            pass
                        else:
                            if isDATA: raise("Can not combine data and MC")

                            if config["DY_stitching_applied"] and (
                                    "DYJetsToLL_M-50" in data_name or
                                    "DY1JetsToLL_M-50" in data_name or
                                    "DY2JetsToLL_M-50" in data_name or
                                    "DY3JetsToLL_M-50" in data_name or
                                    "DY4JetsToLL_M-50" in data_name ):
                                # print("Stitching:", data_name)
                                _hist_data.Scale(config["luminosity"])
                            elif config["W_stitching_applied"] and (
                                    ("WJetsToLNu" in data_name and (not "TTWJets" in data_name)) or
                                    "W1JetsToLNu" in data_name or
                                    "W2JetsToLNu" in data_name or
                                    "W3JetsToLNu" in data_name or
                                    "W4JetsToLNu" in data_name ):
                                # print("Stitching:", data_name)
                                _hist_data.Scale(config["luminosity"])
                            else:
                                # N = cutflow[_histogram_data]["all"]["NanDrop"] #After Nan dropper
                                N = cutflow[data_name]["all"]["BeforeCuts"]
                                _hist_data.Scale( (xsec[data_name] * config["luminosity"]) / N)

                       
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

            # ---- Part to assign signal histogram
            if config["prediction_hist"]["plot_signal"]:
                for _group_idx, _group_name in enumerate(config["Signal_samples"]):
                    # Accumulate the dataset for the particular data group as specified in config "Labels".
                    for _dataset_idx, _histogram_data in enumerate(config["Labels"][_group_name]):
                        print("Adding signal dataset:", _histogram_data)
                        if isinstance(data_bin, str):
                            open_tag = _histogram_data+"/"+data_bin+"/hist"
                            _hist = file_n_pass_sig.Get(open_tag)
                        elif isinstance(data_bin, int):
                            _hist = file_n_pass_sig.Get(_histogram_data)
                            _hist = _hist.ProjectionX(_histogram_data+"_proj", data_bin, data_bin)
                        _hist.SetDirectory(0)
                        N = cutflow[_histogram_data]["all"]["BeforeCuts"]
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
                    line_style = color_setup[2]
                    signal_hists[-1].SetLineStyle(line_style)
                    signal_hists[-1].SetMarkerSize(0)
                    signal_hists[-1].SetLineWidth(2)
                    signal_hists[-1].SetLineColor(line_color)
                    signal_hists[-1].SetFillColor(fill_color)
                    signal_hists[-1].SetTitle(_group_name)

            # print(signal_hists)

            # remove $ char from the title name
            x_axis_title = hist_prediction.GetXaxis().GetTitle()
            x_axis_title = x_axis_title.replace("$", "")

            root_plot1D(
                l_hist = hists_main,
                l_hist_overlay = signal_hists,
                # l_hist_overlay = [hist_data] if config["prediction_hist"]["plot_unblind"] else [],
                outfile = output_path + "/" + hist + "_" + prediction_bin + ".pdf",
                xrange = [x_min, x_max],
                yrange = (0.01,  1000*hist_prediction.GetMaximum()),
                logx = False, logy = True,
                logx_ratio = False, logy_ratio = False,
                include_overflow = overflow,
                xtitle = x_axis_title,
                ytitle = "Events",
                xtitle_ratio = x_axis_title,
                ytitle_ratio = "Data/Pred.",
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
                legendheightscale = 4.0,
                lumiText = "2018 (13 TeV)",
                yrange_ratio = (0.0, 2.0),
                ndivisionsy_ratio = (5, 5, 0),
                signal_to_background_ratio = True,
                draw_errors = True
            )

def plot_predict2D(dirname, config, xsec, cutflow, output_path):
    
    for prediction_bin, data_bin in zip(config["prediction_hist2D"]["predictions"], config["prediction_hist2D"]["bin_data"]):
        for hist in config["prediction_hist2D"]["hists"]:

            # ~~~~~~~~~~~~~~~ Prediction 
            print("Prediction bin:", prediction_bin, "Histogram:", hist)
            # Extracting histogram information
            axis_rebin, overflow, bin_labels = \
                config["prediction_hist2D"]["hists"][hist]
            cut = config["prediction_hist2D"]["cut"]
            path_predict = dirname+"/"+cut+"_"+hist+"_"+prediction_bin+".root"
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
            
            x_axis = axis_rebin["x_axis"]
            y_axis =  np.array(axis_rebin["y_axis"], dtype=np.double)
            xmin, xmax = float(x_axis[0][0]), float(x_axis[0][-1])
    
            hist_prediction_nonunif = ROOT.Histogram_2D("nonunif_pred", y_axis, xmin, xmax)
            for i in range(len(x_axis)):
                hist_prediction_nonunif.add_x_binning_by_index(i, np.array(x_axis[i], dtype=np.double))
            
            try:
                hist_prediction_nonunif.th2d_add(hist_prediction)
            except:
                print("Error is inside project nonunif histograms")
                print("Error: can not add histogram because of the inconsistency")
                del hist_prediction_nonunif
                exit()
            hist_prediction_nonunif.print("")
            hist_prediction_nonunif_th2 = hist_prediction_nonunif.get_weights_th2d_simpl("hist_prediction_nonunif_th2",
                                                                                         "hist_prediction_nonunif_th2")
            del hist_prediction_nonunif
            
            ## ~~~~~~~~~~~~~~~ Signal
            path_data = dirname+"/"+cut+"_"+hist+"_pass.root"
            print(path_data)
            file_n_pass_sig = ROOT.TFile.Open(path_data, 'read')
            _signal_name = config["prediction_hist2D"]["signal_model"]
            signal_hists = None
            for _dataset_idx, _histogram_data in enumerate(config["Labels"][_signal_name]):
                print("Adding signal dataset:", _histogram_data)
                if isinstance(data_bin, str):
                    _hist = file_n_pass_sig.Get(_histogram_data+"_"+data_bin)
                elif isinstance(data_bin, int):
                    _hist = file_n_pass_sig.Get(_histogram_data)
                    _hist = _hist.ProjectionX(_histogram_data+"_proj", data_bin, data_bin)
                _hist.SetDirectory(0)
                N = cutflow[_histogram_data]["all"]["Before cuts"]
                scale =  xsec[_histogram_data] * config["luminosity"] / N
                _hist.Scale(scale)
                if signal_hists == None:
                    signal_hists = _hist
                else:
                    signal_hists.Add(_hist)
                    
            hist_sig_nonunif = ROOT.Histogram_2D("nonunif_sig", y_axis, xmin, xmax)
            for i in range(len(x_axis)):
                hist_sig_nonunif.add_x_binning_by_index(i, np.array(x_axis[i], dtype=np.double))
            
            try:
                hist_sig_nonunif.th2d_add(signal_hists)
            except:
                print("Error is inside project nonunif histograms")
                print("Error: can not add histogram because of the inconsistency")
                del hist_sig_nonunif
                exit()
            hist_sig_nonunif.print("")
            hist_sig_nonunif_th2 = hist_sig_nonunif.get_weights_th2d_simpl("hist_sig_nonunif_th2",
                                                                           "hist_sig_nonunif_th2")
            del hist_sig_nonunif
            
        root_plots2D_simple(
                hist_prediction_nonunif_th2,
                outfile = output_path + "/" + hist + "_" + prediction_bin +"_predict"+ ".pdf",
                yrange=(y_axis[0], y_axis[-1]),
                xrange=(x_axis[0][0], x_axis[0][-1]),
                logx = False, logy = False, logz = False,
                title = "",
                # xtitle = "p^{miss}_{T} [GeV]", ytitle = "m_{T2} [GeV]",
                # xtitle = "p^{miss}_{T} [GeV]", ytitle = "#Sigma m_{T} [GeV]",
                xtitle = "m_{T2} [GeV]", ytitle = "#Sigma m_{T} [GeV]",
                centertitlex = True, centertitley = True,
                centerlabelx = False, centerlabely = False,
                gridx = False, gridy = False,
                CMSextraText = "Private work (Prediction)",
                lumiText = "(13 TeV)",
            )
        
        root_plots2D_simple(
                hist_sig_nonunif_th2,
                outfile = output_path + "/" + hist + "_" + prediction_bin +"_signal"+ ".pdf",
                yrange=(y_axis[0], y_axis[-1]),
                xrange=(x_axis[0][0], x_axis[0][-1]),
                logx = False, logy = False, logz = False,
                title = "",
                # xtitle = "p^{miss}_{T} [GeV]", ytitle = "m_{T2} [GeV]",
                # xtitle = "p^{miss}_{T} [GeV]", ytitle = "#Sigma m_{T} [GeV]",
                xtitle = "m_{T2} [GeV]", ytitle = "#Sigma m_{T} [GeV]",
                centertitlex = True, centertitley = True,
                centerlabelx = False, centerlabely = False,
                gridx = False, gridy = False,
                CMSextraText = f"Private work (Signal 250 GeV,10 cm)",
                lumiText = "(13 TeV)",
            )
        
        hist_prediction_nonunif_th2.Add(hist_sig_nonunif_th2)
        hist_sig_nonunif_th2.Divide(hist_prediction_nonunif_th2)
        root_plots2D_simple(
                hist_sig_nonunif_th2,
                outfile = output_path + "/" + hist + "_" + prediction_bin +"_RATIO"+ ".pdf",
                yrange=(y_axis[0], y_axis[-1]),
                xrange=(x_axis[0][0], x_axis[0][-1]),
                logx = False, logy = False, logz = False,
                title = "",
                # xtitle = "p^{miss}_{T} [GeV]", ytitle = "m_{T2} [GeV]",
                # xtitle = "p^{miss}_{T} [GeV]", ytitle = "#Sigma m_{T} [GeV]",
                xtitle = "m_{T2} [GeV]", ytitle = "#Sigma m_{T} [GeV]",
                centertitlex = True, centertitley = True,
                centerlabelx = False, centerlabely = False,
                gridx = False, gridy = False,
                CMSextraText = f"Private work ( S/(S+B) )",
                lumiText = "(13 TeV)",
            )


def plot1D(histfiles, histnames, config, xsec, cutflow, output_path, isData):

    categories_list = list(itertools.product(*config["Categories"]))
    # categories_list = [f"{cat1}_{cat2}_{cat3}" for cat1,cat2,cat3 in categories_list]
    categories_list = ["_".join(cat) for cat in categories_list]
    # categories_list = [""]
    
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
                    # hist = file.Get(_histogram_data)
                    if not _categ=="":
                        hist = file.Get(_histogram_data + "_" + _categ)
                    else:
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
                        line_style = color_setup[2]
                        _histograms[isSignal][-1].SetLineStyle(line_style)
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
                    yrange = (0.001,  1000*y_max),
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
                    yrange = (0.001,  1000*y_max),
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
                    signal_to_background_ratio = True,
                    ratio_mode = "SB",
                    yrange_ratio = (0.0, 2.0),
                    draw_errors = True
                )

def plotBrMC(hist_path, config, xsec, cutflow, output_path, is_per_flavour=False):
    
    '''
    This code is used to plot the MC branching ratio reading 2D histogram
    and projecting it to the x-axis for bin 1 bin2 and bin 3 which coorespond to
    0-tight, 1--tight and 2-tight region
    '''
    print(hist_path)
    
    _histograms = {}
    if is_per_flavour:
        flavours = config["mixing_hists"]["flavours"]
        file = ROOT.TFile.Open(hist_path+"_0.root", 'read') 
        file2 = ROOT.TFile.Open(hist_path+"_1.root", 'read')
        # print(file.ls())
        # print(file2.ls())
    else:
        file = ROOT.TFile.Open(hist_path+".root", 'read') 
        flavours = [None]

    for flav in flavours:
        for _group_idx, _group_name in enumerate(config["Labels"].keys()):

            if _group_name in config["Signal_samples"]:
                raise ValueError("Signal samples are not allowed in this plotter")
            elif _group_name in config["Data"]:
                raise ValueError("Data samples are not allowed in this plotter")

            # has_group_entry = False
            # Accumulate the dataset for the particular data group as specified in config “Labels”.
            for _idx, _histogram_data in enumerate(config["Labels"][_group_name]):
                
                if is_per_flavour:
                    print("Reading hist:", f"{_histogram_data}/{flav}/hist")
                    hist = file.Get(f"{_histogram_data}/{flav}/hist")
                    hist2 = file2.Get(f"{_histogram_data}/{flav}/hist")
                else:
                    print("Reading hist:", _histogram_data)
                    hist = file.Get(_histogram_data+"/hist")
                
                if not hist:
                    print("Warning: Histogram not found! ", end='')
                    print("Histogram->", file, _histogram_data)
                    continue
                
                if is_per_flavour:
                    hist.Add(hist2)
                # hist = hist2

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
                    N = cutflow[_histogram_data]["all"]["BeforeCuts"]
                    hist.Scale( (xsec[_histogram_data] * config["luminosity"]) / N)

                if is_per_flavour:
                    _group_name = flav

                if not _group_name in _histograms.keys():
                    print("Append:", _histogram_data)
                    _histograms[_group_name] = hist
                else:
                    print("Add:", _histogram_data)
                    _histograms[_group_name].Add(hist)
                    
                if _group_name in _histograms.keys(): # if there is at least one histogram in input
                
                    color_setup = config["MC_bkgd"][_group_name]  
                    line_color = color_setup[1]
                    fill_color = color_setup[0]
                    _histograms[_group_name].SetMarkerSize(0)
                    _histograms[_group_name].SetLineWidth(4)
                    
                    _histograms[_group_name].SetLineColor(line_color)
                    _histograms[_group_name].SetFillColor(fill_color)
                    _histograms[_group_name].SetTitle(_group_name)
                    print("Set title:", _group_name)
            
            # _histograms[isSignal][-1].SetLineColor(line_color)
            # _histograms[isSignal][-1].SetFillColor(fill_color)
            # _histograms[isSignal][-1].SetLineWidth(2)
            # _histograms[isSignal][-1].SetMarkerSize(0)
            # _histograms[isSignal][-1].SetTitle(_group_name)
    
    
    # get maximum for the y-scale
    _histograms_values = list(_histograms.values())
    print(_histograms_values)
    y_max =_histograms_values[0].GetMaximum()
    for _h in _histograms_values:
        y_max = max(y_max,_h.GetMaximum())

    # sort histogram from min to max
    _histograms_background_entries = []
    _histograms_background_sorted = []
    for _h in _histograms_values:
        _histograms_background_entries.append(_h.Integral())
    _sorted_hist = np.argsort(_histograms_background_entries)
    for _idx in _sorted_hist:
        _histograms_background_sorted.append(_histograms_values[_idx])

    
    # read the binning if available:
    # if _histname in config["SetupBins"]:
    #     xrange_min = config["SetupBins"][_histname][0]
    #     xrange_max = config["SetupBins"][_histname][1]
    #     overflow =  bool(config["SetupBins"][_histname][3])
    # else:
    xrange_min =_histograms_values[0].GetXaxis().GetXmin()
    # xrange_max = _histograms["background"][0].GetXaxis().GetXmax()
    xrange_max = 53
    overflow =  True

    score_pass_finebin = config["mixing_hists"]["labels"]
    
    for bin_filled in [0, 1, 2]:
        # two loose region
        two_loose_hists = []
        for _h in _histograms_background_sorted:
            for bin_i, label in enumerate(score_pass_finebin):
                _h.GetXaxis().SetBinLabel(bin_i+1, ">"+str(label))
                _h.GetXaxis().SetTitle("")
            two_loose_hists.append(_h.ProjectionX(_h.GetTitle()+"_proj", 2+bin_filled, 2+bin_filled))
            print(two_loose_hists[-1].GetTitle(), two_loose_hists[-1].Integral())
        
        root_plot1D(
                l_hist = two_loose_hists,
                l_hist_overlay = two_loose_hists,
                outfile = output_path + "/" +  f"mixMC_plot_bin{bin_filled}.png",
                xrange = [xrange_min, xrange_max],
                yrange = (0.001,  1000*y_max),
                # yrange = (0.0,  1.5*y_max),
                logx = False, logy = True,
                include_overflow = overflow,
                xtitle = _histograms_values[0].GetXaxis().GetTitle(),
                ytitle = "events",
                xtitle_ratio = _histograms_values[0].GetXaxis().GetTitle(),
                ytitle_ratio = "%",
                centertitlex = True, centertitley = True,
                centerlabelx = False, centerlabely = False,
                gridx = False, gridy = False,
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
                ratio_mode = "percentage",
                logx_ratio = False, logy_ratio = False,
                yrange_ratio = (0, 1),
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
