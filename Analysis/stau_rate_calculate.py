import os
import json
from argparse import ArgumentParser
import re
import numpy as np

import ROOT
ROOT.gROOT.SetBatch(True)

from pepper import Config
import utils.utils as utils
from utils.hist_rebin import TH3Histogram, th3_to_cumulative
ROOT.gInterpreter.Declare('#include "utils/histogram2d.cpp"')

def OverflowIntegralTHN(hist_th3):
    # Calculate integral of all TH3/TH2/TH1 histogram bins including overflow bins
    Integral = 0
    if hist_th3.GetDimension() == 1:
        for i in range(0, hist_th3.GetNbinsX()+2):
            Integral += hist_th3.GetBinContent(i)
        return Integral
    elif hist_th3.GetDimension() == 2:
        for i in range(0, hist_th3.GetNbinsX()+2):
            for j in range(0, hist_th3.GetNbinsY()+2):
                Integral += hist_th3.GetBinContent(i, j)
        return Integral
    elif hist_th3.GetDimension() == 3:
        for i in range(0, hist_th3.GetNbinsX()+2):
            for j in range(0, hist_th3.GetNbinsY()+2):
                for k in range(0, hist_th3.GetNbinsZ()+2):
                    Integral += hist_th3.GetBinContent(i, j, k)
        return Integral
    else:
        raise ValueError("Wrong dimension of histogram")
    
parser = ArgumentParser(
    description="The following script calculate fake rate for stau analysis.")
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


if config["fake_rate"]["mode"] == "ratio":
    
    nominator = config["fake_rate"]["nominator"]
    denominator = config["fake_rate"]["denominator"]

    dirname = os.path.dirname(args.histfile[0])
    isDATA = False
    hist_fake = {}
    for region, name in zip([nominator, denominator], ["nom", "denom"]):
        file_path = dirname + "/" + region[0] + ".root"
        file = ROOT.TFile.Open(str(file_path), 'read')
        hist_fake[name] = None
        for _group_idx, _group_name in enumerate(config["Labels"].keys()):
            # if (not _group_name in config["MC_bkgd"]) and (not _group_name in config["Signal_samples"]):
            #     continue
            # if not( _group_name in config["Data"]): continue
            # Accumulate the dataset for the particular data group as specified in config "Labels".
            for _idx, _histogram_data in enumerate(config["Labels"][_group_name]):
                print("File:",file_path)
                print("Open:",_histogram_data+region[1])
                hist = file.Get(_histogram_data+region[1])
                if not hist:
                    print("Warning: Histogram not found! ", end='')
                    print("Histogram->", file, _histogram_data+region[1])
                    continue
                hist.SetDirectory(0)

                # print("No scaling!")
                # N = cutflow[_histogram_data]["all"]["Before cuts"]
                # hist.Scale( (crosssections[_histogram_data] * config["luminosity"]) / N)

                if _group_name in config["Data"]:
                    isDATA = True
                    pass
                else:
                    if isDATA: raise("Can not combine data and MC")

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
                        hist.Scale( (crosssections[_histogram_data] * config["luminosity"]) / N)

                        
                if hist_fake[name] is None:
                    hist_fake[name] = hist
                else:
                    hist_fake[name].Add(hist)
        file.Close()
        # print(hist_fake[region].Integral())   

    nominator = config["fake_rate"]["nominator"][0]
    denominator = config["fake_rate"]["denominator"][0]

    print("Integral pre-rebin nom:", OverflowIntegralTHN(hist_fake["nom"]))
    print("Integral pre-rebin den:", OverflowIntegralTHN(hist_fake["denom"]))
    print("Divide histogram:")
    nominator_th3 = TH3Histogram(hist_fake["nom"], *list(config['fake_rate']["rate_bins"].values()))
    nominator_th3_hist = nominator_th3.get_rebinned_histogram()
    denominator_th3 = TH3Histogram(hist_fake["denom"], *list(config['fake_rate']["rate_bins"].values()))
    denominator_th3_hist = denominator_th3.get_rebinned_histogram()
    print("Integral post-rebin nom:", OverflowIntegralTHN(nominator_th3_hist))
    print("Integral post-rebin den:", OverflowIntegralTHN(denominator_th3_hist))
    fake_sf = nominator_th3_hist.Clone()
    fake_sf.SetDirectory(0)
    fake_sf.Divide(denominator_th3_hist)
    output = ROOT.TFile(args.outdir + f"/fake_rate_{'_'.join(list(config['fake_rate']['rate_bins'].keys()))}.root", "RECREATE")
    fake_sf.SetName(f"fake_rate_{'_'.join(list(config['fake_rate']['rate_bins'].keys()))}")
    fake_sf.Write()
    output.Close()

    # Plots
    # X - jet_pt
    # Y - jet_eta
    # Z - jet_dxy
    options = ""
    if config["fake_rate"]["NOF"]: options+="_NOF"
    if config["fake_rate"]["NUF"]: options+="_NUF"
    
    for name in config['fake_rate']["sf_project"]:
        project = config['fake_rate']["sf_project"][name][0]
        rebin_non_unifor = config['fake_rate']["sf_project"][name][1]
        print(name, project)
        hist_projection_nom = nominator_th3_hist.Project3D("nom_" + project + options)
        hist_projection_den = denominator_th3_hist.Project3D("den_" + project + options)
        print("Integral projection nom:", OverflowIntegralTHN(hist_projection_nom))
        print("Integral projection den:", OverflowIntegralTHN(hist_projection_den))
        # if hist_projection_nom[-1].GetDimension() == 1:
        #     hist_projection_nom[-1].Print("all")
        #     hist_projection_den[-1].Print("all")
        if rebin_non_unifor:
            assert(isinstance(rebin_non_unifor, dict))
            assert("y_axis" in rebin_non_unifor)
            assert("x_axis" in rebin_non_unifor)
            
            x_axis = rebin_non_unifor["x_axis"]
            y_axis =  np.array(rebin_non_unifor["y_axis"], dtype=np.double)
            xmin, xmax = float(x_axis[0][0]), float(x_axis[0][-1])
            nom_nonunif = ROOT.Histogram_2D("nom_nonunif_" + project + options, y_axis, xmin, xmax)
            den_nonunif = ROOT.Histogram_2D("den_nonunif_" + project + options, y_axis, xmin, xmax)
            for i in range(len(x_axis)):
                nom_nonunif.add_x_binning_by_index(i, np.array(x_axis[i], dtype=np.double))
                den_nonunif.add_x_binning_by_index(i, np.array(x_axis[i], dtype=np.double))
            
            try:
                nom_nonunif.th2d_add(hist_projection_nom)
                den_nonunif.th2d_add(hist_projection_den)
            except:
                print("Error is inside project nonunif histograms")
                print("Error: can not add histogram because of the inconsistency")
                del nom_nonunif
                del den_nonunif
                exit()
                
            nom_nonunif.divide(den_nonunif)
            weight_hist = nom_nonunif.get_weights_th2d_simpl("hist_weight","hist_weight")
            # weight_hist.Print("all")
            
            output = ROOT.TFile(args.outdir+f"/fake_rate_{name}.root", "RECREATE")
            weight_hist.SetName(f"fake_rate_{name}")
            weight_hist.Write()
            output.Close()
            
            del nom_nonunif
            del den_nonunif
            del weight_hist

        else:
            hist_projection_nom.Divide(hist_projection_den)
            output = ROOT.TFile(args.outdir+f"/fake_rate_{name}.root", "RECREATE")
            hist_projection_nom.SetName(f"fake_rate_{name}")
            hist_projection_nom.Write()
            output.Close()


if config["fake_rate"]["mode"] == "score-dim": # calculated for data
    
    dirname = os.path.dirname(args.histfile[0])
    
    hist_fake = None

    file_path = dirname + "/" + config["fake_rate"]["histogram"] + ".root"
    file = ROOT.TFile.Open(str(file_path), 'read')
    print(file_path)
    for _group_idx, _group_name in enumerate(config["Labels"].keys()):
        
    #     if not ((_group_name in config["MC_bkgd"]) or
    #             (_group_name in config["Signal_samples"]) or
    #             (_group_name in config["Data"])):
    #         continue
    
        # Only data is taken for the calculation of fake rate
        if not( _group_name in config["Data"]): continue
        
        # Accumulate the dataset for the particular data group as specified in config “Labels”.
        for _idx, _histogram_data in enumerate(config["Labels"][_group_name]):
            print("SF data:", _histogram_data)
            hist = file.Get(_histogram_data)
            hist.SetDirectory(0)
            hist.Sumw2()
            if not _group_name in config["Data"]:
                N = cutflow[_histogram_data]["all"]["Before cuts"]
                hist.Scale( (crosssections[_histogram_data] * config["luminosity"]) / N)
            if hist_fake is None:
                hist_fake = hist
            else:
                hist_fake.Add(hist)
    file.Close()
    
    print("Integral pre-rebin nom:", OverflowIntegralTHN(hist_fake))
    th3_hist = hist_fake.Clone()
    th3_to_cumulative = th3_to_cumulative(th3_hist, axis_to_integrate=0)
    th3_to_cumulative_rate = th3_to_cumulative.Clone("th3_to_cumulative_rate")
    print("Integral post-rebin nom:", OverflowIntegralTHN(th3_to_cumulative_rate))
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~test
    # prob_projection = th3_hist.Project3D("proj_x")
    # print("From TH3")
    # print(prob_projection.Print("all"))
    
    # projected_hist_path = "Cut-008_two_valid_jets_fake_rate_hist_probonly.root"
    # file_path = dirname + "/" + projected_hist_path
    # file = ROOT.TFile.Open(str(file_path), 'read')
    # proj_origin = file.Get("SingleMuon")
    # print("From prob hist")
    # print(proj_origin.Print("all"))
    
    # print("Comulative!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    
    # prob_projection = th3_to_cumulative.Project3D("proj_x")
    # print("From TH3")
    # print(prob_projection.Print("all"))
    
    # projected_hist_path = "Cut-008_two_valid_jets_fake_rate_hist_probonly.root"
    # file_path = dirname + "/" + projected_hist_path
    # proj_origin_comul = proj_origin.GetCumulative(False)
    # print("From prob hist")
    # print(proj_origin_comul.Print("all"))
    
    # exit()
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~test
    
    # Perform bin-by-bin division by the first bin in x for all bins in y and z
    n_bins_x = th3_to_cumulative_rate.GetNbinsX()
    n_bins_y = th3_to_cumulative_rate.GetNbinsY()
    n_bins_z = th3_to_cumulative_rate.GetNbinsZ()
    
    for bin_y in range(0, n_bins_y + 2):
        for bin_z in range(0, n_bins_z + 2):
            for bin_x in range(0, n_bins_x + 2):
                content_nom = th3_to_cumulative.GetBinContent(bin_x, bin_y, bin_z)
                content_nom_err = th3_to_cumulative.GetBinError(bin_x, bin_y, bin_z)
                content_den = th3_to_cumulative.GetBinContent(1, bin_y, bin_z)
                content_den_err = th3_to_cumulative.GetBinError(1, bin_y, bin_z)
                rate = (content_nom/content_den) if content_den != 0 else 0
                th3_to_cumulative_rate.SetBinContent(bin_x, bin_y, bin_z, rate)
                th3_to_cumulative_rate.SetBinError(bin_x, bin_y, bin_z, 0)
                if content_den != 0 and content_nom != 0:
                    error = rate * np.sqrt((content_nom_err/content_nom)**2+(content_den_err/content_den)**2)
                    th3_to_cumulative_rate.SetBinError(bin_x, bin_y, bin_z, error)
                if content_nom > content_den and content_den!=0:
                    raise("Detected rate > 1.0:", bin_x, bin_y, bin_z,
                          "nom:",content_nom, "den:", content_den)
                
    print("3D scale factors")
    output = ROOT.TFile(args.outdir + f"/fake_rate_score_{'_'.join(list(config['fake_rate']['rate_bins'].keys()))}.root", "RECREATE")
    th3_to_cumulative_rate.SetName(f"fake_rate_score_{'_'.join(list(config['fake_rate']['rate_bins'].keys()))}")
    th3_to_cumulative_rate.Write()
    output.Close()
    
    print("Histogram ranges")
    print("X:", th3_to_cumulative_rate.GetXaxis().GetXmin(), th3_to_cumulative_rate.GetXaxis().GetXmax())
    print("Y:", th3_to_cumulative_rate.GetYaxis().GetXmin(), th3_to_cumulative_rate.GetYaxis().GetXmax())
    print("Z:", th3_to_cumulative_rate.GetZaxis().GetXmin(), th3_to_cumulative_rate.GetZaxis().GetXmax())
    
    for name in config['fake_rate']["sf_project"]:
        print("Creatintg projection for", name)
        hist_projection = []
        hist_projection_rate = []
        project = config['fake_rate']["sf_project"][name]
        # underflow/overflow are included by default in the projection 
        # To exclude underflow and/or overflow (for both axis in case of a projection to a 1D histogram) 
        # use option "NUF" and/or "NOF" With SetRange() you can have all bins except underflow/overflow 
        # only if you set the axis bit range as following after having called SetRange: axis->SetRange(1, axis->GetNbins());
        hist_projection.append(th3_to_cumulative.Project3D("nom_"+project))
        hist_projection_rate.append(hist_projection[-1].Clone("rate_"+project))
        if hist_projection[-1].GetDimension() == 2:
            n_bins_x = hist_projection[-1].GetNbinsX()
            n_bins_y = hist_projection[-1].GetNbinsY()
            for bin_x in range(0, n_bins_x + 2):
                for bin_y in range(0, n_bins_y + 2):
                    content_nom = hist_projection[-1].GetBinContent(bin_x, bin_y)
                    content_den = hist_projection[-1].GetBinContent(1, bin_y)
                    content_nom_error = hist_projection[-1].GetBinError(bin_x, bin_y)
                    content_den_error = hist_projection[-1].GetBinError(1, bin_y)
                    rate = (content_nom/content_den) if content_den != 0 else 0
                    hist_projection_rate[-1].SetBinContent(bin_x, bin_y, rate)
                    hist_projection_rate[-1].SetBinError(bin_x, bin_y, 0)
                    th3_to_cumulative_rate.SetBinError(bin_x, bin_y, 0)
                    if content_den != 0 and content_nom != 0:
                        error = rate * np.sqrt((content_nom_err/content_nom)**2+(content_den_err/content_den)**2)
                        hist_projection_rate[-1].SetBinError(bin_x, bin_y, error)
        elif hist_projection[-1].GetDimension() == 1:
            n_bins_x = hist_projection[-1].GetNbinsX()
            for bin_x in range(0, n_bins_x + 2):
                content_nom = hist_projection[-1].GetBinContent(bin_x)
                content_den = hist_projection[-1].GetBinContent(1)
                rate = (content_nom/content_den) if content_den != 0 else 0
                content_nom_error = hist_projection[-1].GetBinError(bin_x)
                content_den_error = hist_projection[-1].GetBinError(1)
                hist_projection_rate[-1].SetBinContent(bin_x, rate)
                th3_to_cumulative_rate.SetBinError(bin_x, 0)
                if content_den != 0 and content_nom != 0:
                    error = rate * np.sqrt((content_nom_err/content_nom)**2+(content_den_err/content_den)**2)
                    hist_projection_rate[-1].SetBinError(bin_x, error)
        else:
            raise ValueError("Wrong dimension of the histogram")

        output = ROOT.TFile(args.outdir+f"/fake_rate_{name}.root", "RECREATE")
        hist_projection_rate[-1].SetName(f"fake_rate_{name}")
        hist_projection_rate[-1].Write()
        output.Close()