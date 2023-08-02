import os
import json
from argparse import ArgumentParser
import re

import ROOT
ROOT.gROOT.SetBatch(True)

from pepper import Config
import utils.utils as utils
from utils.hist_rebin import TH3Histogram, th3_to_cumulative

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

    hist_fake = {}
    for region in [nominator, denominator]:
        file_path = dirname + "/" + region + ".root"
        file = ROOT.TFile.Open(str(file_path), 'read')
        hist_fake[region] = None
        for _group_idx, _group_name in enumerate(config["Labels"].keys()):
            if (not _group_name in config["MC_bkgd"]) and (not _group_name in config["Signal_samples"]):
                continue
            # Accumulate the dataset for the particular data group as specified in config "Labels".
            for _idx, _histogram_data in enumerate(config["Labels"][_group_name]):
                print(file_path, _histogram_data)
                hist = file.Get(_histogram_data)
                hist.SetDirectory(0)
                # print("No scaling!")
                N = cutflow[_histogram_data]["all"]["Before cuts"]
                hist.Scale( (crosssections[_histogram_data] * config["luminosity"]) / N)
                if hist_fake[region] is None:
                    hist_fake[region] = hist
                else:
                    hist_fake[region].Add(hist)
        file.Close()
        # print(hist_fake[region].Integral())   

    print("Integral pre-rebin nom:", OverflowIntegralTHN(hist_fake[nominator]))
    print("Integral pre-rebin den:", OverflowIntegralTHN(hist_fake[denominator]))
    print("Divide histogram:")
    nominator_th3 = TH3Histogram(hist_fake[nominator], *list(config["rate_bins"].values()))
    nominator_th3_hist = nominator_th3.get_rebinned_histogram()
    denominator_th3 = TH3Histogram(hist_fake[denominator], *list(config["rate_bins"].values()))
    denominator_th3_hist = denominator_th3.get_rebinned_histogram()
    print("Integral post-rebin nom:", OverflowIntegralTHN(nominator_th3_hist))
    print("Integral post-rebin den:", OverflowIntegralTHN(denominator_th3_hist))
    fake_sf = nominator_th3_hist.Clone()
    fake_sf.SetDirectory(0)
    fake_sf.Divide(denominator_th3_hist)
    output = ROOT.TFile(args.outdir + f"/fake_rate_{'_'.join(list(config['rate_bins'].keys()))}.root", "RECREATE")
    fake_sf.SetName(f"fake_rate_{'_'.join(list(config['rate_bins'].keys()))}")
    fake_sf.Write()
    output.Close()

    # Plots
    # X - jet_pt
    # Y - jet_eta
    # Z - jet_dxy
    hist_projection_nom = []
    hist_projection_den = []

    for name in config["sf_project"]:
        project = config["sf_project"][name]
        print(name, project)
        hist_projection_nom.append(nominator_th3_hist.Project3D("nom_"+project))
        hist_projection_den.append(denominator_th3_hist.Project3D("den_"+project))
        print("Integral projection nom:", OverflowIntegralTHN(hist_projection_nom[-1]))
        print("Integral projection den:", OverflowIntegralTHN(hist_projection_den[-1]))
        if hist_projection_nom[-1].GetDimension() == 1:
            hist_projection_nom[-1].Print("all")
            hist_projection_den[-1].Print("all")
        hist_projection_nom[-1].Divide(hist_projection_den[-1])
        output = ROOT.TFile(args.outdir+f"/fake_rate_{name}.root", "RECREATE")
        hist_projection_nom[-1].SetName(f"fake_rate_{name}")
        hist_projection_nom[-1].Write()
        output.Close()


if config["fake_rate"]["mode"] == "score-dim": # calculated for data
    
    dirname = os.path.dirname(args.histfile[0])
    
    hist_fake = None

    file_path = dirname + "/" + config["fake_rate"]["histogram"] + ".root"
    file = ROOT.TFile.Open(str(file_path), 'read')
    for _group_idx, _group_name in enumerate(config["Labels"].keys()):
        if not ((_group_name in config["MC_bkgd"]) or
                (_group_name in config["Signal_samples"]) or
                (_group_name in config["Data"])):
            continue
        # Accumulate the dataset for the particular data group as specified in config “Labels”.
        for _idx, _histogram_data in enumerate(config["Labels"][_group_name]):
            print(file_path, _histogram_data)
            hist = file.Get(_histogram_data)
            hist.SetDirectory(0)
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
    
    # Perform bin-by-bin division by the first bin in x for all bins in y and z
    n_bins_x = th3_to_cumulative_rate.GetNbinsX()
    n_bins_y = th3_to_cumulative_rate.GetNbinsY()
    n_bins_z = th3_to_cumulative_rate.GetNbinsZ()
    
    for bin_y in range(0, n_bins_y + 1):
        for bin_z in range(0, n_bins_z + 1):
            for bin_x in range(0, n_bins_x + 1):
                content_nom = th3_to_cumulative.GetBinContent(bin_x, bin_y, bin_z)
                content_den = th3_to_cumulative.GetBinContent(1, bin_y, bin_z)
                th3_to_cumulative_rate.SetBinContent(bin_x, bin_y, bin_z, (content_nom/content_den) if content_den != 0 else 0)
                
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
        hist_projection.append(th3_to_cumulative.Project3D("nom_"+project))
        hist_projection_rate.append(hist_projection[-1].Clone("rate_"+project))
        if hist_projection[-1].GetDimension() == 2:
            n_bins_x = hist_projection[-1].GetNbinsX()
            n_bins_y = hist_projection[-1].GetNbinsY()
            for bin_x in range(0, n_bins_x + 1):
                for bin_y in range(0, n_bins_y + 1):
                    content_nom = hist_projection[-1].GetBinContent(bin_x, bin_y)
                    content_den = hist_projection[-1].GetBinContent(1, bin_y)
                    hist_projection_rate[-1].SetBinContent(bin_x, bin_y, content_nom/content_den if content_den != 0 else 0)
        elif hist_projection[-1].GetDimension() == 1:
            n_bins_x = hist_projection[-1].GetNbinsX()
            for bin_x in range(0, n_bins_x + 1):
                content_nom = hist_projection[-1].GetBinContent(bin_x)
                content_den = hist_projection[-1].GetBinContent(1)
                hist_projection_rate[-1].SetBinContent(bin_x, content_nom/content_den if content_den != 0 else 0)
        else:
            raise ValueError("Wrong dimension of the histogram")

        output = ROOT.TFile(args.outdir+f"/fake_rate_{name}.root", "RECREATE")
        hist_projection_rate[-1].SetName(f"fake_rate_{name}")
        hist_projection_rate[-1].Write()
        output.Close()