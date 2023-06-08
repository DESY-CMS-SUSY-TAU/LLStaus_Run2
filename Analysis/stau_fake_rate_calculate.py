import os
import json
from argparse import ArgumentParser
import re

import ROOT
ROOT.gROOT.SetBatch(True)

from pepper import Config
import utils.utils as utils
from utils.hist_rebin import TH3Histogram

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

nominator = "Cut-008_n_valid_jets_jet_fake_tight"
denominator = "Cut-008_n_valid_jets_jet_fake_loose"

fake_rate_bins = {
    "jet_pt" : [20, 25, 30, 35, 40, 50, 70, 120, 1000],
    "jet_eta" : [-2.4, -1.5, 0, 1.5, 2.4],
    "jet_dxy" : [-20.0, -16.0, -10.4, -6.8, -3.2, 0, 3.2, 6.8, 10.4, 16.0, 20.0]
}

dirname = os.path.dirname(args.histfile[0])

hist_fake = {}
for region in [nominator, denominator]:
    file_path = dirname + "/" + region + ".root"
    file = ROOT.TFile.Open(str(file_path), 'read')
    hist_fake[region] = None
    for _group_idx, _group_name in enumerate(config["Labels"].keys()):
        if not _group_name in config["MC_bkgd"]:
            continue
        # Accumulate the dataset for the particular data group as specified in config “Labels”.
        for _idx, _histogram_data in enumerate(config["Labels"][_group_name]):
            # print(region, _histogram_data)
            hist = file.Get(_histogram_data)
            hist.SetDirectory(0)
            N = cutflow[_histogram_data]["all"]["Before cuts"]
            hist.Scale( (crosssections[_histogram_data] * config["luminosity"]) / N)
            if hist_fake[region] is None:
                hist_fake[region] = hist
            else:
                hist_fake[region].Add(hist)
    file.Close()
    # print(hist_fake[region].Integral())   
# exit()
# Divide histogram:
nominator_th3 = TH3Histogram(hist_fake[nominator],
                             fake_rate_bins["jet_pt"],
                             fake_rate_bins["jet_eta"],
                             fake_rate_bins["jet_dxy"])
nominator_th3_hist = nominator_th3.get_th3_histogram()
denominator_th3 = TH3Histogram(hist_fake[denominator],
                               fake_rate_bins["jet_pt"],
                               fake_rate_bins["jet_eta"],
                               fake_rate_bins["jet_dxy"])
denominator_th3_hist = denominator_th3.get_th3_histogram()

fake_sf = nominator_th3_hist
fake_sf.SetDirectory(0)
fake_sf.Divide(denominator_th3_hist)
output = ROOT.TFile(args.outdir+"/fake_rate_pt_eta_dxy.root", "RECREATE")
fake_sf.Write()
output.Close()

# Plots
# X - jet_pt
# Y - jet_eta
# Z - jet_dxy

hist_projection_XY_nom = nominator_th3_hist.Project3D("nom_xy")
hist_projection_XY_den = denominator_th3_hist.Project3D("den_xy")
hist_projection_XY_nom.Divide(hist_projection_XY_den)
output_xy = ROOT.TFile(args.outdir+"/fake_rate_pt_eta.root", "RECREATE")
hist_projection_XY_nom.Write()
output_xy.Close()

hist_projection_XZ_nom = nominator_th3_hist.Project3D("nom_xz")
hist_projection_XZ_den = denominator_th3_hist.Project3D("den_xz")
hist_projection_XZ_nom.Divide(hist_projection_XZ_den)
output_xz = ROOT.TFile(args.outdir+"/fake_rate_pt_dxy.root", "RECREATE")
hist_projection_XZ_nom.Write()
output_xz.Close()

# def print_hist(histogram):
#     for xbin in range(1, histogram.GetNbinsX() + 1):
#         for ybin in range(1, histogram.GetNbinsY() + 1):
#             for zbin in range(1, histogram.GetNbinsZ() + 1):
#                 bin_content = histogram.GetBinContent(xbin, ybin, zbin)
#                 print(f"Bin ({xbin}, {ybin}, {zbin}): {bin_content}")
            
# print_hist(nominator_th3_hist)
# print_hist(denominator_th3_hist)

hist_projection_X_nom = nominator_th3_hist.Project3D("nom_x")
hist_projection_X_den = denominator_th3_hist.Project3D("den_x")
# print(hist_projection_X_nom.Print("all"))
# print(hist_projection_X_den.Print("all"))
hist_projection_X_nom.Divide(hist_projection_X_den)
# print(hist_projection_X_nom.Print("all"))
output_x = ROOT.TFile(args.outdir+"/fake_rate_pt.root", "RECREATE")
hist_projection_X_nom.Write()
output_x.Close()

hist_projection_Y_nom = nominator_th3_hist.Project3D("nom_y")
hist_projection_Y_den = denominator_th3_hist.Project3D("den_y")
hist_projection_Y_nom.Divide(hist_projection_Y_den)
output_y = ROOT.TFile(args.outdir+"/fake_rate_eta.root", "RECREATE")
hist_projection_Y_nom.Write()
output_y.Close()

hist_projection_Z_nom = nominator_th3_hist.Project3D("nom_z")
hist_projection_Z_den = denominator_th3_hist.Project3D("den_z")
hist_projection_Z_nom.Divide(hist_projection_Z_den)
output_z = ROOT.TFile(args.outdir+"/fake_rate_dxy.root", "RECREATE")
hist_projection_Z_nom.Write()
output_z.Close()
