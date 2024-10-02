import os
import json
from argparse import ArgumentParser
import re
import ROOT
import itertools
import numpy as np
import array

from pepper import Config

import sys
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(current_dir, '..'))
sys.path.append(parent_dir)
# from utils.plotter import plot1D, plot2D
# from utils.utils import ColorIterator, root_plot1D, root_plot2D

output_dir = "plots_spectrum_deeptau"
ext = "png"
spectrum_path_init = "/nfs/dust/cms/user/mykytaua/dataDeepTau/full_tuples_spectrum"
spectrum_files_init = [
    "DY1JetsToLL_M50.root",                                                       
    "DY2JetsToLL_M50.root",                                                          
    "DY3JetsToLL_M50.root",                                                         
    "DY4JetsToLL_M50.root",                                                       
    "DYJetsToLL_M-50-amcatnloFXFX_ext2.root",                                
    "DYJetsToLL_M-50-amcatnloFXFX.root",                              
    "DYJetsToLL_M-50_HT-100to200.root",                
    "DYJetsToLL_M-50_HT-1200to2500.root",       
    "DYJetsToLL_M-50_HT-200to400.root",                                           
    "DYJetsToLL_M-50_HT-2500toInf.root",                                          
    "DYJetsToLL_M-50_HT-400to600_ext2.root",                                          
    "DYJetsToLL_M-50_HT-400to600.root",                                                
    "DYJetsToLL_M-50_HT-600to800.root",                                                 
    "DYJetsToLL_M-50_HT-70to100.root",                                            
    "DYJetsToLL_M-50_HT-800to1200.root",                                          
    "DYJetsToLL_M-50.root",     
    # "Embedded_MuTau_Run2018A.root",
    # "Embedded_MuTau_Run2018B.root",
    # "Embedded_MuTau_Run2018C.root",
    # "Embedded_MuTau_Run2018D.root",
    # "GluGluHToTauTau_M125.root",
    # "Muminus_Pt1000-gun_ext1.root",
    # "Muplus_Pt1000-gun_ext1.root",
    "QCD_HT1000to1500.root",
    "QCD_HT100to200.root",
    "QCD_HT1500to2000.root",
    "QCD_HT2000toInf.root",
    "QCD_HT200to300.root",
    "QCD_HT300to500.root",
    "QCD_HT500to700.root",
    "QCD_HT50to100.root",                     
    "QCD_HT700to1000.root",
    "QCD_Pt_1000to1400.root",
    "QCD_Pt_120to170.root",
    "QCD_Pt_1400to1800_ext1.root",
    "QCD_Pt_1400to1800.root",
    "QCD_Pt_15to30_ext1.root",
    "QCD_Pt_15to7000_Flat2018_ext1.root",
    "QCD_Pt_15to7000_Flat_ext1.root",
    "QCD_Pt_170to300.root",
    "QCD_Pt_1800to2400_ext1.root",
    "QCD_Pt_1800to2400.root",
    "QCD_Pt_2400to3200_ext1.root",
    "QCD_Pt_2400to3200.root",
    "QCD_Pt_300to470.root",
    "QCD_Pt_30to50.root",
    "QCD_Pt_3200toInf.root",                               
    "QCD_Pt_470to600_ext1.root",
    "QCD_Pt_470to600.root",
    "QCD_Pt_50to80_ext1.root",
    "QCD_Pt_50to80.root",
    "QCD_Pt_600to800.root",
    "QCD_Pt_800to1000_ext1.root",
    "QCD_Pt_80to120.root",
    "ZprimeToEE_M-1500.root",
    "ZprimeToEE_M-2000.root",
    "ZprimeToEE_M-2500.root",
    "ZprimeToEE_M-3000.root",
    "ZprimeToEE_M-3500.root",
    "ZprimeToEE_M-4000.root",
    "ZprimeToEE_M-1000.root",
    # "ZprimeToMuMu_M-5000_ext1.root",
    # "ZprimeToTauTau_M-1000.root",
    # "ZprimeToTauTau_M-1500.root",
    # "ZprimeToTauTau_M-2000.root",
    # "ZprimeToTauTau_M-2500.root",
    # "ZprimeToTauTau_M-3000.root",
    # "ZprimeToTauTau_M-3500.root",
    # "ZprimeToTauTau_M-4000.root",
    "W1JetsToLNu.root",
    "W2JetsToLNu.root",
    "W3JetsToLNu.root",
    "W4JetsToLNu.root",
    "WJetsToLNu_HT-100To200.root",
    "WJetsToLNu_HT-1200To2500.root",
    "WJetsToLNu_HT-200To400.root",
    "WJetsToLNu_HT-2500ToInf.root",
    "WJetsToLNu_HT-400To600.root",
    "WJetsToLNu_HT-600To800.root",
    "WJetsToLNu_HT-800To1200.root",
    "WJetsToLNu.root",
    "WminusHToTauTau_M125.root",
    "WplusHToTauTau_M125.root",
    "ZHToTauTau_M125.root",
    "TauGunPt-15to3000.root",
    "VBFHToTauTau.root",
    "GluGluToHHTo2B2Tau.root",
    # "ttHToTauTau_M125.root",
    "TTJets_ext1.root",
    "TTTo2L2Nu.root",
    "TTToHadronic_ext2.root",
    "TTToHadronic.root",
    "TTToSemiLeptonic_ext3.root",
    # # "TTToSemiLeptonic.root",
]

spectrum_file_final = "/nfs/dust/cms/user/mykytaua/dataDeepTau/DeepTauTraining/ShuffleMergeSpectral_TrainingSpectrum/ShuffleMergeSpectral_trainingSamples-2_rerun.root"

pt_binning = [20, 30, 40, 50, 60, 70, 100, 200, 500, 1000]
eta_binning = [0.0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.5]

left_margin = 0.13
bottom_margin = 0.13

# sum up all init histograms over 4 types of tau .Add() them together
initial_hists = {
    "eta_pt_hist_e" : None,
    "eta_pt_hist_mu" : None,
    "eta_pt_hist_tau" : None,
    "eta_pt_hist_jet" : None,
} # for each type of tau

ROOT.gStyle.SetTitleSize(0.07, "yz")
ROOT.gStyle.SetLabelSize(0.07, "yz")
ROOT.gStyle.SetTitleSize(0.09, "x")
ROOT.gStyle.SetLabelSize(0.09, "x")
ROOT.gStyle.SetLegendTextSize(0.07)

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

for file in spectrum_files_init:
    f = ROOT.TFile.Open(spectrum_path_init + "/" + file, "read")
    for key in initial_hists:
        hist = f.Get(key)
        if hist:
            hist.SetDirectory(0)
            if initial_hists[key] is None:
                initial_hists[key] = hist.Clone()
                initial_hists[key].SetDirectory(0)
            else:
                initial_hists[key].Add(hist)
            
# show the number of entried for each type:
for key in initial_hists:
    print(f"{key}: {initial_hists[key].GetEntries()}")


########################################## Pt ##########################################

# first project the 2D histograms to 1D histograms
# plot pt, eta for each type of tau (tau: orange, mu: red, e: blue, jet: green)
projected_pt_init = {}
projected_eta_init = {}
for key in initial_hists:
    if "pt" in key:
        projected_pt_init[key] = initial_hists[key].ProjectionY()
        projected_pt_init[key].SetDirectory(0)
    if "eta" in key:
        projected_eta_init[key] = initial_hists[key].ProjectionX()
        projected_eta_init[key].SetDirectory(0)
# plot the pt and eta for each type of tau  (tau: orange, mu: red, e: blue, jet: green)
color_hists = {
        "eta_pt_hist_jet" : 3,
        "eta_pt_hist_e" : 4,
        "eta_pt_hist_mu" : 2,
        "eta_pt_hist_tau" : 92
    }

# create a canvas
c = ROOT.TCanvas("c", "c", 1400, 700)
c.SetLogy()
c.SetLogx()
c.SetLeftMargin(left_margin)
c.SetBottomMargin(bottom_margin)
for i, key in enumerate(projected_pt_init):
    projected_pt_init[key].SetLineColor(color_hists[key])
    projected_pt_init[key].GetXaxis().SetRangeUser(20, 1000)
    projected_pt_init[key].GetYaxis().SetRangeUser(10, 5E7)
    # projected_pt_init[key].GetYaxis().SetRangeUser(0, 1000)
    new_bins = array.array('d', pt_binning)
    # new_bins = array.array('d', [20, 30, 40, 50, 60, 70, 200, 300, 400, 500, 600, 700, 800, 900, 1000])
    # new_bins = array.array('d', [i for i in range(20, 1000, 20)])
    projected_pt_init[key] = projected_pt_init[key].Rebin(len(new_bins)-1, f"{key}_rebin", new_bins)
    # divide bin content by bin width
    for i in range(1, projected_pt_init[key].GetNbinsX()+1):
        bin_content = projected_pt_init[key].GetBinContent(i)
        bin_width = projected_pt_init[key].GetBinWidth(i)
        projected_pt_init[key].SetBinContent(i, bin_content/bin_width)
        projected_pt_init[key].SetBinError(i, 0)
    projected_pt_init[key].SetDirectory(0)
    projected_pt_init[key].Print("all")
    projected_pt_init[key].GetXaxis().SetTitle("p_{T} [GeV]")
    # projected_pt_init[key].GetYaxis().SetNdivisions(505)
    projected_pt_init[key].GetYaxis().SetTitle("#Delta N-entries/#Delta p_{T} [1/GeV]")
    projected_pt_init[key].SetTitle("Initial Data")
    projected_pt_init[key].SetStats(0)
    projected_pt_init[key].SetLineWidth(2)
    projected_pt_init[key].SetFillColor(0)
    projected_pt_init[key].SetFillStyle(0)
    projected_pt_init[key].GetXaxis().SetTitleOffset(1.3)
    projected_pt_init[key].GetYaxis().SetTitleOffset(0.8)
    if i == 0:
        projected_pt_init[key].Draw()
    else:
        projected_pt_init[key].Draw("same")
    projected_pt_init[key].SetDirectory(0)

# add legend
leg = ROOT.TLegend(0.7, 0.6, 0.9, 0.9)
leg.AddEntry(projected_pt_init["eta_pt_hist_jet"],"#tau_{j}", "l")
leg.AddEntry(projected_pt_init["eta_pt_hist_e"],"#tau_{e}", "l")
leg.AddEntry(projected_pt_init["eta_pt_hist_mu"], "#tau_{#mu}", "l")
leg.AddEntry(projected_pt_init["eta_pt_hist_tau"], "#tau_{h}", "l")
leg.Draw()

c.SaveAs(f"{output_dir}/{key}_pt.{ext}")
# do the same plot for spectrum_file_final
initial_hists = {
    "eta_pt_hist_e" : None,
    "eta_pt_hist_mu" : None,
    "eta_pt_hist_tau" : None,
    "eta_pt_hist_jet" : None,
} # for each type of tau

f = ROOT.TFile.Open(spectrum_file_final, "read")
for key in initial_hists:
    hist = f.Get(key)
    if hist:
        hist.SetDirectory(0)
        if initial_hists[key] is None:
            initial_hists[key] = hist.Clone()
            initial_hists[key].SetDirectory(0)
        else:
            initial_hists[key].Add(hist)


# show the number of entried for each type:
for key in initial_hists:
    print(f"{key}: {initial_hists[key].GetEntries()}")

# first project the 2D histograms to 1D histograms
# plot pt, eta for each type of tau (tau: orange, mu: red, e: blue, jet: green)
projected_pt = {}
projected_eta = {}
for key in initial_hists:
    if "pt" in key:
        projected_pt[key] = initial_hists[key].ProjectionY()
    if "eta" in key:
        projected_eta[key] = initial_hists[key].ProjectionX()
        

# create a canvas
c = ROOT.TCanvas("c", "c", 1400, 700)
c.SetLogy()
c.SetLogx()
c.SetLeftMargin(left_margin)
c.SetBottomMargin(bottom_margin)
for i, key in enumerate(projected_pt):
    projected_pt[key].SetLineColor(color_hists[key])
    projected_pt[key].GetXaxis().SetRangeUser(20, 1000)
    projected_pt[key].GetYaxis().SetRangeUser(10, 5E7)
    # projected_pt[key].GetYaxis().SetRangeUser(0, 1000)
    new_bins = array.array('d', pt_binning)
    # new_bins = array.array('d', [20, 30, 40, 50, 60, 70, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])
    # new_bins = array.array('d', [i for i in range(20, 1000, 20)])
    projected_pt[key] = projected_pt[key].Rebin(len(new_bins)-1, f"{key}_rebin", new_bins)
    # divide bin content by bin width
    for i in range(1, projected_pt[key].GetNbinsX()+1):
        bin_content = projected_pt[key].GetBinContent(i)
        bin_width = projected_pt[key].GetBinWidth(i)
        projected_pt[key].SetBinContent(i, bin_content/bin_width)
        projected_pt[key].SetBinError(i, 0)
    projected_pt[key].SetDirectory(0)
    projected_pt[key].Print("all")
    projected_pt[key].GetXaxis().SetTitle("p_{T} [GeV]")
    # projected_pt[key].GetYaxis().SetNdivisions(505)
    projected_pt[key].GetYaxis().SetTitle("#DeltaN-entries/#Deltap_{T} [1/GeV]")
    projected_pt[key].SetTitle("Data after S&M sampling")
    projected_pt[key].SetStats(0)
    projected_pt[key].SetLineWidth(2)
    projected_pt[key].SetFillColor(0)
    projected_pt[key].SetFillStyle(0)
    projected_pt[key].GetXaxis().SetTitleOffset(1.3)
    projected_pt[key].GetYaxis().SetTitleOffset(0.8)
    if i == 0:
        projected_pt[key].Draw()
    else:
        projected_pt[key].Draw("same")
    projected_pt[key].SetDirectory(0)

# add legend
leg = ROOT.TLegend(0.7, 0.6, 0.9, 0.9)
leg.AddEntry(projected_pt["eta_pt_hist_jet"],"#tau_{j}", "l")
leg.AddEntry(projected_pt["eta_pt_hist_e"],"#tau_{e}", "l")
leg.AddEntry(projected_pt["eta_pt_hist_mu"], "#tau_{#mu}", "l")
leg.AddEntry(projected_pt["eta_pt_hist_tau"], "#tau_{h}", "l")
leg.Draw()

c.SaveAs(f"{output_dir}/{key}_pt_final.{ext}")

########################################## Eta ##########################################

c2 = ROOT.TCanvas("c2", "c2", 1400, 700)
c2.SetLogy()
c2.SetLeftMargin(left_margin)
c2.SetBottomMargin(bottom_margin)
for i, key in enumerate(projected_eta_init):
    projected_eta_init[key].SetLineColor(color_hists[key])
    # projected_eta_init[key].GetXaxis().SetRangeUser(-2.4, 2.4)
    projected_eta_init[key].GetYaxis().SetRangeUser(1E6, 1E10)
    new_bins = array.array('d', eta_binning)
    projected_eta_init[key] = projected_eta_init[key].Rebin(len(new_bins)-1, f"{key}_rebin", new_bins)
    # divide bin content by bin width
    for i in range(1, projected_eta_init[key].GetNbinsX()+1):
        bin_content = projected_eta_init[key].GetBinContent(i)
        bin_width = projected_eta_init[key].GetBinWidth(i)
        projected_eta_init[key].SetBinContent(i, bin_content/bin_width)
        projected_eta_init[key].SetBinError(i, 0)
    projected_eta_init[key].SetDirectory(0)
    projected_eta_init[key].Print("all")
    projected_eta_init[key].GetXaxis().SetTitle("#eta")

    projected_eta_init[key].GetYaxis().SetTitle("#DeltaN-entries/#Delta#eta")
    projected_eta_init[key].SetTitle("Initial Data")
    projected_eta_init[key].SetStats(0)
    projected_eta_init[key].SetLineWidth(2)
    projected_eta_init[key].SetFillColor(0)
    projected_eta_init[key].SetFillStyle(0)
    projected_eta_init[key].GetXaxis().SetTitleOffset(1.3)
    projected_eta_init[key].GetYaxis().SetTitleOffset(0.8)
    if i == 0:
        projected_eta_init[key].Draw()
    else:
        projected_eta_init[key].Draw("same")
    projected_eta_init[key].SetDirectory(0)

# add legend
leg = ROOT.TLegend(0.7, 0.6, 0.9, 0.9)
leg.AddEntry(projected_eta_init["eta_pt_hist_jet"],"#tau_{j}", "l")
leg.AddEntry(projected_eta_init["eta_pt_hist_e"],"#tau_{e}", "l")
leg.AddEntry(projected_eta_init["eta_pt_hist_mu"], "#tau_{#mu}", "l")
leg.AddEntry(projected_eta_init["eta_pt_hist_tau"], "#tau_{h}", "l")
leg.Draw()

c2.SaveAs(f"{output_dir}/{key}_eta.{ext}")

# create a canvas
c2 = ROOT.TCanvas("c2", "c2", 1400, 700)
c2.SetLogy()
c2.SetLeftMargin(left_margin)
c2.SetBottomMargin(bottom_margin)
for i, key in enumerate(projected_eta):
    projected_eta[key].SetLineColor(color_hists[key])
    # projected_eta[key].GetXaxis().SetRangeUser(-2.4, 2.4)
    projected_eta[key].GetYaxis().SetRangeUser(1E6, 1E10)
    new_bins = array.array('d', eta_binning)
    projected_eta[key] = projected_eta[key].Rebin(len(new_bins)-1, f"{key}_rebin", new_bins)
    # divide bin content by bin width
    for i in range(1, projected_eta[key].GetNbinsX()+1):
        bin_content = projected_eta[key].GetBinContent(i)
        bin_width = projected_eta[key].GetBinWidth(i)
        projected_eta[key].SetBinContent(i, bin_content/bin_width)
        projected_eta[key].SetBinError(i, 0)
    projected_eta[key].SetDirectory(0)
    projected_eta[key].Print("all")
    projected_eta[key].GetXaxis().SetTitle("#eta")
    # projected_eta[key].GetYaxis().SetNdivisions(505)
    projected_eta[key].GetYaxis().SetTitle("#DeltaN-entries/#Delta#eta")
    projected_eta[key].SetTitle("Data after S&M sampling")
    projected_eta[key].SetStats(0)
    projected_eta[key].SetLineWidth(2)
    projected_eta[key].SetFillColor(0)
    projected_eta[key].SetFillStyle(0)
    projected_eta[key].GetXaxis().SetTitleOffset(1.3)
    projected_eta[key].GetYaxis().SetTitleOffset(0.8)
    if i == 0:
        projected_eta[key].Draw()
    else:
        projected_eta[key].Draw("same")
    projected_eta[key].SetDirectory(0)

# add legend
leg = ROOT.TLegend(0.7, 0.6, 0.9, 0.9)
leg.AddEntry(projected_eta["eta_pt_hist_jet"],"#tau_{j}", "l")
leg.AddEntry(projected_eta["eta_pt_hist_e"],"#tau_{e}", "l")
leg.AddEntry(projected_eta["eta_pt_hist_mu"], "#tau_{#mu}", "l")
leg.AddEntry(projected_eta["eta_pt_hist_tau"], "#tau_{h}", "l")
leg.Draw()

c2.SaveAs(f"{output_dir}/{key}_eta_final.{ext}")