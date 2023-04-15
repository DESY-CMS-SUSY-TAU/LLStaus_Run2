#!/usr/bin/env python3

import argparse
import json
import os
import yaml
import functools
import operator

import ROOT

import utils.commonutils as cmut


def main() :
    
    # Argument parser
    parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument(
        "--config",
        help = "Configuration file",
        type = str,
        required = True,
    )
    
    # Parse arguments
    args = parser.parse_args()
    #d_args = vars(args)
    
    
    d_config = cmut.load_config(args.config)
    print(d_config)
    
    d_xsec = cmut.load_config(d_config["xsecfile_bkg"])
    #print(d_xsec)
    
    d_cutflows = cmut.load_config(d_config["cutflowsfile"])
    #print(d_cutflows)
    
    d_neventtot = {}
    
    samples = []
    
    for proctype in ["procs_sig", "procs_bkg"] :
        
        for proc in d_config[proctype] :
            
            samples.extend(d_config[proctype][proc])
            
            # Parse the signal sample name to the the stau mass
            # Get the corresponding xsec
            # Update the xsec dictionary
            if (proctype == "procs_sig") :
                
                d_xsec.update(cmut.get_stau_xsec_dict(l_samplestr = d_config[proctype][proc], xsecfile = d_config["xsecfile_sig"]))
    
    print(samples)
    
    for sample in samples :
        
        d_neventtot[sample] = functools.reduce(operator.getitem, [sample]+d_config["neventkey"].split("."), d_cutflows)
    
    print(d_neventtot)
    
    inhistfile = ROOT.TFile.Open(d_config["inhistfile"])
    os.system(f"mkdir -p {d_config['prephistdir']}")
    
    for proctype in ["procs_sig", "procs_bkg"] :
        
        for proc in d_config[proctype] :
            
            outfilename = f"{d_config['prephistdir']}/{proc}.root"
            outhistfile = ROOT.TFile.Open(outfilename, "RECREATE")
            
            hist_proc = None
            
            for sample in d_config[proctype][proc] :
                
                histname = f"{sample}{d_config['histnametag']}"
                hist_sample = inhistfile.Get(histname)
                hist_sample.Sumw2()
                
                scale = d_config["lumi"] * d_xsec[sample] / d_neventtot[sample]
                hist_sample.Scale(scale)
                
                # Fix negative weights
                nbins = hist_sample.GetNbinsX()
                for bin in range(1, nbins+1) :
                    
                    if (hist_sample.GetBinContent(bin) < 0) :
                        
                        hist_sample.SetBinContent(bin, 0.0)
                        hist_sample.SetBinError(bin, 0.0)
                
                if (hist_proc is None) :
                    
                    hist_proc = hist_sample.Clone()
                    hist_proc.SetDirectory(0)
                
                else :
                    
                    hist_proc.Add(hist_sample)
            
            hist_proc.SetName(proc)
            hist_proc.SetTitle(proc)
            
            outhistfile.cd()
            hist_proc.Write()
            outhistfile.Close()
    
    inhistfile.Close()


if (__name__ == "__main__") :
    
    main()