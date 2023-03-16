#!/usr/bin/env python3

import argparse
import json
import os
import yaml
import functools
import operator

import CombineHarvester.CombineTools.ch as ch
import CombineHarvester.CombineTools.systematics.SMLegacy as SMLegacySysts
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
    
    era = d_config["era"]
    channel = "tauhtauh"
    nbins = 24
    binstart = 2
    categories = [(_binnum, f"bin{_binnum}") for _binnum in range(binstart, nbins+binstart)]
    bin_ids = [_cat[0] for _cat in categories]
    bin_names = [_cat[1] for _cat in categories]
    procs_sig = list(d_config["procs_sig"].keys())
    procs_bkg = list(d_config["procs_bkg"].keys())
    
    print(categories)
    
    cb = ch.CombineHarvester()
    
    cb.AddProcesses(
        mass = procs_sig,
        analysis = ["llstau"],
        era = [era],
        channel = [channel],
        procs = ["sig"],
        bin = categories,
        signal = True,
    )
    
    cb.AddProcesses(
        mass = ["*"],
        analysis = ["llstau"],
        era = [era],
        channel = [channel],
        procs = procs_bkg,
        bin = categories,
        signal = False,
    )
    
    for proc in procs_sig :
        
        infilename = f"{d_config['prephistdir']}/{proc}.root"
        infile = ROOT.TFile.Open(infilename)
        histname = f"{proc}"
        hist = infile.Get(histname)
        
        cb.cp().mass([proc]).channel([channel]).process(["sig"]).era([era]).ForEachProc(
          lambda x : x.set_rate(
            hist.GetBinContent(x.bin_id()) if hist.GetBinContent(x.bin_id()) else 1e-3
        ))
        
        for bin_id, bin_name in categories :
            
            val = hist.GetBinContent(bin_id)
            err = hist.GetBinError(bin_id)
            err_rel = 1.0
            
            if (val) :
                
                err_rel = 1.0+(err/val)
            
            else :
                
                err_rel = 2.0
            
            cb.cp().mass([proc]).bin([bin_name]).process(["sig"]).AddSyst(
                target = cb,
                name = "stat_$PROCESS_$BIN",
                type = "lnN",
                valmap = ch.SystMap()
                    (err_rel)
            )
        
        infile.Close()
    
    for proc in procs_bkg :
        
        infilename = f"{d_config['prephistdir']}/{proc}.root"
        infile = ROOT.TFile.Open(infilename)
        histname = f"{proc}"
        hist = infile.Get(histname)
        
        cb.cp().channel([channel]).process([proc]).era([era]).ForEachProc(
          lambda x : x.set_rate(
            hist.GetBinContent(x.bin_id()) if hist.GetBinContent(x.bin_id()) else 1e-3
        ))
        
        for bin_id, bin_name in categories :
            
            val = hist.GetBinContent(bin_id)
            err = hist.GetBinError(bin_id)
            err_rel = 1.0
            
            if (val) :
                
                err_rel = 1.0+(err/val)
            
            else :
                
                err_rel = 2.0
            
            cb.cp().bin([bin_name]).process([proc]).AddSyst(
                target = cb,
                name = "stat_$PROCESS_$BIN",
                type = "lnN",
                valmap = ch.SystMap()
                    (err_rel)
            )
        
        infile.Close()
    
    cb.cp().AddSyst(
        target = cb,
        name = "dummysyst",
        type = "lnN",
        valmap = ch.SystMap("channel", "era", "bin_id")
            ([channel], [era], bin_ids, 1.3)
    )
    
    #cb.cp().bin({"A"}).AddSyst(cb, "scale_$BIN", "rateParam", SystMapFunc<>::init
    #    ("(@0*@1/@2)", "scale_B,scale_C,scale_D")
    #);
    
    #cb.cp().bin(bin_names).AddSyst(
    #    target = cb,
    #    name = "stat_$BIN",
    #    type = "lnN",
    #    valmap = ch.SystMap("channel", "era")
    #        ([channel], [era], 1.123)
    #)
    
    #writer = ch.CardWriter(
    #    '$TAG/$MASS/$ANALYSIS_$CHANNEL_$BINID_$ERA.txt',
    #    '$TAG/common/$ANALYSIS_$CHANNEL.input_$ERA.root'
    #)
    
    writer = ch.CardWriter(
        '$TAG/$MASS/$ANALYSIS_$CHANNEL_$ERA.txt',
        '$TAG/common/$ANALYSIS_$CHANNEL.input_$ERA.root'
    )
    
    writer.WriteCards(d_config["carddir"], cb)


if (__name__ == "__main__") :
    
    main()