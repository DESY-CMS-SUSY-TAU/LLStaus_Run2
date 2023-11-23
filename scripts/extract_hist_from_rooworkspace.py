import argparse
import ROOT
import os



def main() :
    
    # Argument parser
    parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument(
        "--infilename",
        help = "Input ROOT file",
        type = str,
        required = True,
    )
    
    parser.add_argument(
        "--wspace",
        help = "Workspace name in the file",
        type = str,
        required = True,
    )
    
    parser.add_argument(
        "--histname",
        help = "Histogram name in the workspace",
        type = str,
        required = True,
    )
    
    parser.add_argument(
        "--outfilename",
        help = "Output file name",
        type = str,
        required = False,
        default = None
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    #infile = ROOT.TFile.Open("/nfs/dust/cms/user/sobhatta/work/stopPairToTau/analysis/CMSSW_10_5_0/src/stopPair/resources/htt_scalefactors_legacy_2018.root")
    #infile = ROOT.TFile.Open("/home/soham/nfs_dust/user/sobhatta/work/stopPairToTau/analysis/CMSSW_10_5_0/src/stopPair/resources/htt_scalefactors_legacy_2018.root")
    
    infile = ROOT.TFile.Open(args.infilename)
    
    wspace = infile.Get(args.wspace)
    #hist = wspace.genobj("hist_zptmass_weight_nom")
    hist = wspace.genobj(args.histname).Clone()
    hist.SetDirectory(0)
    
    infile.Close()
    
    if (args.outfilename is None) :
        args.outfilename = f"{args.histname}.root"
    
    outdir = os.path.dirname(args.outfilename)
    
    if len(outdir) :
        
        os.system(f"mkdir -p {outdir}")
    
    outfile = ROOT.TFile.Open(args.outfilename, "RECREATE")
    outfile.cd()
    
    hist.Write()
    
    print("Done")


if (__name__ == "__main__") :
    
    main()