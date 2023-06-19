#!/usr/bin/env python3

import argparse
import glob
import re
import uproot

#import ROOT


def natural_sort(l):
    
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


def main() :
    
    # Argument parser
    parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument(
        "--path",
        help = "Glob path. For e.g. a/b/*/*/.root",
        type = str,
        required = True,
    )
    
    parser.add_argument(
        "--tree",
        help = "Tree path in the file",
        type = str,
        default = "Events",
        required = False,
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    l_filename = natural_sort(glob.glob(args.path))
    nfile = len(l_filename)
    
    nevent_total = 0
    l_skipped = []
    
    print(f"Reading {nfile} files from: {args.path}")
    
    for ifile, filename in enumerate(l_filename) :
        
        try :
            with uproot.open(f"{filename}:{args.tree}") as tree :
                
                nevent = tree.num_entries
                print(f"    Read file [{ifile+1}/{nfile}] [{filename}] [{nevent} entries]")
                nevent_total += nevent
        
        except Exception as exc :
            
            print(f"    Skipped file [{ifile+1}/{nfile}] [{filename}]")
            l_skipped.append(filename)
    
    print(f"Path: {args.path}")
    print(f"Tree: {args.tree}")
    print(f"Total entries: {nevent_total}")
    print(f"Skipped: {len(l_skipped)} files")
    print("\n".join(l_skipped)+"\n")
    
    return 0

if __name__ == "__main__" :
    
    main()