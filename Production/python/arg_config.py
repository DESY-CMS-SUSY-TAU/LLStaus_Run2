import os

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing


############################## Parse arguments ##############################

def get_args() :
    
    args = VarParsing.VarParsing("analysis")
    
    args.register("sourceFile",
        "", # Default value
        VarParsing.VarParsing.multiplicity.singleton, # singleton or list
        VarParsing.VarParsing.varType.string, # string, int, or float
        "File containing list of input files" # Description
    )
    
    args.register("outFile",
        "nanoaod.root", # Default value
        VarParsing.VarParsing.multiplicity.singleton, # singleton or list
        VarParsing.VarParsing.varType.string, # string, int, or float
        "Output file base name (w/o extension): [base name].root" # Description
    )
    
    args.register("eventRange",
        [], # Default value
        VarParsing.VarParsing.multiplicity.list, # singleton or list
        VarParsing.VarParsing.varType.string, # string, int, or float
        "Syntax: Run1:Event1-Run2:Event2 Run3:Event3-Run4:Event4(includes both)" # Description
    )
    
    args.register("isMC",
        1, # Default value
        VarParsing.VarParsing.multiplicity.singleton, # singleton or list
        VarParsing.VarParsing.varType.int, # string, int, or float
        "Whether Data or MC" # Description
    )
    
    args.register("debug",
        0, # Default value
        VarParsing.VarParsing.multiplicity.singleton, # singleton or list
        VarParsing.VarParsing.varType.int, # string, int, or float
        "Print debug statements" # Description
    )
    
    args.register("debugEDM",
        0, # Default value
        VarParsing.VarParsing.multiplicity.singleton, # singleton or list
        VarParsing.VarParsing.varType.int, # string, int, or float
        "Create EDM file for debugging collection content" # Description
    )
    
    args.parseArguments()
    
    if (len(args.sourceFile)) :
        
        sourceFile = args.sourceFile
    
    
    fNames = []
    
    if (len(args.inputFiles)) :
        
        fNames = args.inputFiles
    
    elif (len(args.sourceFile)) :
        
        with open(sourceFile) as f:
            
            args.inputFiles = f.readlines()
            args.inputFiles = [_fname.strip() for _fname in args.inputFiles]
            args.inputFiles = [_fname for _fname in args.inputFiles if _fname[0] != "#"]
    
    
    if ("/" in args.outFile) :
        
        outDir = args.outFile[0: args.outFile.rfind("/")]
        
        os.system("mkdir -p %s" %(outDir))
    
    
    return args
