import os

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing


############################## Parse arguments ##############################

def get_args() :
    
    args = VarParsing("analysis")
    
    args.register("sourceFile",
        "", # Default value
        VarParsing.multiplicity.singleton, # singleton or list
        VarParsing.varType.string, # string, int, or float
        "File containing list of input files" # Description
    )

    # args.register('fileList',
    #     '[]',
    #     VarParsing.multiplicity.list,
    #     VarParsing.varType.string,
    #     "List of root files to process (alternative to sourceFile option).")

    args.register("outFile",
        "nanoaod.root", # Default value
        VarParsing.multiplicity.singleton, # singleton or list
        VarParsing.varType.string, # string, int, or float
        "Output file base name (w/o extension): [base name].root" # Description
    )

    args.register('fileNamePrefix',
        '',
        VarParsing.multiplicity.singleton,
        VarParsing.varType.string,
        "Prefix to add to input file names.")

    args.register('lumiFile',
        '', 
        VarParsing.multiplicity.singleton, 
        VarParsing.varType.string,
        "JSON file with lumi mask.")
    
    args.register("eventRange",
        '', # Default value
        VarParsing.multiplicity.singleton, # singleton
        VarParsing.varType.string, # string, int, or float
        "Syntax: Run1:Event1-Run2:Event2 Run3:Event3-Run4:Event4(includes both)" # Description
    )
    
    args.register("isMC",
        True, # Default value
        VarParsing.multiplicity.singleton, # singleton or list
        VarParsing.varType.bool, # string, int, or float
        "Whether Data or MC" # Description
    )
    
    # args.register("debug",
    #     0, # Default value
    #     VarParsing.VarParsing.multiplicity.singleton, # singleton or list
    #     VarParsing.VarParsing.varType.int, # string, int, or float
    #     "Print debug statements" # Description
    # )
    
    args.register("debugEDM",
        False, # Default value
        VarParsing.multiplicity.singleton, # singleton or list
        VarParsing.varType.bool, # string, int, or float
        "Create EDM file for debugging collection content" # Description
    )
    
    args.parseArguments()
    
    fNames = []
    
    if(len(args.inputFiles) and len(args.sourceFile)):
        raise ValueError(
            'Error: fileList and sourceFile are interchangeable, \
             only one should be specified.'
            )

    if (len(args.inputFiles)) :
        fNames = args.inputFiles
    
    elif (len(args.sourceFile)) :
        with open(args.sourceFile) as f:
            args.fileList = f.readlines()
            args.fileList = [_fname.strip() for _fname in args.inputFiles]
            args.fileList = [_fname for _fname in args.inputFiles if _fname[0] != "#"]
    
    if ("/" in args.outFile) :
        outDir = args.outFile[0: args.outFile.rfind("/")]
        os.system("mkdir -p %s" %(outDir))
    
    return args
