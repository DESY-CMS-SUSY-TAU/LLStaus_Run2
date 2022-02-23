import argparse
import os


# Argument parser
parser = argparse.ArgumentParser(formatter_class = argparse.RawTextHelpFormatter)


parser.add_argument(
    "--getCount",
    help = "Get the number of events",
    default = False,
    action = "store_true",
)

parser.add_argument(
    "--getFiles",
    help = "Get the list of files",
    default = False,
    action = "store_true",
)

parser.add_argument(
    "--replace",
    help = "Will replace \"str1\" with \"str2\": \"str1\" \"str2\" (e.g. \"/store\" \"root://dcache-cms-xrootd.desy.de://pnfs/desy.de/cms/tier2/store\")",
    type = str,
    nargs = 2,
    required = False,
    default = ["/store", "root://cms-xrd-global.cern.ch//store"],
)

parser.add_argument(
    "--instance",
    help = "Choose the DBS instance (such as \"prod/phys03\" for USER produced datasets)",
    type = str,
    default = "",
    choices = ["prod/phys03"],
)

parser.add_argument(
    "--createRucioRule",
    help = "Create Rucio rules to copy datasets",
    default = False,
    action = "store_true",
)

parser.add_argument(
    "--tier2site",
    help = "T2 site name for Rucio rules (e.g. T2_DE_DESY)",
    type = str,
    required = False,
    default = None,
)

parser.add_argument(
    "--sampleNames",
    help = "Sample names (if provided, will not use the internal list)",
    type = str,
    required = False,
    nargs = "*",
    default = None,
)



# Parse arguments
args = parser.parse_args()
d_args = vars(args)


outDir = "sourceFiles"

prefix = args.replace[1]
toReplace = args.replace[0]


l_sampleName = [
    
    #"/ZprimeToTT_M1000_W10_TuneCP2_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    #"/ZprimeToTT_M1500_W15_TuneCP2_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    #"/ZprimeToTT_M2000_W20_TuneCP2_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    
    #"/TT_Mtt-700to1000_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v8/MINIAODSIM",
    #"/TT_Mtt-1000toInf_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM",
    
    #"/TT_Mtt-700to1000_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    #"/TT_Mtt-1000toInf_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
    
    #"/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    #"/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    
    #"/QCD_bEnriched_HT100to200_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    #"/QCD_bEnriched_HT200to300_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
    #"/QCD_bEnriched_HT300to500_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
    #"/QCD_bEnriched_HT500to700_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
    #"/QCD_bEnriched_HT700to1000_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
    #"/QCD_bEnriched_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
    #"/QCD_bEnriched_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
    #"/QCD_bEnriched_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",

    "/ZprimeToTT_M1000_W10_TuneCP2_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    "/ZprimeToTT_M1000_W100_TuneCP2_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    "/ZprimeToTT_M1000_W300_TuneCP2_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    "/ZprimeToTT_M1250_W12p5_TuneCP2_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    "/ZprimeToTT_M1250_W125_TuneCP2_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    "/ZprimeToTT_M1250_W375_TuneCP2_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    "/ZprimeToTT_M1500_W15_TuneCP2_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    "/ZprimeToTT_M1500_W150_TuneCP2_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    "/ZprimeToTT_M1500_W450_TuneCP2_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    "/ZprimeToTT_M2000_W20_TuneCP2_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    "/ZprimeToTT_M2000_W200_TuneCP2_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    "/ZprimeToTT_M2000_W600_TuneCP2_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    "/ZprimeToTT_M3000_W30_TuneCP2_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    "/ZprimeToTT_M3000_W300_TuneCP2_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    "/ZprimeToTT_M3000_W900_TuneCP2_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    "/ZprimeToTT_M3500_W35_TuneCP2_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    "/ZprimeToTT_M3500_W350_TuneCP2_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    "/ZprimeToTT_M3500_W1050_TuneCP2_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    "/ZprimeToTT_M4000_W40_TuneCP2_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    "/ZprimeToTT_M4000_W400_TuneCP2_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    "/ZprimeToTT_M4000_W1200_TuneCP2_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    
    #"/QCD_bEnriched_HT500to700_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
    #"/QCD_bEnriched_HT700to1000_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
    #"/QCD_bEnriched_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
    #"/QCD_bEnriched_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
    #"/QCD_bEnriched_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
    
    ##"/QCD_Pt_15to30_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-PREMIX_RECODEBUG_102X_upgrade2018_realistic_v15_ext2-v2/MINIAODSIM",
    ##"/QCD_Pt_30to50_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-PREMIX_RECODEBUG_102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    ##"/QCD_Pt_50to80_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-PREMIX_RECODEBUG_102X_upgrade2018_realistic_v15_ext2-v1/MINIAODSIM",
    ##"/QCD_Pt_80to120_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-PREMIX_RECODEBUG_102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    ##"/QCD_Pt_120to170_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-PREMIX_RECODEBUG_102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    "/QCD_Pt_170to300_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-PREMIX_RECODEBUG_102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    "/QCD_Pt_300to470_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-PREMIX_RECODEBUG_102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    "/QCD_Pt_470to600_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-PREMIX_RECODEBUG_102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM",
    "/QCD_Pt_600to800_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-PREMIX_RECODEBUG_102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    "/QCD_Pt_800to1000_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-PREMIX_RECODEBUG_102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM",
    "/QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-PREMIX_RECODEBUG_102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    "/QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-PREMIX_RECODEBUG_102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM",
    "/QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-PREMIX_RECODEBUG_102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM",
    "/QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8/RunIIAutumn18MiniAOD-PREMIX_RECODEBUG_102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM",

]


if (args.sampleNames is not None) :
    
    l_sampleName = args.sampleNames


# Reconfirm rucio rule creation request with user input
rucioConfirmed = False

if (args.createRucioRule and args.tier2site is not None) :
    
    print("Samples:")
    print("\n".join(l_sampleName))
    print("")
    
    print("REALLY create Rucio rules?")
    inputStr = str(raw_input("Enter CONFIRM to confirm: ")).strip()
    
    rucioConfirmed = (inputStr == "CONFIRM")
    
    if (not rucioConfirmed) :
        
        print("Rucio rule creation not confirmed. Exiting...")
        exit()


for iSample, sampleName in enumerate(l_sampleName) :
    
    print "\n"
    print "*"*50
    print "Sample %d/%d: %s" %(iSample+1, len(l_sampleName), sampleName)
    print "*"*50
    
    instance_str = "instance=%s" %(args.instance) if len(args.instance) else ""
    
    if (args.getCount) :
        
        command = "dasgoclient -query=\"file dataset=%s %s | sum(file.nevents)\"" %(sampleName, instance_str)
        os.system(command)
    
    
    if (args.getFiles) :
        
        sampleName_mod = sampleName[1:].replace("/", "_")
        
        outDir_mod = "%s/%s" %(outDir, sampleName_mod)
        
        command = "mkdir -p %s" %(outDir_mod)
        print "Command:", command
        print ""
        os.system(command)
        
        outFile = "%s/%s.txt" %(outDir_mod, sampleName_mod)
        
        command = "dasgoclient -query=\"file dataset=%s %s\" > %s" %(sampleName, instance_str, outFile)
        print "Command:", command
        print ""
        os.system(command)
        
        fileContent = ""
        
        print "Replacing \"%s\" with \"%s\" in file." %(toReplace, prefix)
        print ""
        
        print "Number of lines:"
        os.system("wc -l %s" %(outFile))
        print ""
        
        with open(outFile, "r") as f :
            
            fileContent = f.read()
        
        fileContent = fileContent.replace(toReplace, prefix)
        
        with open(outFile, "w") as f :
            
            f.write(fileContent)
    
    # https://twiki.cern.ch/twiki/bin/view/CMS/Rucio
    # https://twiki.cern.ch/twiki/bin/view/CMSPublic/RucioUserDocsQuotas
    # https://twiki.cern.ch/twiki/bin/view/CMSPublic/RucioUserDocsRules
    if (rucioConfirmed) :
        
        #command = "rucio add-rule cms:%s 1 %s" %(sampleName, args.tier2site)
        command = "rucio add-rule --ask-approval cms:%s 1 %s" %(sampleName, args.tier2site)
        print "Command:", command
        print ""
        os.system(command)
