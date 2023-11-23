import ROOT as r
import json

#infile = r.TFile.Open("/nfs/dust/cms/user/sobhatta/work/stopPairToTau/analysis/CMSSW_10_5_0/src/stopPair/resources/htt_scalefactors_legacy_2018.root")
infile = r.TFile.Open("/home/soham/nfs_dust/user/sobhatta/work/stopPairToTau/analysis/CMSSW_10_5_0/src/stopPair/resources/htt_scalefactors_legacy_2018.root")
wspace = infile.Get("w")
roohistfunc = wspace.function("zptmass_weight_nom")
hist = wspace.genobj("hist_zptmass_weight_nom")

# Assuming you have already created a RooHistFunc named 'roohistfunc'

# Extract histogram information
hist = roohistfunc.dataHist()
#bin_contents = [hist.getBinContent(i) for i in range(1, hist.numEntries() + 1)]
#bin_errors = [hist.getBinError(i) for i in range(1, hist.numEntries() + 1)]
#bin_edges = [hist.getHist().GetXaxis().GetBinUpEdge(i) for i in range(1, hist.numEntries() + 2)]

# Extract function parameters
params = roohistfunc.getParameters(r.RooArgSet())
param_names = [param.GetName() for param in params]
param_values = [param.getVal() for param in params]

# Get the function formula
#function_formula = roohistfunc.formulaString(r.RooArgSet())

# Create a JSON representation
json_data = {
    "histogram": {
        #"bin_contents": bin_contents,
        #"bin_errors": bin_errors,
        #"bin_edges": bin_edges,
    },
    "function": {
        "parameters": dict(zip(param_names, param_values)),
        #"formula": function_formula,
    }
}

# Serialize to JSON string
json_string = json.dumps(json_data, indent=2)

print(json_string)


