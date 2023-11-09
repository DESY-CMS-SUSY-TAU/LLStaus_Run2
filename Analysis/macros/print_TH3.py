from argparse import ArgumentParser
import ROOT
ROOT.gROOT.SetBatch(True)

def print_TH3_with_borders(hists_l):
    # Get the dimensions of the histogram
    nx = hists_l[0].GetNbinsX()
    ny = hists_l[0].GetNbinsY()
    nz = hists_l[0].GetNbinsZ()

    # Loop over all bins, including the overflow bins
    for i in range(nx + 2):
        for j in range(ny + 2):
            for k in range(nz + 2):
                # Get the bin content
                # content = hist.GetBinContent(i, j, k)

                # Get the bin borders
                xLow = hists_l[0].GetXaxis().GetBinLowEdge(i)
                xUp = hists_l[0].GetXaxis().GetBinUpEdge(i)
                yLow = hists_l[0].GetYaxis().GetBinLowEdge(j)
                yUp = hists_l[0].GetYaxis().GetBinUpEdge(j)
                zLow = hists_l[0].GetZaxis().GetBinLowEdge(k)
                zUp = hists_l[0].GetZaxis().GetBinUpEdge(k)

                contents = [hist.GetBinContent(i, j, k) for hist in hists_l]
                # Print the bin information
                # print("Bin ({}, {}, {}):".format(i, j, k))
                print("pt: [{}, {}]".format(xLow, xUp), "dz: [{}, {}]".format(yLow, yUp), "dxy: [{}, {}]".format(zLow, zUp), "cont:", contents)
                # print()
                
                

# def print_TH3_histogram():
#     # Assuming you have already created and filled your TH3 histogram
#     # Replace 'hist' with the actual name of your histogram object
#     hist = ROOT.TH3D("hist", "3D Histogram", 10, 0, 10, 10, 0, 10, 10, 0, 10)
#     hist.Fill(5, 5, 5)  # Example of filling a bin with some content

#     # Print the histogram content with bin borders
#     print_TH3_with_borders(hist)

parser = ArgumentParser(
    description="The following script prints the TH3 histogram.")
parser.add_argument(
    "hist_files",nargs="+", help="Path to histogram file")
args = parser.parse_args()
# print(args.hist_files)
# Call the function to print the histogram
# print_TH3_histogram()
files = []
hists = []
for path in args.hist_files:
    files.append( ROOT.TFile.Open(path, 'read') )
    hists.append(files[-1].Get("fake_rate_jet_pt_jet_dz_jet_dxy"))

print_TH3_with_borders(hists)