import numpy as np
import ROOT

class TH3Histogram:
    def __init__(self, existing_hist, new_bin_edges_x, new_bin_edges_y, new_bin_edges_z):
        self.bin_edges_x = np.array(new_bin_edges_x, dtype=np.float32)
        self.bin_edges_y = np.array(new_bin_edges_y, dtype=np.float32)
        self.bin_edges_z = np.array(new_bin_edges_z, dtype=np.float32)
        self.bins = self._rebin(existing_hist, self.bin_edges_x, self.bin_edges_y, self.bin_edges_z)
    
    def _rebin(self, existing_hist, new_bin_edges_x, new_bin_edges_y, new_bin_edges_z):
        existing_bin_edges_x = np.array([existing_hist.GetXaxis().GetBinLowEdge(bin) for bin in range(1, existing_hist.GetNbinsX() + 2)], dtype=np.float32)
        existing_bin_edges_y = np.array([existing_hist.GetYaxis().GetBinLowEdge(bin) for bin in range(1, existing_hist.GetNbinsY() + 2)], dtype=np.float32)
        existing_bin_edges_z = np.array([existing_hist.GetZaxis().GetBinLowEdge(bin) for bin in range(1, existing_hist.GetNbinsZ() + 2)], dtype=np.float32)
        
        if not (np.isin(new_bin_edges_x, existing_bin_edges_x).all()):
            print("OLD bins:", existing_bin_edges_x)
            print("NEW bins:", new_bin_edges_x)
            print(np.isin(new_bin_edges_x, existing_bin_edges_x))
            raise ValueError("The new bin edges for X-axis are not equal to one of the existing histogram's bin edges.")
        
        if not (np.isin(new_bin_edges_y, existing_bin_edges_y).all()):
            print("OLD bins:", existing_bin_edges_y)
            print("NEW bins:", new_bin_edges_y)
            print(np.isin(new_bin_edges_y, existing_bin_edges_y))
            raise ValueError("The new bin edges for Y-axis are not equal to one of the existing histogram's bin edges.")
        
        if not (np.isin(new_bin_edges_z, existing_bin_edges_z).all()):
            print("OLD bins:", existing_bin_edges_z)
            print("NEW bins:", new_bin_edges_z)
            print(np.isin(new_bin_edges_z, existing_bin_edges_z))
            raise ValueError("The new bin edges for Z-axis are not equal to one of the existing histogram's bin edges.")
        
        
        n_bins_x = len(new_bin_edges_x) - 1
        n_bins_y = len(new_bin_edges_y) - 1
        n_bins_z = len(new_bin_edges_z) - 1
        
        new_bins = np.zeros((n_bins_x, n_bins_y, n_bins_z))
        
        for i in range(n_bins_x):
            for j in range(n_bins_y):
                for k in range(n_bins_z):
                    x_indices = range(existing_hist.GetXaxis().FindBin(new_bin_edges_x[i]), existing_hist.GetXaxis().FindBin(new_bin_edges_x[i + 1]))
                    y_indices = range(existing_hist.GetYaxis().FindBin(new_bin_edges_y[j]), existing_hist.GetYaxis().FindBin(new_bin_edges_y[j + 1]))
                    z_indices = range(existing_hist.GetZaxis().FindBin(new_bin_edges_z[k]), existing_hist.GetZaxis().FindBin(new_bin_edges_z[k + 1]))
                    
                    # Sum up the bin contents from the existing bins
                    new_bins[i, j, k] = np.sum(existing_hist.GetBinContent(x, y, z) for x in x_indices for y in y_indices for z in z_indices)
        
        return new_bins
    
    def get_th3_histogram(self):
        n_bins_x = len(self.bin_edges_x) - 1
        n_bins_y = len(self.bin_edges_y) - 1
        n_bins_z = len(self.bin_edges_z) - 1
        
        th3_hist = ROOT.TH3D("hist", "Rebinned TH3 Histogram", n_bins_x, self.bin_edges_x, n_bins_y, self.bin_edges_y, n_bins_z, self.bin_edges_z)
        
        for i in range(n_bins_x):
            for j in range(n_bins_y):
                for k in range(n_bins_z):
                    th3_hist.SetBinContent(i + 1, j + 1, k + 1, self.bins[i, j, k])
        
        return th3_hist
    
if __name__ == "__main__":
    
    # Define the existing TH3D histogram
    existing_hist = ROOT.TH3D("existing_hist", "Existing TH3 Histogram", 10, 0, 10, 10, 0, 10, 10, 0, 10)

    # Fill the existing histogram with some values
    for i in range(1, existing_hist.GetNbinsX() + 1):
        for j in range(1, existing_hist.GetNbinsY() + 1):
            for k in range(1, existing_hist.GetNbinsZ() + 1):
                existing_hist.SetBinContent(i, j, k, i + j + k)
    
    # Define the new bin edges for rebinning
    new_bin_edges_x = [0, 2, 4, 6, 10]
    new_bin_edges_y = [0, 3, 6, 9, 10]
    new_bin_edges_z = [0, 2, 4, 6, 8, 10]

    # Create an instance of the TH3Histogram class
    rebinned_hist = TH3Histogram(existing_hist, new_bin_edges_x, new_bin_edges_y, new_bin_edges_z)

    # Get the rebinned TH3 histogram
    rebinned_th3_hist = rebinned_hist.get_th3_histogram()

    # Verify the rebinned histogram properties
    assert rebinned_th3_hist.GetNbinsX() == len(new_bin_edges_x) - 1
    assert rebinned_th3_hist.GetNbinsY() == len(new_bin_edges_y) - 1
    assert rebinned_th3_hist.GetNbinsZ() == len(new_bin_edges_z) - 1

    # Verify the bin contents in the rebinned histogram
    for i in range(1, rebinned_th3_hist.GetNbinsX() + 1):
        for j in range(1, rebinned_th3_hist.GetNbinsY() + 1):
            for k in range(1, rebinned_th3_hist.GetNbinsZ() + 1):
                bin_content_rebinned = rebinned_th3_hist.GetBinContent(i, j, k)
                
                # Calculate the expected bin content based on the rebinning
                expected_bin_content = np.sum(
                    existing_hist.GetBinContent(x, y, z)
                    for x in range(existing_hist.GetXaxis().FindBin(new_bin_edges_x[i - 1]), existing_hist.GetXaxis().FindBin(new_bin_edges_x[i]))
                    for y in range(existing_hist.GetYaxis().FindBin(new_bin_edges_y[j - 1]), existing_hist.GetYaxis().FindBin(new_bin_edges_y[j]))
                    for z in range(existing_hist.GetZaxis().FindBin(new_bin_edges_z[k - 1]), existing_hist.GetZaxis().FindBin(new_bin_edges_z[k]))
                )
                
                assert bin_content_rebinned == expected_bin_content
    
    print("X projection of rebinned histogram:")
    print(rebinned_th3_hist.Project3D("x").Print('all'))
    
    print("X projection of initial histogram:")
    print(existing_hist.Project3D("x").Print('all'))
    
    print("Test passed successfully!")