import numpy as np
import ROOT
import unittest

class TH3Histogram:
    def __init__(self, hist, new_bin_edges_x, new_bin_edges_y, new_bin_edges_z):
        self._hist = hist
        self._new_bin_edges_x = np.array(new_bin_edges_x, dtype=np.float32)
        self._new_bin_edges_y = np.array(new_bin_edges_y, dtype=np.float32)
        self._new_bin_edges_z = np.array(new_bin_edges_z, dtype=np.float32)
        self._rebin()

    def _rebin(self):

        old_bin_edges_x = np.array([self._hist.GetXaxis().GetBinLowEdge(bin) for bin in range(1, self._hist.GetNbinsX() + 2)], dtype=np.float32)
        old_bin_edges_y = np.array([self._hist.GetYaxis().GetBinLowEdge(bin) for bin in range(1, self._hist.GetNbinsY() + 2)], dtype=np.float32)
        old_bin_edges_z = np.array([self._hist.GetZaxis().GetBinLowEdge(bin) for bin in range(1, self._hist.GetNbinsZ() + 2)], dtype=np.float32)
        
        if not (np.isin(self._new_bin_edges_x, old_bin_edges_x).all()):
            print("OLD bins:", old_bin_edges_x)
            print("NEW bins:", self._new_bin_edges_x)
            print(np.isin(self._new_bin_edges_x, old_bin_edges_x))
            raise ValueError("The new bin edges for X-axis are not equal to one of the existing histogram's bin edges.")
        
        if not (np.isin(self._new_bin_edges_y, old_bin_edges_y).all()):
            print("OLD bins:", old_bin_edges_y)
            print("NEW bins:", self._new_bin_edges_y)
            print(np.isin(self._new_bin_edges_y, old_bin_edges_y))
            raise ValueError("The new bin edges for Y-axis are not equal to one of the existing histogram's bin edges.")
        
        if not (np.isin(self._new_bin_edges_z, old_bin_edges_z).all()):
            print("OLD bins:", old_bin_edges_z)
            print("NEW bins:", self._new_bin_edges_z)
            print(np.isin(self._new_bin_edges_z, old_bin_edges_z))
            raise ValueError("The new bin edges for Z-axis are not equal to one of the existing histogram's bin edges.")
        

        # Create the rebinned histogram
        rebinned_hist = ROOT.TH3D(
            f"{self._hist.GetName()}_rebinned",
            f"Rebinned TH3D Histogram",
            len(self._new_bin_edges_x) - 1,
            self._new_bin_edges_x,
            len(self._new_bin_edges_y) - 1,
            self._new_bin_edges_y,
            len(self._new_bin_edges_z) - 1,
            self._new_bin_edges_z
        )

        # Fill the rebinned histogram
        for i in range(0, self._hist.GetNbinsX() + 2):
            for j in range(0, self._hist.GetNbinsY() + 2):
                for k in range(0, self._hist.GetNbinsZ() + 2):
                    content = self._hist.GetBinContent(i, j, k)
                    bin_x = self._hist.GetXaxis().GetBinCenter(i)
                    bin_y = self._hist.GetYaxis().GetBinCenter(j)
                    bin_z = self._hist.GetZaxis().GetBinCenter(k)
                    rebinned_hist.Fill(bin_x, bin_y, bin_z, content)

        self._rebinned_hist = rebinned_hist

    def get_rebinned_histogram(self):
        return self._rebinned_hist
    
class TH3HistogramTest(unittest.TestCase):
    def test_rebinning(self):
        # Create a TH3D histogram for testing
        hist = ROOT.TH3D("hist", "Original TH3D Histogram", 10, 0, 10, 10, 0, 10, 10, 0, 10)
        
        # Fill the histogram with some data
        for i in range(100):
            x = ROOT.gRandom.Uniform(0, 10)
            y = ROOT.gRandom.Uniform(0, 10)
            z = ROOT.gRandom.Uniform(0, 10)
            hist.Fill(x, y, z)
        
        # [print(hist.GetXaxis().GetBinLowEdge(bin)) for bin in range(1, hist.GetNbinsX() + 2)]
        
        # Fill the underflow bin
        hist.Fill(-1, -1, -1)
        
        # Fill the overflow bin
        hist.Fill(11, 11, 11)
        
        # Define the new bin edges for rebinning
        new_bin_edges_x = [0, 2, 4, 6, 8, 9]
        new_bin_edges_y = [0, 2, 4, 6, 8, 10]
        new_bin_edges_z = [0, 2, 4, 6, 8, 10]
        
        # Rebin the histogram
        rebinned_hist = TH3Histogram(hist, new_bin_edges_x, new_bin_edges_y, new_bin_edges_z)
        
        # Get the rebinned TH3D histogram
        rebinned_th3_hist = rebinned_hist.get_rebinned_histogram()
        
        # Check the dimensions of the rebinned histogram
        self.assertEqual(rebinned_th3_hist.GetNbinsX(), len(new_bin_edges_x) - 1)
        self.assertEqual(rebinned_th3_hist.GetNbinsY(), len(new_bin_edges_y) - 1)
        self.assertEqual(rebinned_th3_hist.GetNbinsZ(), len(new_bin_edges_z) - 1)
        

        projection_before = hist.Project3D("x")
        projection_after = rebinned_th3_hist.Project3D("x")
        
        print("X projection of initial histogram:")
        print(projection_before.Print('all'))

        print("X projection of rebinned histogram:")
        print(projection_after.Print('all'))
        
        self.assertAlmostEqual(projection_before.Integral() + projection_before.GetBinContent(0) +  projection_before.GetBinContent(projection_before.GetNbinsX()+1),
                               projection_after.Integral() + projection_after.GetBinContent(0) +  projection_after.GetBinContent(projection_after.GetNbinsX()+1)
                               )
        
        

if __name__ == '__main__':
    
    unittest.main()