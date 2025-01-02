import os
import json
import uproot
import numpy as np
import coffea.util
from coffea import hist as hi
from itertools import product
import sys

class HistCollection:
    def __init__(self, path):
        self.path = path
        self._content = self._scan_histograms()

    def _scan_histograms(self):
        # Auto-detect files and parse keys based on filenames
        content = {}
        for filename in os.listdir(self.path):
            if filename.endswith('.coffea'):
                # Assuming filename structure: key1_key2_key3.coffea
                key = tuple(filename[:-7].split('_'))
                content[key] = filename
        return content

    def load(self, key):
        path = os.path.join(self.path, self._content[key])
        if path.endswith(".root"):
            with uproot.open(path) as f:
                return self.rootdir_to_hist(f)
        else:
            return coffea.util.load(path)

    def rootdir_to_hist(self, rootdir, infokey="hist_info"):
        # Assume all necessary information to rebuild hist from ROOT is stored in infokey
        info = json.loads(rootdir[infokey].member("fTitle"))
        axes = [hi.axis.Regular(**ax) for ax in info["axes"]]
        hist = hi.Hist(*axes, storage=hi.storage.Weight())
        # Assume the directory's structure and naming follows a particular pattern
        for k, v in rootdir.items():
            hist.fill(**{ax.name: v.member("fXaxis").edges[:-1] for ax in axes}, weight=v.values())
        return hist

    def save_all_to_root(self, root_filename):
        with uproot.recreate(root_filename) as f:
            for key, _ in self._content.items():
                hist = self.load(key)
                # Convert hist to a format suitable for ROOT and save
                data = hist.values(flow=True)
                edges = [ax.edges for ax in hist.axes]
                root_key = self.root_key(key)
                f[root_key] = (data, edges)

    @staticmethod
    def root_key(key):
        return "/".join(map(str, key)) + "/hist"

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <hist_directory> <output_root_file>")
        sys.exit(1)

    hist_directory = sys.argv[1]
    output_root_file = sys.argv[2]

    hist_collection = HistCollection(hist_directory)
    hist_collection.save_all_to_root(output_root_file)

    print(f"All histograms have been successfully saved to {output_root_file}")
