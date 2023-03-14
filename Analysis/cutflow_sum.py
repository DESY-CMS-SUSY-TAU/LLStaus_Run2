import os
import json
from argparse import ArgumentParser
import re

parser = ArgumentParser(
    description="Cut flow sum of signal hists")
parser.add_argument(
    "xsec", help="Crosssection file path")


args = parser.parse_args()

try:
    with open(args.xsec) as f:
        crosssections = json.load(f)
except:
    raise ValueError('Error reading/open file with crosssections')

for data in crosssections.keys():
    print(data, crosssections[data]["all"]["NanDrop"])
    print("categories: ", "(PP:", crosssections[data]["PP"]["S_INCL"]["b_tagged_0_cut"],
                          "| PDtag:",crosssections[data]["PDtag"]["S_INCL"]["b_tagged_0_cut"],
                          "| DDtag:", crosssections[data]["DDtag"]["S_INCL"]["b_tagged_0_cut"],")",
                          "sum:",crosssections[data]["PP"]["S_INCL"]["b_tagged_0_cut"]+crosssections[data]["PDtag"]["S_INCL"]["b_tagged_0_cut"]+crosssections[data]["DDtag"]["S_INCL"]["b_tagged_0_cut"]
                          )
    print(crosssections[data]["PP"]["S_INCL"]["b_tagged_0_cut"] / crosssections[data]["all"]["NanDrop"])
#     print(data)