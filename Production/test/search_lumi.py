import subprocess
from multiprocessing import Pool
from tqdm import tqdm
from functools import partial
import argparse

parser = argparse.ArgumentParser(description=\
'''
The following script is searching for the file with luminosity block.
Example of using:
> python ./search_lumi.py -lu 358813 -ds "/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM"
''')
parser.add_argument('-lu','--lumi', help='luminosity block (e.g: 358813)', required=True)
parser.add_argument('-ds','--dataset', help='dataset DAS name in which to search for lumi-block', required=True)
args = parser.parse_args()

def get_files(datset):
    result = subprocess.getoutput(f'dasgoclient --query="file dataset={datset}"')
    list_res = result.split()
    return list_res
    
def get_lumis(file, lumi=89019):
    result = subprocess.getoutput(f'dasgoclient --query="lumi file={file}"')
    check = result.find(str(lumi)+"\n")
    if check != -1:
        print("Found! > ", lumi)
        print(check)
        print("In file:")
        print(file)

if __name__ == "__main__":
    
    files = get_files(args.dataset)
    with Pool(processes=20) as pool:
        results = list(tqdm(pool.imap(partial(get_lumis, lumi=args.lumi), files), total=len(files)))
    print("Search is done!")
