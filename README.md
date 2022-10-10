# Long-lived stau repository

This repository is to perform the long-lives stau searches at DESY CMS.

# Environment installation
1. Clone package from the github without loading any additional environment (like CMSSW):
    ```sh
    > git clone -o DESY-CMS-SUSY-TAU -b master --recurse-submodules git@github.com:DESY-CMS-SUSY-TAU/LLStaus_Run2.git
    > cd LLStaus_Run2
    > git submodule update --init --recursive # Update the submodules
    ```
    Update git to >=2.13 if possible, for `--recurse-submodules` to work.<br>
    Otherwise, use `--recursive` instead of `--recurse-submodules` with 1.9<=git<=2.12.<br>
    Note that this is needed to clone the submodules (such as `pepper`).
2. Go to the directory and setup appropriate environment with setup.sh script:
    ```sh
    > source setup.sh ENVIRONMENT
    ```
    where supported `ENVIRONMENT` are:
    - `NanoAOD_UL2018`: step for production of the customized NanoAOD tuples. 
    - `conda`: anaconda environment - used for all the analysis except the production of NanoAOD.

# Nano Tuple production
See documentation [here](Analysis/README.md)

# Analysis
See documentation [here](Production/README.md)
