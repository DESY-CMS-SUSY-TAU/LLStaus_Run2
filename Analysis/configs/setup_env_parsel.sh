#!/bin/bash

eval "$(/nfs/dust/cms/user/mykytaua/softML/miniconda3_v2/bin/conda shell.bash hook)"
conda activate llstau

# export mamba_nfs_prefix=/nfs/dust/cms/user/mykytaua/softML/mamba_v0/
# export MAMBA_ROOT_PREFIX=${mamba_nfs_prefix}  # optional, defaults to ~/micromamba
# eval "$(${mamba_nfs_prefix}/bin/micromamba shell hook -s posix)"
# micromamba activate llstau-new