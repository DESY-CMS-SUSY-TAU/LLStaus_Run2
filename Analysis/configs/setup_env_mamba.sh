#!/bin/bash

export mamba_nfs_prefix=/nfs/dust/cms/user/mykytaua/softML/mamba_v0/
export MAMBA_ROOT_PREFIX=${mamba_nfs_prefix}  # optional, defaults to ~/micromamba
eval "$(${mamba_nfs_prefix}/bin/micromamba shell hook -s posix)"
micromamba activate llstau-new