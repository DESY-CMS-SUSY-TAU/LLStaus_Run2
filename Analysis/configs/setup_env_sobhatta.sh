#!/bin/bash
PEPPER_PATH="/nfs/dust/cms/user/sobhatta/work/LongLivedStaus/LLStaus_Run2/pepper"

cd /nfs/dust/cms/user/sobhatta/work/LongLivedStaus/LLStaus_Run2/Analysis
source /afs/desy.de/user/s/sobhatta/.local/bin/init-python3p8.sh

if [[ ${PYTHONPATH} != *"$PEPPER_PATH"* ]]; then
    echo "Adding $PEPPER_PATH to PYTHONPATH."
    export PYTHONPATH=$PYTHONPATH:$PEPPER_PATH
fi


