#!/bin/bash

if [ $# -ne 1 ] ; then
    echo "To run: source setup.sh mode"
    echo "e.g: source setup.sh NanoAOD_UL2018"
    echo "Supported modes:"
    echo "1) CMSSW_12_4_0 (for NanoAOD UL 2018 production)"
    echo "2) conda (for analysis)"
    return 1
fi

MODE=$1
BASE_PATH=$PWD

function run_cmd {
    "$@"
    RESULT=$?
    if [ $RESULT -ne 0 ] ; then
        echo "Error while running '$@'"
        kill -INT $$
    fi
}

if [[ $MODE = "NanoAOD_UL2018" || "$MODE" == *"CMSSW"* ]] ; then

    if [[ "$MODE" == *"CMSSW"* ]] ; then
        CMSSW_VER=$MODE
        export SCRAM_ARCH=slc7_amd64_gcc10
    else
        echo "Mode "$MODE" is not supported."
        exit
    fi

    if ! [ -f soft/$CMSSW_VER/.installed ] ; then
        run_cmd mkdir -p soft
        run_cmd cd soft
        if [ -d $CMSSW_VER ] ; then
            echo "Removing incomplete $CMSSW_VER installation..."
            run_cmd rm -rf $CMSSW_VER
        fi
        echo "Creating new $CMSSW_VER area..."
        run_cmd scramv1 project CMSSW $CMSSW_VER
        run_cmd cd $CMSSW_VER/src
        run_cmd eval `scramv1 runtime -sh`

        run_cmd mkdir LLStaus_Run2
        run_cmd cd LLStaus_Run2
        run_cmd ln -s ../../../../Production Production
        run_cmd mkdir -p ../data 
        run_cmd cp -rf ../../../../Production/data/models/* ../data/
        run_cmd mkdir -p ../src/data
        run_cmd cp -rf ../../../../Production/data/models/* ../src/data
        run_cmd scram b -j8
        run_cmd touch ../../.installed
        run_cmd cd ../../../..
    else
        run_cmd cd soft/$CMSSW_VER/src
        run_cmd eval `scramv1 runtime -sh`
        run_cmd cd ../../..
    fi
elif [ $MODE = "conda" ] ; then
    CONDA=$(which conda 2>/dev/null)
    if [ x$CONDA = "x" -o x$CONDA = "x/usr/bin/conda" ] ; then
        PRIVATE_CONDA_INSTALL=$BASE_PATH/soft/conda
        if ! [ -f "$PRIVATE_CONDA_INSTALL/.installed" ] ; then
            if [ -d $PRIVATE_CONDA_INSTALL ] ; then
                echo "Removing incomplete private conda installation..."
                run_cmd rm -rf $PRIVATE_CONDA_INSTALL
            fi
            echo "Installing conda..."
            run_cmd mkdir -p soft
            run_cmd cd soft
            run_cmd curl https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh -o Miniconda2-latest-Linux-x86_64.sh
            run_cmd bash Miniconda2-latest-Linux-x86_64.sh -b -p $PRIVATE_CONDA_INSTALL
            run_cmd touch "$PRIVATE_CONDA_INSTALL/.installed"
        fi
        __conda_setup="$($PRIVATE_CONDA_INSTALL/bin/conda shell.${SHELL##*/} hook)"
        if [ $? -eq 0 ]; then
            eval "$__conda_setup"
        else
            if [ -f "$PRIVATE_CONDA_INSTALL/etc/profile.d/conda.sh" ]; then
                . "$PRIVATE_CONDA_INSTALL/etc/profile.d/conda.sh"
            else
                export PATH="$PRIVATE_CONDA_INSTALL/bin:$PATH"
            fi
        fi
        unset __conda_setup
    fi
    env_found=$(conda env list | grep -E '^llstau .*' | wc -l)
    if [ $env_found -ne 1 ]; then
        echo "Creating llstau environment..."
        run_cmd conda env create -f $BASE_PATH/conda-env.yaml
    fi

    run_cmd conda activate llstau

else
    echo echo "Mode "$MODE" is not supported."
    exit
fi
run_cmd cd $BASE_PATH
echo "$MODE environment is successfully loaded."
