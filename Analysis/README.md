# Analysis (STau 2018 UL setup) - Manual

## Main steps of fake rate estimation (MC-based):

1. For producing the histograms for WJets enriched region:
```sh
> python -m pepper.runproc stau_fake_rate.py ./configs/stau2018_fake_rate.json -o ./output_fakes/fake_rate_v1/ --statedata ./output_fakes/pepper_condor_fake_rate.coffea -i ./configs/setup_env.sh --condor 2000 --retries 10
```

2. To plot control plots
```sh
DIR_PATH=./output_fakes/fake_rate_v1/; python ./stau_plotter.py ./configs/stau2018_fake_rate_plotter.json ${DIR_PATH}/hists/hists.json --outdir ${DIR_PATH}/hists/plots_signal --cutflow ${DIR_PATH}/cutflows.json --mode 1D
```

3. For the ratio calculation (the binning is hardcoded inside the script)
```sh
python ./stau_rate_calculate.py ./configs/stau2018_fake_rate_plotter.json ./output_fakes/fake_rate_v1/hists/hists.json --outdir ./output_fakes/fake_rate_v1/fake_sf --cutflow ./output_fakes/fake_rate_v1/cutflows.json
```

### The list of related config files:
- stau2018_fake_rate_hist.json - histograms definition for fake rate estimate
- stau2018_fake_rate_plotter.json - config file for plotting script
- stau2018_fake_rate.json - config file for main stau_fake_rate processor


## CMDs

DIR_PATH=./output_fakes/fake_rate_v1/; python ./stau_plotter.py ./configs/stau2018_fake_rate_plotter.json ${DIR_PATH}/hists/hists.json --outdir ${DIR_PATH}/hists/plots_signal --cutflow ${DIR_PATH}/cutflows.json --mode 1D

DIR_PATH=./output_prompt/prompt_rate_v2/; python ./stau_plotter.py ./configs/stau2018_prompt_rate_plotter.json ${DIR_PATH}/hists/hists.json --outdir ${DIR_PATH}/hists/plots_signal --cutflow ${DIR_PATH}/cutflows.json --mode 1D

python ./stau_rate_calculate.py ./configs/stau2018_fake_rate_plotter.json ./output_fakes/fake_rate_v1/hists/hists.json --outdir ./output_fakes/fake_rate_v1/fake_sf --cutflow ./output_fakes/fake_rate_v1/cutflows.json     

python ./stau_rate_calculate.py ./configs/stau2018_prompt_rate_plotter.json ./output_prompt/prompt_rate_v2/hists/hists.json --outdir ./output_prompt/prompt_rate_v2/fake_sf --cutflow ./output_prompt/prompt_rate_v2/cutflows.json

## Data-driven backgraund prediction:

mumu region

```sh
python -m pepper.runproc stau_processor_ztomumu_FR.py ./configs/proc_2017/stau2017_ztomumu.json -o ./output_iteration_4/2017/output_zmumu/zmumu_v1/ --statedata ./output_iteration_4/2017/output_zmumu/zmumu_v1.coffea -i ./configs/setup_env_mamba.sh --metadata pepper_metadata_mamba_v2.pepper --condorlogdir pepper_logs_new --condor 400 --retries 20 -m 8 -R

DIR_MC=./output_iteration_4/2017/output_zmumu/zmumu_v1/; python ./stau_plotter.py ./configs/proc_2018/stau2018_ztomumu_plot_config.json ${DIR_MC}/hists/hists.json --outdir ${DIR_MC}/yield_prediction --data --cutflow ${DIR_MC}/cutflows.json -m prediction

```

## Macros:

Script to compare 1D fake rates in DY+jets and W+jets for 3 years:
```sh
python ./Analysis/macros/plot_fake_compare_dxy_dxy_run2.py
```

## Example of 2018UL yield prediction:

### Step 1: measure fake rates and calculate scalefactors

set up configuration in `stau2018_wjets.json`

```sh
python -m pepper.runproc stau_processor_wjets.py ./configs/stau2018_wjets.json -o ./output_iteration_2/output_wjet/wjet_v12 --statedata ./output_iteration_2/output_wjet/wjet_v12.coffea -i ./configs/setup_env_parsel.sh --condor 1000 --retries 10
```

setup binning and variables (`fake_rate`) and configuration in `stau2018_wjets_plotter.json`
```sh
DIR_MC=./output_iteration_4/2018/output_wjet/wjet_fake_v1/; python ./stau_rate_calculate.py ./configs/proc_2018/stau2018_wjets_plotter.json ${DIR_MC}/hists/hists.json --outdir ${DIR_MC}/fake_rate_ext --cutflow ${DIR_MC}/cutflows.json
```

### Step 2: run the processor with signal

setup the path to scale factors (`jet_fake_rate`), selections, data (`exp_datasets`) and signal path (`mc_datasets`) in `stau2018_signal.json`

```sh
python -m pepper.runproc stau_processor_signal.py ./configs/proc_2018_new/stau2018_signal.json -o ./output_iteration_4/2018/output_signal/signal_v1_signal_pass_v3/ --statedata ./output_iteration_4/2018/output_signal/signal_v1_signal_pass_v3.coffea -i ./configs/setup_env_mamba.sh --metadata pepper_metadata_mamba_v2.pepper --condorlogdir pepper_logs_new --condor 400 --retries 20 -m 8 -R

```

### Step 3: plot the prediction histograms

configuration of the prediction histograms can be done checked in the group `prediction_hist` in `stau2018_signal_plot_config.json` file
```sh
DIR_MC=./output_iteration_3/output_signal/signal_v23/; python ./stau_plotter.py ./configs/stau2018_signal_plot_config.json ${DIR_MC}/hists/hists.json --outdir ${DIR_MC}/plots_predict_unblind --cutflow ${DIR_MC}/cutflows.json -m prediction_sys
```