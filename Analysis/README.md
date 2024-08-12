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