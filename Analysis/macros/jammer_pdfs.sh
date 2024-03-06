
DIR_MC=./output_iteration_2/output_wjet/wjet_v12_predict_fromWjet/; python ./stau_plotter.py ./configs/stau2018_wjets_plotter.json ${DIR_MC}/hists/hists.json --outdir ${DIR_MC}/plots_predict --cutflow ${DIR_MC}/cutflows.json -m prediction;
DIR_MC=./output_iteration_2/output_wjet/wjet_v12_predict_1d/; python ./stau_plotter.py ./configs/stau2018_wjets_plotter.json ${DIR_MC}/hists/hists.json --outdir ${DIR_MC}/plots_predict --cutflow ${DIR_MC}/cutflows.json -m prediction;
DIR_MC=./output_iteration_2/output_zmumu/zmumu_v4_predict_2d; python ./stau_plotter.py ./configs/stau2018_ztomumu_plot_config.json ${DIR_MC}/hists/hists.json --outdir ${DIR_MC}/plots_predict --cutflow ${DIR_MC}/cutflows.json -m prediction

DIR_MC=./output_iteration_2/output_signal/signal_v11; python ./stau_plotter.py ./configs/stau2018_signal_plot_config.json ${DIR_MC}/hists/hists.json --outdir ${DIR_MC}/plots_predict_unblind --data --cutflow ${DIR_MC}/cutflows.json -m prediction

PRED="bin0to1";
pdfjam ~/sshfs/output_iteration_2/empty_canvas.pdf ./mt_sum_${PRED}.pdf ./jet1_pt_${PRED}.pdf ./jet1_dxy_${PRED}.pdf ./METpt_${PRED}.pdf ./mt2_j1_j2_MET_${PRED}.pdf ./jet2_pt_${PRED}.pdf ./jet2_dxy_${PRED}.pdf --nup 4x2 --papersize '{20cm,9cm}' --outfile merge_${PRED}.pdf;
pdftoppm -r 300 -png -cropbox merge_bin0to1.pdf merge_bin0tp1;

PRED="bin1to2";
pdfjam ~/sshfs/output_iteration_2/empty_canvas.pdf ./mt_sum_${PRED}.pdf ./jet1_pt_${PRED}.pdf ./jet1_dxy_${PRED}.pdf ./METpt_${PRED}.pdf ./mt2_j1_j2_MET_${PRED}.pdf ./jet2_pt_${PRED}.pdf ./jet2_dxy_${PRED}.pdf --nup 4x2 --papersize '{60cm,28cm}' --outfile merge_${PRED}.pdf;
PRED="bin0to2";
pdfjam ~/sshfs/output_iteration_2/empty_canvas.pdf ./mt_sum_${PRED}.pdf ./jet1_pt_${PRED}.pdf ./jet1_dxy_${PRED}.pdf ./METpt_${PRED}.pdf ./mt2_j1_j2_MET_${PRED}.pdf ./jet2_pt_${PRED}.pdf ./jet2_dxy_${PRED}.pdf --nup 4x2 --papersize '{60cm,28cm}' --outfile merge_${PRED}.pdf;

pdftoppm -r 300 -png -cropbox merge_bin1to2.pdf merge_bin1to2;
pdftoppm -r 300 -png -cropbox merge_bin0to2.pdf merge_bin0to2;



