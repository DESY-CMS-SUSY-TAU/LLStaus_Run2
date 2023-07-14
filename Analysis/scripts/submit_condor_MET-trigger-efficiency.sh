python -u -m pepper.runproc \
evaluate_MET-trigger-efficiency.py ./configs/config_trigger_eff_MET.json \
-o ./output/pepper_condor/MET-trigger-efficiency/ \
--statedata ./output/pepper_condor/MET-trigger-efficiency/pepper_condor_fake_rate.coffea \
-i /nfs/dust/cms/user/sobhatta/work/LongLivedStaus/LLStaus_Run2/Analysis/configs/setup_env_sobhatta.sh \
--condor 2000 \
--retries 10