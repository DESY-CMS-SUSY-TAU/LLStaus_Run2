import awkward as ak
import numpy as np
import numba as nb
import scipy
import glob, os
# import plotext as plt
import matplotlib.pyplot as plt
from tqdm import tqdm

signal_factor = 365.7 * 59.74 / 2584850.0
loose_thr = 0.17
# score_thrs = 0.9972
score_thrs = 0.9900
dxy_min = 0.2

class FakeRate():

    def __init__(self, jets, thr):
        
        # self.pt_edges = np.array([30, 50,  70, 100, 200, 300, 10000])
        self.pt_edges = np.array([30, 40, 50, 60, 70, 100, 10000])
        # self.pt_edges = np.array([30, 10000])
        # self.eta_edges = np.array([-2.5, -1.5, -0.5, 0.5, 1.5, 2.5])
        self.eta_edges = np.array([-2.5, 2.5])
        self.threshold = thr
        
        pt = ak.flatten(jets.pt).to_numpy()
        eta = ak.flatten(jets.eta).to_numpy()
        
        self.hist_ll, _, _ = np.histogram2d(pt, eta, bins=(self.pt_edges, self.eta_edges))
        
        jets_upd = jets[jets.disTauTag_score1 >= thr]
        pt = ak.flatten(jets_upd.pt).to_numpy()
        eta = ak.flatten(jets_upd.eta).to_numpy()
        
        self.hist_tt, _, _ = np.histogram2d(pt, eta, (self.pt_edges, self.eta_edges))
        
        
        self.rate = self.hist_tt / self.hist_ll
        
    def print(self):
        
        print(self.rate)
        
    def get(self, pt, eta, sys="nom", method="list"):
        
        if method == "list":
            starts = ak.layout.Index64(np.array(pt.layout.starts))
            stops = ak.layout.Index64(np.array(pt.layout.stops))
            pt = pt.layout.content
            eta = eta.layout.content
        
        if sys=="up" or sys=="down":
            xidx = np.clip(np.digitize(pt, self.pt_edges), 0, self.rate.shape[0]-1) - 1
            yidx = np.clip(np.digitize(eta, self.eta_edges), 0, self.rate.shape[1]-1) - 1
            H = self.rate[xidx, yidx]
            if sys=="up":
                H += 0.3*H
            else:
                H -= 0.3*H
        else:
            xidx = np.clip(np.digitize(pt, self.pt_edges), 0, self.rate.shape[0]-1) - 1
            yidx = np.clip(np.digitize(eta, self.eta_edges), 0, self.rate.shape[1]-1) - 1
            H = self.rate[xidx, yidx]
        
        if method == "list":
            content = ak.layout.NumpyArray(np.array(H))
            return ak.Array(ak.layout.ListArray64(starts, stops, content))
        return H
    
    def add(self, jets):
        
        pt = ak.flatten(jets.pt).to_numpy()
        eta = ak.flatten(jets.eta).to_numpy()
        
        if len(jets.pt) == 0: return
        
        add_hist_ll, _, _ = np.histogram2d(pt, eta, bins=(self.pt_edges, self.eta_edges))
        jets_upd = jets[jets.disTauTag_score1 >= self.threshold]
        
        if len(jets_upd.pt) == 0: return
        
        pt = ak.flatten(jets_upd.pt).to_numpy()
        eta = ak.flatten(jets_upd.eta).to_numpy()
        add_hist_tt, _, _ = np.histogram2d(pt, eta, (self.pt_edges, self.eta_edges))
        
        self.hist_tt += add_hist_tt
        self.hist_ll += add_hist_ll
        
        self.rate = self.hist_tt / self.hist_ll

     
# Collect the data:
# arrays = []
# for file in glob.glob("jets_skims/DATA_MET/*.parquet"):
#     akw_arr = ak.from_parquet(file)
#     akw_save = ak.Array([{ "disTauTag_score1": akw_arr.disTauTag_score1, "pt": akw_arr.pt ,"eta": akw_arr.eta, "dxy": akw_arr.dxy}])
#     arrays.append(akw_save)
# jets_info = ak.concatenate(arrays, axis=0)
# score_p = jets_info.disTauTag_score1
# n_jets = ak.num(score_p, axis=-1)
# print("collected")
# ak.to_parquet(jets_info, "jets_skims/DATA_MET.parquet")
# exit()


score_plot = {"pass": [], "nopass": []}
# ------------------------------------ Zmumu ------------------------------------

jets = ak.from_parquet("./jets_skims/SingleMuon_data.parquet")

# plt.hist(ak.flatten(jets.disTauTag_score1), bins=50, label = f"Zmumu", range=[0,1], density=True, histtype="step")
jets = jets[abs(jets.dxy) > dxy_min]
jets = jets[(jets.disTauTag_score1 > loose_thr)]
n_jets = ak.num(jets.disTauTag_score1, axis=1)
jets = jets[n_jets == 2]
n_jets = ak.num(jets.disTauTag_score1, axis=1)
score_p = jets.disTauTag_score1

print("number of jets:", ak.sum(n_jets))
print("mean n-jets:", np.mean(n_jets))

# jets that pass the score
# score_thrs = 0.01 # to be compatitive with our fake_rate
# score_thrs = 1.0 - score_thrs
print("score threshold:", score_thrs)
pass_ = (score_p >= (score_thrs))
jets_pass = score_p[pass_]
jets_nopass = score_p[~pass_]

n_jets_pass = ak.num(jets_pass, axis=1)
n_pass = ak.sum(n_jets_pass)
n_jets_nopass = ak.num(jets_nopass, axis=1)
n_notpass = ak.sum(n_jets_nopass)
assert(n_pass != 0)

f = n_pass / ak.sum(n_jets)
print("Number of passing jets:", n_pass)
print("Number of not passing jets:", n_notpass)
print("Fake rate:", f)

jets_pass = jets[jets.disTauTag_score1 >= score_thrs]
n_jets_pass = ak.num(jets_pass.pt, axis=1)

jets_nopass = jets[jets.disTauTag_score1 < score_thrs]

score_plot["pass"].append(ak.flatten(jets_pass.pt))
score_plot["nopass"].append(ak.flatten(jets_nopass.pt))


# count number of passes:
events_bin0 = (n_jets_pass == 0)
events_bin1 = (n_jets_pass == 1)
events_bin2 = (n_jets_pass == 2)
events_overlow = (n_jets_pass > 2)


T1 = (jets.disTauTag_score1[:, 0] > score_thrs)
T2 = (jets.disTauTag_score1[:, 1] > score_thrs)

L1 = (jets.disTauTag_score1[:, 0] > loose_thr)
L2 = (jets.disTauTag_score1[:, 1] > loose_thr)


# print("f(LL->TL)", ak.sum(T1)*ak.sum(L2)/(ak.sum(L1)*ak.sum(L2)))
# print("f(TL->TT)", ak.sum(T1)*ak.sum(T2)/(ak.sum(T1)*ak.sum(L2)))

print("f(LL->TL)", ak.sum(T1&L2)/(ak.sum(L1&L2)))
print("f(TL->TT)", ak.sum(T1&T2)/(ak.sum(T1&L2)))

jets_events_bin0 = jets[ events_bin0 ]
jets_events_bin1 = jets[ events_bin1 ]
jets_events_bin2 = jets[ events_bin2 ]
jets_events_overlow = jets[ events_overlow ]


print("events_bin0:", ak.sum(events_bin0))
print("events_bin1:", ak.sum(events_bin1))
print("events_bin2:", ak.sum(events_bin2))
print("events_overlow:", ak.sum(events_overlow))

rate_zmumu = FakeRate(jets, score_thrs)
rate_zmumu.print()


# prediction from bin 0 to 1
# per_event_sfs_0to1 = ak.num(jets_events_bin0) * f
fake0 = rate_zmumu.get(jets_events_bin0.pt, jets_events_bin0.eta)
# print(jets_events_bin0.pt)
# print(jets_events_bin0.eta)
# print(fake0)
predict_0to1 = ak.sum( fake0, axis=1 )
print("predict 0->1:", ak.sum(predict_0to1))

# prediction from bin 0 to 2
# fake_rates_0 = ak.full_like(jets_events_bin0, f)
combinations = ak.combinations(fake0, 2, axis=1) # to have all possible combinations of 2 jets
combinations_unzipped = ak.unzip(combinations)
products = combinations_unzipped[0] * combinations_unzipped[1]
per_event_sfs_0to2 = ak.sum(products, axis=-1)
print("predict 0->2:", ak.sum(per_event_sfs_0to2, axis=0))

# prediction from bin 1 to 2
jets_events_bin1_notag = jets_events_bin1[ (jets_events_bin1.disTauTag_score1 < score_thrs) ] # events with 1 tag
fake1 = rate_zmumu.get(jets_events_bin1_notag.pt, jets_events_bin1_notag.eta)
per_event_sfs_1to2 = ak.sum(fake1, axis=-1)
predict_1to2 = ak.sum( per_event_sfs_1to2 )
print("predict 1->2:", predict_1to2, "corrected:", predict_1to2/2.0)

# prediction from bin 0 to 1
per_event_sfs_0to1 = ak.num(jets_events_bin0.pt) * f
predict_0to1 = ak.sum( per_event_sfs_0to1 )
print("f > predict 0->1:", predict_0to1)

# prediction from bin 0 to 2
fake_rates_0 = ak.full_like(jets_events_bin0.pt, f)
combinations = ak.combinations(fake_rates_0, 2, axis=1) # to have all possible combinations of 2 jets
combinations_unzipped = ak.unzip(combinations)
products = combinations_unzipped[0] * combinations_unzipped[1]
per_event_sfs_0to2 = ak.sum(products, axis=1)
predict_0to2 = ak.sum(per_event_sfs_0to2, axis=0)
print("f > predict 0->2:", predict_0to2)

# prediction from bin 1 to 2
jets_events_bin1_notag = jets_events_bin1[ (jets_events_bin1.disTauTag_score1 < score_thrs) ] # events with 1 tag
per_event_sfs_1to2 = ak.num(jets_events_bin1_notag.pt, axis=-1) * f
predict_1to2 = ak.sum( per_event_sfs_1to2 )
print("f > predict 1->2:", predict_1to2, "corrected:", predict_1to2/2.0)

# print("Prediction for two jet events only:")
# alternative_from0to1 = ak.sum(events_bin0) * 2 * f / (1-f)
# alternative_from1to2 = ak.sum(events_bin1) * f / (2 * (1-f))
# alternative_from0to2 = ak.sum(events_bin0) * f**2 / ((1-f)**2)
# print("alternative 0->1:", alternative_from0to1)
# print("alternative 1->2:", alternative_from1to2)
# print("alternative 0->2:", alternative_from0to2)


# from all
# pt_1, pt_2 = jets.pt[:,0], jets.pt[:,1]
# eta_1, eta_2 = jets.eta[:,0], jets.eta[:,1]
# f_1 = rate_zmumu.get(pt_1, eta_1)
# f_2 = rate_zmumu.get(pt_2, eta_2)

#from bin 0 to 1/2
pt_1, pt_2 = jets_events_bin0.pt[:,0], jets_events_bin0.pt[:,1]
eta_1, eta_2 = jets_events_bin0.eta[:,0], jets_events_bin0.eta[:,1]
f_1 = rate_zmumu.get(pt_1, eta_1, method="none")
f_2 = rate_zmumu.get(pt_2, eta_2, method="none")
from0to1 = ( f_1*(1-f_2) + f_2*(1-f_1) ) / ((1-f_2)*(1-f_1))
from0to2 = ( f_1*f_2 ) / ((1-f_2)*(1-f_1))

# from bin 1 to 2
pt_1, pt_2 = jets_events_bin1.pt[:,0], jets_events_bin1.pt[:,1]
eta_1, eta_2 = jets_events_bin1.eta[:,0], jets_events_bin1.eta[:,1]
f_1 = rate_zmumu.get(pt_1, eta_1, method="none")
f_2 = rate_zmumu.get(pt_2, eta_2, method="none")
from1to2 = ( f_1*f_2 ) / (f_1*(1-f_2) + f_2*(1-f_1))

# print("from all:", ak.sum(fromall))
print("fake 0->1:", ak.sum(from0to1))
print("fake 1->2:", ak.sum(from1to2))
print("fake 0->2:", ak.sum(from0to2))


plt.hist(ak.flatten(jets[pass_].pt), bins=50, label = f"pass", range=[0, 500], density=False, histtype="step")
plt.hist(ak.flatten(jets[~pass_].pt), bins=50, label = f"nopass", range=[0, 500], density=False, histtype="step")
plt.yscale('log')

plt.legend(loc="upper right")
plt.show()
plt.savefig('pt_pass_zmumu.png')

# exit()
# ------------------------------------ MET -------------------------------------

event_counter = {
    "bin0" : 0,
    "bin1" : 0,
    "bin2" : 0,
    "overlow" : 0
}
prediction_counter = {
    "from0to1" : [0, 0, 0],
    "from1to2" : [0, 0, 0],
    "from0to2" : [0, 0, 0]
}
data_MET = ak.from_parquet("./jets_skims/DATA_MET.parquet")
# print("Analyzing MET data")
rate_hist = None
# rate_hist = rate_zmumu

print("Analyzing MET data: counting rate")
for jets in tqdm(data_MET):
    jets = jets[abs(jets.dxy) > dxy_min]
    jets_selected = jets[jets.disTauTag_score1 > loose_thr]
    n_jets = ak.num(jets_selected.disTauTag_score1, axis=1)
    jets_selected = jets_selected[n_jets == 2]
    if rate_hist == None:
        rate_hist = FakeRate(jets_selected, score_thrs)
    else:
        rate_hist.add(jets_selected)
rate_hist.print()

# exit()

print("Analyzing MET data: counting events")
for jets in tqdm(data_MET):
    
    jets = jets[abs(jets.dxy) > dxy_min]
    if(len(jets.pt) == 0): continue
    
    jets_selected = jets[jets.disTauTag_score1 > loose_thr]
    n_jets = ak.num(jets_selected.disTauTag_score1, axis=1)
    if(len(jets_selected.pt) == 0): continue
    
    jets_selected = jets_selected[n_jets == 2]
    n_jets = ak.num(jets_selected.disTauTag_score1, axis=1)
    if(len(jets_selected.pt) == 0): continue

    jets_pass = jets_selected[jets_selected.disTauTag_score1 >= score_thrs]
    n_jets_pass = ak.num(jets_pass.pt, axis=1)
    
    jets_nopass = jets_selected[jets_selected.disTauTag_score1 < score_thrs]

    score_plot["pass"].append(ak.flatten(jets_pass.pt))
    score_plot["nopass"].append(ak.flatten(jets_nopass.pt))
    
    # count number of passes:
    events_bin0 = (n_jets_pass == 0)
    events_bin1 = (n_jets_pass == 1)
    events_bin2 = (n_jets_pass == 2)
    events_overlow = (n_jets_pass > 2)

    jets_events_bin0 = jets_selected[ events_bin0 ]
    jets_events_bin1 = jets_selected[ events_bin1 ]
    jets_events_bin2 = jets_selected[ events_bin2 ]
    jets_events_overlow = jets_selected[ events_overlow ]

    event_counter["bin0"] += ak.sum(events_bin0)
    event_counter["bin1"] += ak.sum(events_bin1)
    event_counter["bin2"] += ak.sum(events_bin2)
    event_counter["overlow"] += ak.sum(events_overlow)
    
    # from bin 0 to 1/2
    from0to1 = [[0], [0], [0]]
    from0to2 = [[0], [0], [0]]
    if len(jets_events_bin0.pt) != 0:
        for i, estim in enumerate(["up", "nom","down"]):
            pt_1, pt_2 = jets_events_bin0.pt[:,0], jets_events_bin0.pt[:,1]
            eta_1, eta_2 = jets_events_bin0.eta[:,0], jets_events_bin0.eta[:,1]
            f_1 = rate_hist.get(pt_1, eta_1, sys=estim, method="none")
            f_2 = rate_hist.get(pt_2, eta_2, sys=estim, method="none")
            from0to1[i] = ( f_1*(1-f_2) + f_2*(1-f_1) ) / ((1-f_2)*(1-f_1))
            from0to2[i] = ( f_1*f_2 ) / ((1-f_2)*(1-f_1))
            
            # fake0 = rate_zmumu.get(jets_events_bin0.pt, jets_events_bin0.eta,  sys=estim, method="list")
            # predict_0to1 = ak.sum( fake0, axis=1 )
            # # print("predict 0->1:", ak.sum(predict_0to1))
            # from0to1[i] =  predict_0to1

            # # prediction from bin 0 to 2
            # combinations = ak.combinations(fake0, 2, axis=1) # to have all possible combinations of 2 jets
            # combinations_unzipped = ak.unzip(combinations)
            # products = combinations_unzipped[0] * combinations_unzipped[1]
            # per_event_sfs_0to2 = ak.sum(products, axis=-1)
            # # print("predict 0->2:", ak.sum(per_event_sfs_0to2, axis=0))
            # from0to2[i] = per_event_sfs_0to2


    # from bin 1 to 2
    from1to2 = [[0], [0], [0]]
    if len(jets_events_bin1.pt) != 0:
        for i, estim in enumerate(["up", "nom","down"]):
            # jets_events_bin1 = jets_events_bin1[jets_events_bin1.disTauTag_score1 < score_thrs]
            
            pt_1, pt_2 = jets_events_bin1.pt[:,0], jets_events_bin1.pt[:,1]
            eta_1, eta_2 = jets_events_bin1.eta[:,0], jets_events_bin1.eta[:,1]
            f_1 = rate_hist.get(pt_1, eta_1, sys=estim, method="none")
            f_2 = rate_hist.get(pt_2, eta_2, sys=estim, method="none")
            from1to2[i] = ( f_1*f_2 ) / (f_1*(1-f_2) + f_2*(1-f_1))
            
            # # prediction from bin 1 to 2
            # jets_events_bin1_notag = jets_events_bin1[ (jets_events_bin1.disTauTag_score1 < score_thrs) ] # events with 1 tag
            # fake1 = rate_zmumu.get(jets_events_bin1_notag.pt, jets_events_bin1_notag.eta,  sys=estim, method="list")
            # per_event_sfs_1to2 = ak.sum(fake1, axis=-1)
            # # predict_1to2 = ak.sum( per_event_sfs_1to2 )
            # # print("predict 1->2:", predict_1to2, "corrected:", predict_1to2/2.0)
            # from1to2[i] = per_event_sfs_1to2/2.0


    for i, estim in enumerate(["up", "nom","down"]): 
        prediction_counter["from0to1"][i] += ak.sum(from0to1[i])
        prediction_counter["from1to2"][i] += ak.sum(from1to2[i])
        prediction_counter["from0to2"][i] += ak.sum(from0to2[i])

print("MET stats:")
print(event_counter)

print("MET prediction:")
print(prediction_counter)

score_plot["pass"] = ak.concatenate(score_plot["pass"])
plt.hist(score_plot["pass"], bins=50, label = f"pass", range=[0, 500], density=True, histtype="step")
score_plot["nopass"] = ak.concatenate(score_plot["nopass"])
plt.hist(score_plot["nopass"], bins=50, label = f"nopass", range=[0, 500], density=True, histtype="step")
plt.yscale('log')

plt.legend(loc="upper right")
plt.show()
plt.savefig('pt_pass2.png')

# exit(0)
# ------------------------------------ Signal -------------------------------------

print("Analyzing signal data:")
event_counter = {
    "bin0" : 0,
    "bin1" : 0,
    "bin2" : 0,
    "overlow" : 0
}
data_signal = ak.from_parquet("./jets_skims/signal.parquet")
for jets in tqdm(data_signal):

    jets = jets[abs(jets.dxy) > dxy_min]

    if(len(jets.pt) == 0): continue

    jets = jets[jets.disTauTag_score1 > loose_thr]
    n_jets = ak.num(jets.disTauTag_score1, axis=1)
    jets = jets[n_jets == 2]
    
    if(len(jets.pt) == 0): continue

    n_jets = ak.num(jets.disTauTag_score1, axis=1)
    score_p = jets.disTauTag_score1

    if(len(jets.pt) == 0): continue

    pass_ = (score_p >= (score_thrs))
    jets_pass = score_p[pass_]
    jets_nopass = score_p[~pass_]

    n_jets_pass = ak.num(jets_pass, axis=1)
    n_pass = ak.sum(n_jets_pass)
    n_jets_nopass = ak.num(jets_nopass, axis=1)
    n_notpass = ak.sum(n_jets_nopass)

    # count number of passes:
    events_bin0 = (n_jets_pass == 0)
    events_bin1 = (n_jets_pass == 1)
    events_bin2 = (n_jets_pass == 2)
    events_overlow = (n_jets_pass > 2)

    # jets_events_bin0 = score_p[ events_bin0 ]
    # jets_events_bin1 = score_p[ events_bin1 ]
    # jets_events_bin2 = score_p[ events_bin2 ]
    # jets_events_overlow = score_p[ events_overlow ]

    event_counter["bin0"] += ak.sum(events_bin0) * signal_factor
    event_counter["bin1"] += ak.sum(events_bin1) * signal_factor
    event_counter["bin2"] += ak.sum(events_bin2) * signal_factor
    event_counter["overlow"] += ak.sum(events_overlow) * signal_factor

print("Signal stats:")
print(event_counter)
