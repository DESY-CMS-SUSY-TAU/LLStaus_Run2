import awkward as ak
import numpy as np
import numba as nb
import scipy
import glob, os
import plotext as plt
# import matplotlib.pyplot as plt


n_events = 1000000
scale_exp = 1.6
# scale_exp = 4.0

# simulate events with n-jets each having the score 0-1
def truncated_exp_OP(a,b,scale_exp,how_many):
    rands = np.random.exponential(scale=scale_exp, size=how_many)
    rands = rands[(rands>a) & (rands<b)]
    return rands


def binomial(p, n, k):
    # p -> scalar
    # n -> list of n
    # k -> dcalar (number of jets)
    n = n.astype(int)
    C = scipy.special.factorial(n) 
    C = C / scipy.special.factorial(k)
    C = C / scipy.special.factorial(n-k)
    prob = np.power(p, k) * np.power(1-p, n-k)
    return C * prob 
    

def gen_experiment():

    print("generating replica ...")


    n_jets = truncated_exp_OP(2,20, scale_exp, n_events).astype(int) # starts from 2 jets to 10
    # n_jets = np.full((n_events), 3)
    print("mean n-jets:", np.mean(n_jets))
    
    # generate score for each event
    score_p = []
    for jets in n_jets:
        # score_p.append( np.random.uniform(0.0,1.0,jets) )
        rands = np.random.triangular(0.0, min((20-jets)/40, 1.0), 1.0, jets)
        # rands = np.random.triangular(0.0, 16.8/40, 1.0, jets)
        score_p.append( rands )
    score_p = ak.Array(score_p)

    
    
    # jets = jets[jets.disTauTag_score1 > 0.17]
    # n_jets = ak.num(jets, axis=1)
    
    # jets = jets[n_jets == 2]
    # n_jets = ak.num(jets, axis=1)
    

    print("number of jets:", ak.sum(n_jets))
    print("mean n-jets:", np.mean(n_jets))
    
    # jets that pass the score
    # score_thrs = 0.01 # to be compatitive with our fake_rate
    # score_thrs = 1.0 - score_thrs
    score_thrs = 0.9972
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
    
    # rate_hist = FakeRate(jets, 0.9972)

    # count number of passes:
    events_bin0 = (n_jets_pass == 0)
    events_bin1 = (n_jets_pass == 1)
    events_bin2 = (n_jets_pass == 2)
    events_bin3 = (n_jets_pass == 3)
    events_overlow = (n_jets_pass > 3)

    jets_events_bin0 = score_p[ events_bin0 ]
    jets_events_bin1 = score_p[ events_bin1 ]
    jets_events_bin2 = score_p[ events_bin2 ]
    jets_events_bin3 = score_p[ events_bin3 ]
    jets_events_overlow = score_p[ events_overlow ]


    print("events_bin0:", ak.sum(events_bin0))
    print("events_bin1:", ak.sum(events_bin1))
    print("events_bin2:", ak.sum(events_bin2))
    print("events_bin3:", ak.sum(events_bin3))
    print("events_overlow:", ak.sum(events_overlow))
    

    # prediction from bin 0 to 1
    per_event_sfs_0to1 = ak.num(jets_events_bin0) * f
    predict_0to1 = ak.sum( per_event_sfs_0to1 )
    print("predict 0->1:", predict_0to1)

    # prediction from bin 0 to 2
    fake_rates_0 = ak.full_like(jets_events_bin0, f)
    combinations = ak.combinations(fake_rates_0, 2, axis=1) # to have all possible combinations of 2 jets
    combinations_unzipped = ak.unzip(combinations)
    products = combinations_unzipped[0] * combinations_unzipped[1]
    per_event_sfs_0to2 = ak.sum(products, axis=1)
    predict_0to2 = ak.sum(per_event_sfs_0to2, axis=0)
    print("predict 0->2:", predict_0to2)
    
    # prediction from bin 0 to 2 in a binomial way
    # fake_rates_0 = ak.full_like(jets_events_bin0, f)
    # fake_rates_0_njets = ak.num(fake_rates_0, axis = -1)
    # per_event_sfs_0to2_binom = binomial(f, ak.to_numpy(fake_rates_0_njets), 2)
    # predict_0to2_binom = ak.sum(per_event_sfs_0to2_binom, axis=0)
    # print("predict 0->2:", predict_0to2_binom, "(calculated in a binomial way)")
    
    # # prediction from bin 0 to 3
    # fake_rates_0 = ak.full_like(jets_events_bin0, f)
    # combinations = ak.combinations(fake_rates_0, 3, axis=1) # to have all possible combinations of 3 jets
    # combinations_unzipped = ak.unzip(combinations)
    # products = combinations_unzipped[0] * combinations_unzipped[1] * combinations_unzipped[2]
    # per_event_sfs_0to3 = ak.sum(products, axis=1)
    # predict_0to3 = ak.sum(per_event_sfs_0to3, axis=0)
    # print("predict 0->3:", predict_0to3)
    
    # prediction from bin 1 to 2
    jets_events_bin1_notag = jets_events_bin1[ (jets_events_bin1 < score_thrs) ] # events with 1 tag
    per_event_sfs_1to2 = ak.num(jets_events_bin1_notag, axis=-1) * f
    predict_1to2 = ak.sum( per_event_sfs_1to2 )
    print("predict 1->2:", predict_1to2, "corrected:", predict_1to2/2.0)
    
    '''
    # prediction from bin 1 to 2 calculated as a ratio of weights
    jets_score_bin1 = ak.full_like(jets_events_bin1, f)
    denominator = ak.sum(jets_score_bin1, axis=-1)
    combinations = ak.combinations(jets_score_bin1, 2, axis=1)
    combinations_unzipped = ak.unzip(combinations)
    products = combinations_unzipped[0] * combinations_unzipped[1]
    numerator = ak.sum(products, axis=1)
    corrected_predict_1to2 = numerator / denominator
    predict_corrected_predict_1to2 = ak.sum( corrected_predict_1to2 )
    print("corrected predict 1->2 ratio:", predict_corrected_predict_1to2)
    
    # prediction from bin 1 to 2 with extra factor two
    jets_events_bin1_notag = jets_events_bin1[ (jets_events_bin1 < (1.0 - score_thrs)) ] # events with 1 tag
    per_event_sfs_1to2 = ak.num(jets_events_bin1_notag, axis=-1) * f
    predict_corrected_predict_1to2_v2 = ak.sum( per_event_sfs_1to2 / 2.0)
    print("predict 1->2:", predict_corrected_predict_1to2_v2)
    
    per_event_sfs_1to2 = ak.num(jets_events_bin1, axis=-1) * f
    predict_corrected_predict_1to2_v3 = ak.sum( per_event_sfs_1to2 / 2.0)
    print("predict 1->2:", predict_corrected_predict_1to2_v3)
    '''
    
    # prediction from bin 1 to 3
    # jets_events_bin1_notag = jets_events_bin1[ (jets_events_bin1 < (1.0 - score_thrs)) ] # events with 1 tag
    # fake_rates_1 = ak.full_like(jets_events_bin1_notag, f)
    # combinations = ak.combinations(fake_rates_1, 2, axis=1) # to have all possible combinations of 3 jets
    # combinations_unzipped = ak.unzip(combinations)
    # products = combinations_unzipped[0] * combinations_unzipped[1]
    # per_event_sfs_1to3 = ak.sum(products, axis=1)
    # predict_1to3 = ak.sum( per_event_sfs_1to3 )
    # print("predict 1->3:", predict_1to3)
    
    # # prediction from bin 2 to 3
    # jets_events_bin2_notag = jets_events_bin2[ (jets_events_bin2 < (1.0 - score_thrs)) ] # events with 1 tag
    # per_event_sfs_2to3 = ak.num(jets_events_bin2_notag, axis=-1) * f
    # predict_2to3 = ak.sum( per_event_sfs_2to3 )
    # print("predict 2->3:", predict_2to3)
    
    print("Prediction for two jet events only:")
    alternative_from0to1 = ak.sum(events_bin0) * 2 * f / (1-f)
    alternative_from1to2 = ak.sum(events_bin1) * f / (2 * (1-f))
    alternative_from0to2 = ak.sum(events_bin0) * f**2 / ((1-f)**2)
    print("alternative 0->1:", alternative_from0to1)
    print("alternative 1->2:", alternative_from1to2)
    print("alternative 0->2:", alternative_from0to2)
    
    # print("Prediction for three jet events only:")
    # alternative_from0to1 = ak.sum(events_bin0) * 3 * f / (1-f)
    # alternative_from1to2 = ak.sum(events_bin1) * f / (1-f)
    # alternative_from0to2 = ak.sum(events_bin0) * 3 * f**2 / ((1-f)**2)
    # print("alternative 0->1:", alternative_from0to1)
    # print("alternative 1->2:", alternative_from1to2)
    # print("alternative 0->2:", alternative_from0to2)

    # return predict_0to2, predict_1to2, ak.sum(events_bin2), \
    #        predict_0to3, predict_2to3, ak.sum(events_bin3)

    return ak.sum(events_bin2), \
           predict_0to2, \
           predict_1to2/2.0

true, from0to2, from1to2 = [], [], []
for exp_replica in range(100):
    print("------------------->", exp_replica)
    pred = gen_experiment()
    true.append(pred[0])
    from0to2.append(pred[1])
    from1to2.append(pred[2])

m_true = np.mean(true)
m_from0to2 = np.mean(from0to2)
m_from1to2 = np.mean(from1to2)


plt.hist(true, bins=30, label = f"true bin 2 ({m_true})")
plt.hist(from0to2, bins=50, label = f"predict 0->2 ({m_from0to2})")
plt.hist(from1to2, bins=50, label = f"predict 1->2 ({m_from1to2})")

plt.legend(loc="upper right")
plt.show()
plt.savefig('predict_cor.png')