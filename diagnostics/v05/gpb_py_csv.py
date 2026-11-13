import sys, os, numpy as np, gpboost as gpb
f, lik, nrounds, slope, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4] == "1", sys.argv[5]
L = lambda n: np.loadtxt(os.path.join(f, n), delimiter=",", ndmin=1)
Xtr, Xte = np.loadtxt(os.path.join(f, "X_tr.csv"), delimiter=",", ndmin=2), np.loadtxt(os.path.join(f, "X_te.csv"), delimiter=",", ndmin=2)
kw = dict(group_data=L("id_tr.csv"), likelihood=lik)
if slope: kw.update(group_rand_coef_data=L("t_tr.csv"), ind_effect_group_rand_coef=[1])
gp = gpb.GPModel(**kw)
bst = gpb.train(params=dict(learning_rate=0.05, max_depth=2, num_leaves=4, min_data_in_leaf=2,
                            feature_fraction=1.0, num_threads=1, verbose=-1),
                train_set=gpb.Dataset(Xtr, L("y_tr.csv")), gp_model=gp, num_boost_round=nrounds)
pk = dict(data=Xte, group_data_pred=L("id_te.csv"), pred_latent=False)
if slope: pk.update(group_rand_coef_data_pred=L("t_te.csv"))
np.savetxt(out, np.asarray(bst.predict(**pk)["response_mean"]), delimiter=",")
