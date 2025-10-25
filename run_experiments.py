
exec(open("test_nevo.py").read())
exec(open("test_BLP.py").read())

n_trial=3
for trial_id in range(n_trial+1):
    #exec(open("test_nevo.py").read())
    exec(open("run_Nevo_est.py").read())
    exec(open("run_BLP_est.py").read())


