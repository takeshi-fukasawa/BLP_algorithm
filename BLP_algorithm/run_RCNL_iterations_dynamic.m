global Pr0_spec
Pr0_spec=0;

n_col_results=5;

results_V_0=zeros(n_market,n_col_results);
results_V_0_spectral=zeros(n_market,n_col_results);
results_V_1=zeros(n_market,n_col_results);
results_V_1_spectral=zeros(n_market,n_col_results);
results_V_2=zeros(n_market,n_col_results);
results_V_2_spectral=zeros(n_market,n_col_results);

results_V_BLP_0=zeros(n_market,n_col_results);
results_V_BLP_0_spectral=zeros(n_market,n_col_results);
results_V_BLP_1=zeros(n_market,n_col_results);
results_V_BLP_1_spectral=zeros(n_market,n_col_results);

TOL_DIST_s_jt=1e-12;
t_dim_id=4;

for m=1:n_market

rng(m);
gpurng(m);

%% Generate data and solve equilibrium
%run DGP.m
run DGP_ABLP.m
mu_ij_est=mu_ijt_true*1;

run solve_equil.m

%V_initial0=-log(S_0t_data.*weight);
delta_initial0=log(S_jt_data)-log(S_0t_data)-rho_est.*log(S_jt_given_g_data);% Initial value of delta

%% Update V (rho==0 & G==1 case only)
if G==1 & rho_est==0
    %%% tune_param==0
    tune_param=0;Newton_spec=0;
    run run_V_update_dynamic.m
    
    
    %%% tune_param==1
    tune_param=1;Newton_spec=0;
    
    run run_V_update_dynamic.m
    
    
    %%% tune_param==1/(1-beta_C)
    tune_param=1/(1-beta_C);Newton_spec=0;
    run run_V_update_dynamic.m
    
    
    %%% Newton iteration
    if 1==0
        tune_param=1;Newton_spec=1;
        run run_V_update_dynamic.m
     end

end


%% Jointly Update V and IV
if G>=2
run run_V_IV_update_dynamic.m
end

%% Joint update of delta and V
%%% tune_param_BLP==0
tune_param_BLP=0;
run run_delta_V_joint_update_dynamic.m

%%% tune_param_BLP==1
tune_param_BLP=1;
run run_delta_V_joint_update_dynamic.m

end

if mistake_spec==1
    mu_ijt=mu_ijt_est;
    run run_Bellman.m %%%%%
end


results_table=[...
    results_V_0;results_V_0_spectral;...
    results_V_1;results_V_1_spectral;...
    results_V_2;results_V_2_spectral;...
    results_V_BLP_0;results_V_BLP_0_spectral;...
    results_V_BLP_1;results_V_BLP_1_spectral];


if 1==0
    filename=append(save_path,"no_nest_results_beta_",...
        string(beta_C),"_",string(mistake_spec),".csv");

writematrix(round(results_table,3),filename)
end

