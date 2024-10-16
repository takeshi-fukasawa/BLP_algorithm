
n_col_results=5;

TOL_DIST_s_jt=1e-12;
t_dim_id=4;

n_dim_V=T+n_grid_IV;

for m=1:n_market
%for m=n_market:n_market%%%%%
%for m=13:13%%%%%
        
market_id=m

rng(m);
gpurng(m);

%% Generate data and solve equilibrium
run DGP.m

run solve_equil.m

%V_initial0=-log(S_0t_data.*weight);
delta_initial0=log(S_jt_data)-log(S_0t_data)-rho_est.*log(S_jt_given_g_data);% Initial value of delta

%% Update V (rho==0 & G==1 case only)
if G==1 & rho_est==0
    %%% tune_param==0
    tune_param=0;Newton_spec=0;
    run run_V_update.m
    
    %%% tune_param==1
    tune_param=1;Newton_spec=0;
    
    run run_V_update.m
    
    
    if T==1
        %%% tune_param==1/(1-beta_C)
        tune_param=1/(1-beta_C);Newton_spec=0;
        run run_V_update.m
    end

    
    %%% Newton iteration
    if 1==0
        tune_param=1;Newton_spec=1;
        run run_V_update.m
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

results_table_V=[...
    results_table_func(results_V_0);...
    results_table_func(results_V_0_Anderson);...
    results_table_func(results_V_0_spectral);...
    results_table_func(results_V_0_SQUAREM);...
    results_table_func(results_V_1);...
    results_table_func(results_V_1_Anderson);...
    results_table_func(results_V_1_spectral);...
    results_table_func(results_V_1_SQUAREM)...
    ];
%if T==1
%    results_table_V=[results_table_V;...
%        results_V_2;results_V_2_spectral;results_V_2_SQUAREM];
%end

results_table_V_BLP=[...
    results_table_func(results_V_BLP_0);...
    results_table_func(results_V_BLP_0_Anderson);...    
    results_table_func(results_V_BLP_0_spectral);...
    results_table_func(results_V_BLP_0_SQUAREM);...
    results_table_func(results_V_BLP_1);...
    results_table_func(results_V_BLP_1_Anderson);...
    results_table_func(results_V_BLP_1_spectral);...
    results_table_func(results_V_BLP_1_SQUAREM)];

results_table=[results_table_V;results_table_V_BLP];

results_table(:,end)=results_table(:,end)*100;%%conv
results_table(:,end-1)=round(results_table(:,end-1),1);%%DIST

results_table(:,end-2)=results_table(:,end-2)*100;%%conv

%%%%%% Test traditional nested loop %%%%%
if 1==0
    run run_BLP_middle_update_dynamic.m
end

%%%%% Test time-dependent or not (perfect foresight case)
if IVS_spec==0 & m==n_market & 1==1
    results_V_delta_t_dep=[results_V_0_spectral(end,:)];

    t_dependent_alpha_spec=0;
    [output_spectral,other_vars,iter_info]=...
        spectral_func(@V_update_func,spec,...
        {V_initial0},...
        weight,mu_ijt_est,...
        S_jt_data,S_0t_data,...
        weight_V,x_V,beta_C,tune_param,Newton_spec);

    V_sol=output_spectral{1};
    delta_sol=other_vars.delta_jt;

    IV=other_vars.IV;
    v_i0t_tilde=other_vars.v_i0t_tilde;
    v_ijt_tilde=other_vars.v_ijt_tilde;
    
    if T==size(V_sol,4)
        s_jt_predict=...
            share_func(v_ijt_tilde,...
            v_i0t_tilde,rho_est,weight);
    else
        s_jt_predict=...
            share_func(v_ijt_tilde(:,:,:,1:T),...
            v_i0t_tilde(:,:,:,1:T),rho_est,weight);
    end
    results_V_delta_t_indep=results_output_func(iter_info,s_jt_predict,S_jt_data);
 
    results_comparison_t_dep_alpha=[results_V_delta_t_dep;results_V_delta_t_indep]

    filename=append(output_path,"dynamic_BLP_IVS_results_beta_",...
        string(beta_C),"_",string(mistake_spec),"_V_mapping_comparison_t_dep_alpha.csv");
    writematrix(results_comparison_t_dep_alpha,filename)

end

