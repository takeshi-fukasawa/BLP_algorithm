global Pr0_spec
Pr0_spec=0; % not use Pr0 (only for dynamic BLP)

%% Parameter settings

ns=I;
Newton_spec=0;

t_dim_id=4;
TOL_DIST_s_jt=1e-12;

results_data=zeros(n_market,1);

for m=1:n_market
market_id=m

rng(m);
gpurng(m);

%run DGP.m
run DGP_ABLP.m


%% Compute product market share, given parameters
run solve_equil.m

summary_market(m,:)=round([S_0t_data],3);

if rho_est==0
run model_conv_stat.m
end

%%%%%%%%%%%%%%
%% BLP contraction mapping (G>=2 case allowed)
tune_param_BLP=0;
run run_BLP_contraction.m

tune_param_BLP=1;
run run_BLP_contraction.m


TOL=1e-12;


%%%%%%%%%%%%%
if rho_est==0 %%%%%%%%%%%%%%%%%%

    %% V_update static
    tune_param=0;
    run run_V_update.m

    tune_param=1;
    run run_V_update.m


    if I==2
    %% Kalouptsidi method

    switch_r_spec=1;%%%%
    %run run_Kalouptsidi_method_static.m
    %results_r_mixed=results_s;
    %%results_r_mixed_spectral=results_s;

    run run_r.m
    results_r_mixed=results_r;
    results_r_mixed_spectral=results_r_spectral;
    results_r_mixed_SQUAREM=results_r_spectral;

    switch_r_spec=0;%%%%
    run run_r.m
    results_r_conservative=results_r;
    results_r_conservative_spectral=results_r_spectral;
    results_r_conservative_SQUAREM=results_r_spectral;
    
    else
        results_r_mixed=[];
    end
    

else %%%%% G>=2

    %% IV_update static
    tune_param=0;
    run run_IV_update_static.m

    tune_param=1;
    run run_IV_update_static.m


end % rho_est==0 or others

end % n_market

if spec.compute_alpha_spec==3
    results_BLP_temp=[...
        results_table_func(results_BLP_0);...
        results_table_func(results_BLP_0_spectral);...
        results_table_func(results_BLP_0_SQUAREM);...
        results_table_func(results_BLP_1);...
        results_table_func(results_BLP_1_spectral);...
        results_table_func(results_BLP_1_SQUAREM)];

else
    results_BLP_temp=[...
        results_table_func(results_BLP_0_spectral);...
        results_table_func(results_BLP_0_SQUAREM);...
        results_table_func(results_BLP_1_spectral);...
        results_table_func(results_BLP_1_SQUAREM)];
end

if rho_est==0

    if spec.compute_alpha_spec==3
        results_V_temp=[...
            results_table_func(results_V_0);...
            results_table_func(results_V_0_spectral);...
            results_table_func(results_V_0_SQUAREM);...
            results_table_func(results_V_1);...
            results_table_func(results_V_1_spectral);...
            results_table_func(results_V_1_SQUAREM)];
    else
        results_V_temp=[...
            results_table_func(results_V_0_spectral);...
            results_table_func(results_V_0_SQUAREM);...
            results_table_func(results_V_1_spectral);...
            results_table_func(results_V_1_SQUAREM)];
    end
    results_temp=[results_BLP_temp;results_V_temp];


if isempty(results_r_mixed)==0
    results_temp=[results_temp;...
        results_table_func(results_r_mixed);...
        results_table_func(results_r_mixed_spectral);...
        results_table_func(results_r_mixed_SQUAREM);...
        results_table_func(results_r_conservative);...
        results_table_func(results_r_conservative_spectral);...
        results_table_func(results_r_mixed_SQUAREM)];
end

elseif rho_est>0
     results_IV_temp=[...
        results_table_func(results_IV_0);...
        results_table_func(results_IV_0_spectral);...
        results_table_func(results_IV_0_SQUAREM);...
        results_table_func(results_IV_1);...
        results_table_func(results_IV_1_spectral);...
        results_table_func(results_IV_1_SQUAREM)];
        results_temp=[results_BLP_temp;results_IV_temp];

end

results=[[rho_est beta_0 G J mean(summary_market,1)].*ones(size(results_temp,1),1),...
    results_temp];

results(:,end-1)=round(results(:,end-1),1);%%DIST


results(:,end)=results(:,end)*100;%%conv
results(:,end-2)=results(:,end-2)*100;%%conv


