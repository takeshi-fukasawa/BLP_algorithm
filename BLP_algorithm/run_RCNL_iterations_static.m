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

summary_market(m,:)=[S_0t_data];

if rho_est==0
run model_conv_stat.m
end

results_data(m,1)=round(median(S_0t_data(:)),3);


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


    if 1==0
    %% Kalouptsidi method

    switch_r_spec=1;%%%%
    run run_Kalouptsidi_method_static.m
    results_r_mixed=results_s;
    %%results_r_mixed_spectral=results_s;

    switch_r_spec=0;%%%%
    run run_Kalouptsidi_method_static.m
    results_r_conservative=results_s;
    %%results_r_conservative_spectral=results_s;
    
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

    results_BLP_temp=[...
        results_table_func(results_BLP_0);...
        results_table_func(results_BLP_0_spectral);...
        results_table_func(results_BLP_0_SQUAREM);...
        results_table_func(results_BLP_1);...
        results_table_func(results_BLP_1_spectral);...
        results_table_func(results_BLP_1_SQUAREM)];

if rho_est==0

    results_V_temp=[...
        results_table_func(results_V_0);...
        results_table_func(results_V_0_spectral);...
        results_table_func(results_V_0_SQUAREM);...
        results_table_func(results_V_1);...
        results_table_func(results_V_1_spectral);...
        results_table_func(results_V_1_SQUAREM)];

    results_temp=[results_BLP_temp;results_V_temp];


if isempty(results_r_mixed)==0
    results_temp=[results_temp;...
    mean(results_r_mixed,1);%mean(results_r_mixed_spectral,1);...
    mean(results_r_conservative,1);%mean(results_r_conservative_spectral,1);...
     ];

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

if ns>=10
    results_temp(:,2)=round(results_temp(:,2),3);
else
    results_temp(:,2)=round(results_temp(:,2),4);
end

results_temp(:,4)=round(results_temp(:,4),2);%log(s_jt_predict)-log(S_jt_data)

results=[[beta_0 G J rho_est mean(results_data,1)].*ones(size(results_temp,1),1),...
    results_temp];

results=...
    [mean(summary_market,1).*ones(size(results,1),1),...
    results];
