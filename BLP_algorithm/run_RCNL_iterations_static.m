global Pr0_spec
Pr0_spec=0; % not use Pr0 (only for dynamic BLP)

%% Parameter settings

ns=I;
TOL=1e-12;
ITER_MAX=1000;
%ITER_MAX=100;
Newton_spec=0;

t_dim_id=4;
TOL_DIST_s_jt=1e-12;

results_data=zeros(n_market,1);

if rho_est==0
    results_V_contraction=zeros(n_market,3);
    results_V_contraction_spectral=zeros(n_market,3);
    results_V_new=zeros(n_market,3);
    results_V_new_spectral=zeros(n_market,3);
else
    results_IV_contraction=zeros(n_market,3);
    results_IV_contraction_spectral=zeros(n_market,3);
    results_IV_new=zeros(n_market,3);
    results_IV_new_spectral=zeros(n_market,3);
end

results_BLP_contraction=zeros(n_market,3);
results_BLP_contraction_spectral=zeros(n_market,3);
results_BLP_new=zeros(n_market,3);
results_BLP_new_spectral=zeros(n_market,3);


results_r_mixed_spectral=zeros(n_market,3);
results_r_mixed=zeros(n_market,3);
results_r_conservative_spectral=zeros(n_market,3);
results_r_conservative=zeros(n_market,3);



for m=1:n_market

rng(m);
gpurng(m);

%run DGP.m
run DGP_ABLP.m


%% Compute product market share, given parameters
run solve_equil.m

if rho_est==0
run model_conv_stat.m
end

DIST_MAT=zeros(ITER_MAX,3);

results_data(m,1)=round(median(S_0t_data(:)),3);


%%%%%%%%%%%%%%
%% BLP contraction mapping (G>=2 case allowed)
tune_param_BLP=0;
run run_BLP_contraction.m

DISTMAT_BLP_0=DIST_MAT_BLP;
ratio_delta_BLP_0=ratio_delta_BLP;

tune_param_BLP=1;
run run_BLP_contraction.m
DISTMAT_BLP_1=DIST_MAT_BLP;
ratio_delta_BLP_1=ratio_delta_BLP;


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
    results_IV_contraction=results_IV;
    results_IV_contraction_spectral=results_IV_spectral;
    
    DISTMAT_IV_0=DIST_MAT_IV;
    ratio_delta_IV_0=ratio_delta_IV;

    tune_param=1;
    run run_IV_update_static.m
    results_IV_new=results_IV;
    results_IV_new_spectral=results_IV_spectral;

    DISTMAT_IV_1=DIST_MAT_IV;
    ratio_delta_IV_1=ratio_delta_IV;


end % rho_est==0 or others

end % n_market

if rho_est==0

results_temp=[mean(results_BLP_contraction,1);mean(results_BLP_contraction_spectral,1);...
    mean(results_BLP_new,1);mean(results_BLP_new_spectral,1);...
    mean(results_V_0,1);mean(results_V_0_spectral,1);...
    mean(results_V_1,1);mean(results_V_1_spectral,1)];

if isempty(results_r_mixed)==0
    results_temp=[results_temp;...
    mean(results_r_mixed,1);%mean(results_r_mixed_spectral,1);...
    mean(results_r_conservative,1);%mean(results_r_conservative_spectral,1);...
     ];

end

elseif rho_est>0
    results_temp=[mean(results_BLP_contraction,1);mean(results_BLP_contraction_spectral,1);...
    mean(results_BLP_new,1);mean(results_BLP_new_spectral,1);...
    mean(results_IV_contraction,1);mean(results_IV_contraction_spectral,1);...
    mean(results_IV_new,1);mean(results_IV_new_spectral,1)];
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
    [S_0t_data*ones(size(results,1),1),...
    log10(max(s_i0t_ccp_true(:)))*ones(size(results,1),1),...
    log10(min(s_i0t_ccp_true(:)))*ones(size(results,1),1),...
    results];
