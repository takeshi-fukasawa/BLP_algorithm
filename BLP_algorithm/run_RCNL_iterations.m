
%% Parameter settings

ns=I;
TOL=1e-12;
ITER_MAX=1000;
%ITER_MAX=100;
Newton_spec=0;

results_data=zeros(n_market,1);

results_V_contraction=zeros(n_market,3);
results_V_contraction_spectral=zeros(n_market,3);

results_V_new=zeros(n_market,3);
results_V_new_spectral=zeros(n_market,3);

results_BLP_contraction=zeros(n_market,3);
results_BLP_contraction_spectral=zeros(n_market,3);

results_BLP_new=zeros(n_market,3);
results_BLP_new_spectral=zeros(n_market,3);

results_s_mixed_spectral=zeros(n_market,3);
results_s_mixed=zeros(n_market,3);

results_s_conservative_spectral=zeros(n_market,3);
results_s_conservative=zeros(n_market,3);



for m=1:n_market

rng(m);
gpurng(m);

%run DGP.m
run DGP_ABLP.m


%% Compute product market share, given parameters
run solve_equil.m

if G==1
run model_conv_stat.m
end

DIST_MAT=zeros(ITER_MAX,3);

results_data(m,1)=round(median(S_0t_data(:)),3);


%%%%%%%%%%%%%%
%% BLP contraction mapping
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
if 1==0
%% joint s update (nest case applicable)
run run_joint_s_update_static.m
end

if G==1 %%%%%%%%%%%%%%%%%%

    %% V_update static
    tune_param=0;
    run run_V_update_static.m
    results_V_contraction=results_V;
    results_V_contraction_spectral=results_V_spectral;
    
    DISTMAT_V_0=DIST_MAT_V;
    ratio_delta_V_0=ratio_delta_V;

    tune_param=1;
    run run_V_update_static.m
    results_V_new=results_V;
    results_V_new_spectral=results_V_spectral;

    DISTMAT_V_1=DIST_MAT_V;
    ratio_delta_V_1=ratio_delta_V;


    if 1==0
    %% Kaloutpsidi method

    switch_r_spec=1;%%%%
    run run_Kalouptsidi_method_static.m
    results_s_mixed=results_s;
    %%results_s_mixed_spectral=results_s;

    switch_r_spec=0;%%%%
    run run_Kalouptsidi_method_static.m
    results_s_conservative=results_s;
    %%results_s_conservative_spectral=results_s;
    
    else
        results_s_mixed=[];
    end
    

end %%%%% G==1
end % n_market

if G==1
%results_temp=[mean(results_BLP);mean(results_BLP_spectral);...
%    mean(results_joint_s);mean(results_joint_s_spectral);...
%    mean(results_s);mean(results_s_spectral);...
%    ];

%results_temp=[mean(results_BLP);mean(results_BLP_spectral);...
%    mean(results_s);mean(results_s_spectral);...
%    mean(results_V);mean(results_V_spectral);...
%    ];

results_temp=[mean(results_BLP_contraction,1);mean(results_BLP_contraction_spectral,1);...
    mean(results_BLP_new,1);mean(results_BLP_new_spectral,1);...
    mean(results_V_contraction,1);mean(results_V_contraction_spectral,1);...
    mean(results_V_new,1);mean(results_V_new_spectral,1)];

if isempty(results_s_mixed)==0
results_temp=[results_temp;...
    mean(results_s_mixed,1);%mean(results_s_mixed_spectral,1);...
    mean(results_s_conservative,1);%mean(results_s_conservative_spectral,1);...
     ];

end

elseif G>=2
    results_temp=[mean(results_BLP_contraction,1);mean(results_BLP_contraction_spectral,1);...
    mean(results_BLP_new,1);mean(results_BLP_new_spectral,1)];
end

if ns>=10
    results_temp(:,2)=round(results_temp(:,2),3);
else
    results_temp(:,2)=round(results_temp(:,2),4);
end

results=[[beta_0 G J rho_est mean(results_data,1)].*ones(size(results_temp,1),1),...
    results_temp];

results=...
    [log10(max(s_i0t_ccp_true(:)))*ones(size(results,1),1),...
    log10(min(s_i0t_ccp_true(:)))*ones(size(results,1),1),...
    results];
