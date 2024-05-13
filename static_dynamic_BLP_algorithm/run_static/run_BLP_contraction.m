%% Compute delta, given nonlinear parameters (sigma)

for method=1:3
    run spec_settings.m

    for kk=1:n_sim
        [output_BLP_spectral,other_output_k,iter_info]=...
    spectral_func(@BLP_update_func,spec,{delta_initial0},...
    weight,mu_ijt_est,rho_est,S_jt_data,tune_param_BLP);

    delta_sol=output_BLP_spectral{1};

    end%% kk

    [s_jt_predict,~]=...
      share_func(delta_sol+mu_ijt_est,zeros(1,ns,1,T),rho_est,weight);%J*1*G*T

    results_BLP_temp=results_output_func(iter_info,s_jt_predict,S_jt_data);
    
    if method==1
        results_BLP=results_BLP_temp;
        iter_info_BLP=iter_info;
    elseif method==2
        results_BLP_spectral=results_BLP_temp;
        iter_info_BLP_spectral=iter_info;
    elseif method==3
        results_BLP_SQUAREM=results_BLP_temp;
        iter_info_BLP_SQUAREM=iter_info;
    end
end% method=1,2,3

if tune_param_BLP==0
    results_BLP_0(m,:)=results_BLP;
    iter_info_BLP_0=iter_info_BLP;
    results_BLP_0_spectral(m,:)=results_BLP_spectral;
    iter_info_BLP_0_spectral=iter_info_BLP_spectral;
    results_BLP_0_SQUAREM(m,:)=results_BLP_SQUAREM;
    iter_info_BLP_0_SQUAREM=iter_info_BLP_SQUAREM;
elseif tune_param_BLP==1
    results_BLP_1(m,:)=results_BLP;
    iter_info_BLP_1=iter_info_BLP;
    results_BLP_1_spectral(m,:)=results_BLP_spectral;
    iter_info_BLP_1_spectral=iter_info_BLP_spectral;
    results_BLP_1_SQUAREM(m,:)=results_BLP_SQUAREM;
    iter_info_BLP_1_SQUAREM=iter_info_BLP_SQUAREM;
end

