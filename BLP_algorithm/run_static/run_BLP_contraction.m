%% Compute delta, given nonlinear parameters (sigma)

if tune_param_BLP==0
    results_BLP=results_BLP_contraction;
    results_BLP_spectral=results_BLP_contraction_spectral;
else
    results_BLP=results_BLP_new;
    results_BLP_spectral=results_BLP_new_spectral;
end


    %% BLP spectral
    delta_initial0=log(S_jt_data)-log(S_0t_data)-rho_est.*log(S_jt_given_g_data);%J by 1; Initial value of delta

    for method=1:2
    spec=spec_default;

    if method==1
        spec.update_spec=0;% fixed point iteration
    end

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
    end
    end% method=1 or 2

    if tune_param_BLP==0
        results_BLP_contraction=results_BLP;
        results_BLP_contraction_spectral=results_BLP_spectral;
        iter_info_BLP_contraction=iter_info_BLP;
        iter_info_BLP_contraction_spectral=iter_info_BLP_spectral;
    else
        results_BLP_new=results_BLP;
        results_BLP_new_spectral=results_BLP_spectral;
        iter_info_BLP_new=iter_info_BLP;
        iter_info_BLP_new_spectral=iter_info_BLP_spectral;
    end

%%%%%%%%%%%%%%
