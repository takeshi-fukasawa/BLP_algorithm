%% Kalouptsidi method
r_initial0=log(S_0t_data.*weight);%1 by ns; Initial value

for method=1:3
    run spec_settings.m

    for kk=1:n_sim
        [output_r_spectral,other_output,iter_info]=...
    spectral_func(@r_update_RCL_func,spec,{r_initial0},...
    weight,mu_ijt_est,rho_est,...
    S_jt_data,S_0t_data,exp(mu_ijt_est),...
    switch_r_spec);

    r_sol=output_r_spectral{1};

    end%% kk

    r_sol=r_normalize_func(r_sol,S_0t_data);%%%%%%%%

    V_sol=log(weight)-r_sol;
    delta_sol=compute_delta_from_V_func(...
        mu_ijt_est,weight,S_jt_data,V_sol);

    [s_jt_predict,~]=...
      share_func(delta_sol+mu_ijt_est,zeros(1,ns,1,T),rho_est,weight);%J*1*G*T

    results_r_temp=results_output_func(iter_info,s_jt_predict,S_jt_data);
    
    if method==1
        results_r(m,:)=results_r_temp;
        iter_info_r=iter_info;
    elseif method==2
        results_r_spectral(m,:)=results_r_temp;
        iter_info_r_spectral=iter_info;
    elseif method==3
        results_r_SQUAREM(m,:)=results_r_temp;
        iter_info_r_SQUAREM=iter_info;
    end
end% method=1,2,3

