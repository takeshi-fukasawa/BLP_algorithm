V_EV_initial0=cat(4,V_initial0,zeros(1,ns,1,T));

%% V_EV_update_func 
for method=1:3
    run spec_settings.m

[output_spectral,other_vars,iter_info]=...
        spectral_func(@V_EV_update_func,spec,...
        {V_EV_initial0},...
        weight,mu_ijt_est,...
        S_jt_data,S_0t_data,...
        weight_V,x_V,beta_C,tune_param,Newton_spec);

    V_EV_sol=output_spectral{1};
    delta_sol=other_vars.delta_jt;

    IV=other_vars.IV;
    v_i0t_tilde=other_vars.v_i0t_tilde;
    v_ijt_tilde=other_vars.v_ijt_tilde;
  
     s_jt_predict=...
            share_func(v_ijt_tilde(:,:,:,1:T),...
            v_i0t_tilde(:,:,:,1:T),rho_est,weight);
      results_V_EV=results_output_func(iter_info,s_jt_predict,S_jt_data);
    ratio_delta=delta_sol./delta_jt_true;


if method==1

    if tune_param==0
        results_V_EV_0(m,:)=results_V_EV;
        ratio_delta_V_EV_0=ratio_delta;
        iter_info_V_EV_0=iter_info;
    elseif tune_param==1
        results_V_EV_1(m,:)=results_V_EV;
        ratio_delta_V_EV_1=ratio_delta;
        iter_info_V_EV_1=iter_info;
    elseif tune_param>1
        results_V_EV_2(m,:)=results_V_EV;
        ratio_delta_V_EV_2=ratio_delta;
        iter_info_V_EV_2=iter_info;
    end

elseif method==2
    if tune_param==0
        results_V_EV_0_spectral(m,:)=results_V_EV;
        ratio_delta_V_EV_0_spectral=ratio_delta;
        iter_info_V_EV_0_spectral=iter_info;
    elseif tune_param==1
        results_V_EV_1_spectral(m,:)=results_V_EV;
        ratio_delta_V_EV_1_spectral=ratio_delta;
        iter_info_V_EV_1_spectral=iter_info;
    elseif tune_param>1
        results_V_EV_2_spectral(m,:)=results_V_EV;
        ratio_delta_V_EV_2_spectral=ratio_delta;
        iter_info_V_EV_2_spectral=iter_info;
    end% tune_param==0 or 1 or others
elseif method==3
    if tune_param==0
        results_V_EV_0_SQUAREM(m,:)=results_V_EV;
        ratio_delta_V_EV_0_SQUAREM=ratio_delta;
        iter_info_V_EV_0_SQUAREM=iter_info;
    elseif tune_param==1
        results_V_EV_1_SQUAREM(m,:)=results_V_EV;
        ratio_delta_V_EV_1_SQUAREM=ratio_delta;
        iter_info_V_EV_1_SQUAREM=iter_info;
    elseif tune_param>1
        results_V_EV_2_SQUAREM(m,:)=results_V_EV;
        ratio_delta_V_EV_2_SQUAREM=ratio_delta;
        iter_info_V_EV_2_SQUAREM=iter_info;
    end% tune_param==0 or 1 or others
end


end % method
