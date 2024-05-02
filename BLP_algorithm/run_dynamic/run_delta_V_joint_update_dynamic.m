

%% BLP_Bellman_joint_update_func 
for method=1:3
    run spec_settings.m

    [output_spectral,other_vars,iter_info]=...
        spectral_func(@BLP_Bellman_joint_update_func,spec,...
        {delta_initial0,V_initial0},...
        weight,mu_ijt_est,rho_est,...
    S_jt_data,weight_V,x_V,beta_C,tune_param_BLP);
    

    delta_sol=output_spectral{1};
    V_sol=output_spectral{2};


results=results_output_func(iter_info,other_vars.s_jt_predict,S_jt_data);
ratio_delta=delta_sol./delta_jt_true;

if method==1

if tune_param_BLP==0
results_V_BLP_0(m,:)=results;
ratio_delta_V_BLP_0=ratio_delta;
elseif tune_param_BLP==1
results_V_BLP_1(m,:)=results;
ratio_delta_V_BLP_1=ratio_delta;
end% tune_param_BLP==0 or 1

elseif method==2
    if tune_param_BLP==0
        results_V_BLP_0_spectral(m,:)=results;
        ratio_delta_V_BLP_0_spectral=ratio_delta;
        iter_info_V_BLP_0_spectral=iter_info;
        elseif tune_param_BLP==1
        results_V_BLP_1_spectral(m,:)=results;
        ratio_delta_V_BLP_1_spectral=ratio_delta;
        iter_info_V_BLP_1_spectral=iter_info;
        end% tune_param_BLP==0 or 1
elseif method==3
    if tune_param_BLP==0
        results_V_BLP_0_SQUAREM(m,:)=results;
        ratio_delta_V_BLP_0_SQUAREM=ratio_delta;
        iter_info_V_BLP_0_SQUAREM=iter_info;
        elseif tune_param_BLP==1
        results_V_BLP_1_SQUAREM(m,:)=results;
        ratio_delta_V_BLP_1_SQUAREM=ratio_delta;
        iter_info_V_BLP_1_SQUAREM=iter_info;
        end% tune_param_BLP==0 or 1

end
end % method
