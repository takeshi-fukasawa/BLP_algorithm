
%% V_update_func 
for method=1:4
    run spec_settings.m

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
        [s_jt_predict,s_ijt_ccp]=...
            share_func(v_ijt_tilde,...
            v_i0t_tilde,rho_est,weight);
    else
        [s_jt_predict,s_ijt_ccp]=...
            share_func(v_ijt_tilde(:,:,:,1:T),...
            v_i0t_tilde(:,:,:,1:T),rho_est,weight);
    end
    results_V=results_output_func(iter_info,s_jt_predict,S_jt_data);
    ratio_delta=delta_sol./delta_jt_true;


if method==1
    if tune_param==0
        results_V_0(m,:)=results_V;
        ratio_delta_V_0=ratio_delta;
        iter_info_V_0=iter_info;
    elseif tune_param==1
        results_V_1(m,:)=results_V;
        ratio_delta_V_1=ratio_delta;
        iter_info_V_1=iter_info;
    elseif tune_param>1
        results_V_2(m,:)=results_V;
        ratio_delta_V_2=ratio_delta;
        iter_info_V_2=iter_info;
    end

elseif method==2
    if tune_param==0
        results_V_0_spectral(m,:)=results_V;
        ratio_delta_V_0_spectral=ratio_delta;
        iter_info_V_0_spectral=iter_info;
    elseif tune_param==1
        results_V_1_spectral(m,:)=results_V;
        ratio_delta_V_1_spectral=ratio_delta;
        iter_info_V_1_spectral=iter_info;
    elseif tune_param>1
        results_V_2_spectral(m,:)=results_V;
        ratio_delta_V_2_spectral=ratio_delta;
        iter_info_V_2_spectral=iter_info;
    end% tune_param==0 or 1 or others
elseif method==3
    if tune_param==0
        results_V_0_SQUAREM(m,:)=results_V;
        ratio_delta_V_0_SQUAREM=ratio_delta;
        iter_info_V_0_SQUAREM=iter_info;
    elseif tune_param==1
        results_V_1_SQUAREM(m,:)=results_V;
        ratio_delta_V_1_SQUAREM=ratio_delta;
        iter_info_V_1_SQUAREM=iter_info;
    elseif tune_param>1
        results_V_2_SQUAREM(m,:)=results_V;
        ratio_delta_V_2_SQUAREM=ratio_delta;
        iter_info_V_2_SQUAREM=iter_info;
    end% tune_param==0 or 1 or others
elseif method==4
    if tune_param==0
        results_V_0_Anderson(m,:)=results_V;
        ratio_delta_V_0_Anderson=ratio_delta;
        iter_info_V_0_Anderson=iter_info;
    elseif tune_param==1
        results_V_1_Anderson(m,:)=results_V;
        ratio_delta_V_1_Anderson=ratio_delta;
        iter_info_V_1_Anderson=iter_info;
    elseif tune_param>1
        results_V_2_Anderson(m,:)=results_V;
        ratio_delta_V_2_Anderson=ratio_delta;
        iter_info_V_2_Anderson=iter_info;
    end% tune_param==0 or 1 or others

end


end % method
