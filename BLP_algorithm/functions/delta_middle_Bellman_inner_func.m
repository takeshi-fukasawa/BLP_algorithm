function [out,other_vars]=...
    delta_middle_Bellman_inner_func(...
    delta_initial,V_initial0,mu_ij,S_j_data,beta_C,L,...
    rho,weight,spectral_V_spec,tune_param_BLP,spec_default,...
    weight_V,x_V,hot_start_V_spec)

global feval_Bellman V_initial_hot_start

spec=spec_default;
if spectral_V_spec==1
    spec.update_spec=4;
else
    spec.update_spec=0;
end

if hot_start_V_spec==1
    if feval_Bellman==0
        V_initial_hot_start=V_initial0;
    end
else
    V_initial_hot_start=V_initial0;
end


    [output_spectral,other_vars_Bellman,iter_info_Bellman]=...
        spectral_func(@Bellman_update_func,spec,{V_initial_hot_start},...
        delta_initial,mu_ij,beta_C,rho,...
        weight_V,x_V);

    V_sol=output_spectral{1};

    feval_Bellman=feval_Bellman+iter_info_Bellman.feval;

    T=size(mu_ijt,4);
    n_dim_V=size(V_sol,4);

    IV=other_vars_Bellman.IV;
    
    if T==n_dim_V
       s_igt_ccp=exp(IV-V_sol);
       s_ijt_given_g_ccp=...
           other_vars_Bellman.numer_1./other_vars_Bellman.denom_1;
    else
       s_igt_ccp=exp(IV(:,:,:,T+1:end)-V_sol(:,:,:,T+1:end));
       s_ijt_given_g_ccp=...
           other_vars_Bellman.numer_1(:,:,:,:,T+1:end)./...
           other_vars_Bellman.denom_1(:,:,:,:,T+1:end);
    end
 
    s_ijt_ccp=s_ijt_given_g_ccp.*s_igt_ccp;

    s_jt_predict=compute_s_jt_func(s_ijt_ccp,weight);

    delta_updated=delta_update_func(...
        delta_initial,S_j_data,s_jt_predict,rho,tune_param_BLP);

    out={delta_updated};

    other_vars.V_sol=V_sol;
    other_vars.s_jt_predict=s_jt_predict;
    
    V_initial_hot_start=V_sol;

end
