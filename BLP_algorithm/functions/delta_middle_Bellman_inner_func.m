function [out,other_vars]=...
    delta_middle_Bellman_inner_func(...
    delta_initial,V_initial0,mu_ij,S_j_data,beta_C,L,...
    rho,weight,spectral_V_spec,tune_param_BLP,spec_default,...
    weight_V,x_V)

global feval_Bellman

spec=spec_default;
if spectral_V_spec==1
    spec.update_spec=4;
else
    spec.update_spec=0;
end

        [output_spectral,other_vars_Bellman,iter_info_Bellman]=...
        spectral_func(@Bellman_update_func,spec,{V_initial0},...
        delta_initial,mu_ij,beta_C,rho,...
        weight_V,x_V);

    V_sol=output_spectral{1};

    feval_Bellman=feval_Bellman+iter_info_Bellman.feval;

    IV=other_vars_Bellman.IV;
    s_igt_ccp=exp(IV(:,:,:,:,1)-V_sol(:,:,:,:,1));
    
    s_ijt_given_g_ccp=...
    other_vars_Bellman.numer_1(:,:,:,:,1)./other_vars_Bellman.denom_1(:,:,:,:,1);
    s_ijt_ccp=s_ijt_given_g_ccp.*s_igt_ccp;

    s_jt_predict=compute_s_jt_func(s_ijt_ccp,weight);

    delta_updated=delta_update_func(...
        delta_initial,S_j_data,s_jt_predict,rho,tune_param_BLP);

    out={delta_updated};

    other_vars.V_sol=V_sol;
    other_vars.s_jt_predict=s_jt_predict;
    
end
