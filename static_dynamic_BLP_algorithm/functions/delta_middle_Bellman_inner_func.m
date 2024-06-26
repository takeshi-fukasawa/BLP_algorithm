function [out,other_vars]=...
    delta_middle_Bellman_inner_func(...
    delta_initial,V_initial0,mu_ijt,S_j_data,beta_C,...
    rho,weight,tune_param_BLP,spec_default,...
    weight_V,x_V,hot_start_V_spec)

global feval_Bellman V_initial_hot_start

spec=spec_default;

if hot_start_V_spec==1
    if feval_Bellman==0
        V_initial_hot_start=V_initial0;
    end
else
    V_initial_hot_start=V_initial0;
end


    [output_spectral,other_vars_Bellman,iter_info_Bellman]=...
        spectral_func(@Bellman_update_func,spec,{V_initial_hot_start},...
        delta_initial,mu_ijt,beta_C,rho,...
        weight_V,x_V);

    V_sol=output_spectral{1};

    %%%%%%%%%%%%%%

    IV=other_vars_Bellman.IV;

    T=size(mu_ijt,4);
    n_dim_V=size(V_sol,4);

    denom=exp(sum(IV(:,:,:,1:T),3,'omitnan'))+...
        exp(other_vars_Bellman.v_i0t_tilde(:,:,:,1:T));
    s_i0t_ccp=exp(other_vars_Bellman.v_i0t_tilde(:,:,:,1:T))./denom;

    %%%%%%%%%%%%

    feval_Bellman=feval_Bellman+iter_info_Bellman.feval;


    IV=other_vars_Bellman.IV;
    
    if T==n_dim_V
       s_igt_ccp=exp(IV-V_sol);
       s_ijt_given_g_ccp=...
           other_vars_Bellman.numer_1./other_vars_Bellman.denom_1;
    else
       s_igt_ccp=exp(IV(:,:,:,1:T)-V_sol(:,:,:,1:T));
       s_ijt_given_g_ccp=...
           other_vars_Bellman.numer_1(:,:,:,1:T)./...
           other_vars_Bellman.denom_1(:,:,:,1:T);
    end
 
    s_ijt_ccp=s_ijt_given_g_ccp.*s_igt_ccp;

    s_jt_predict=compute_s_jt_func(s_ijt_ccp,s_i0t_ccp,weight);

    delta_updated=delta_update_func(...
        delta_initial,S_j_data,s_jt_predict,rho,tune_param_BLP);

    out={delta_updated};

    other_vars.V_sol=V_sol;
    other_vars.s_jt_predict=s_jt_predict;
    
    V_initial_hot_start=V_sol;

end
