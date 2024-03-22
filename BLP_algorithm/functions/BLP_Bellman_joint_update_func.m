function [output,other_vars]=...
    BLP_Bellman_joint_update_func(...
    delta_initial,V_initial,weight,mu_ijt,rho,...
    S_j_data,weight_V,x_V,beta_C,L,tune_param_BLP)

    [out,other_vars]=Bellman_update_func(...
    V_initial,delta_initial,mu_ijt,beta_C,rho,...
    weight_V,x_V);

    V_updated=out{1};
    
    IV=other_vars.IV;

    T=size(mu_ijt,4);
    T_temp=size(V_initial,4);
    if T_temp==T
        s_igt_ccp=exp(IV-V_initial);
        s_ijt_given_g_ccp=...
            other_vars.numer_1./other_vars.denom_1;
            s_ijt_ccp=s_ijt_given_g_ccp.*s_igt_ccp;
            s_jt_predict=compute_s_jt_func(s_ijt_ccp,weight);
    else
        s_igt_ccp=exp(IV(:,:,:,1:T)-V_initial(:,:,:,1:T));
        s_ijt_given_g_ccp=...
            other_vars.numer_1(:,:,:,1:T)./other_vars.denom_1(:,:,:,1:T);
            s_ijt_ccp=s_ijt_given_g_ccp.*s_igt_ccp;
            s_jt_predict=compute_s_jt_func(s_ijt_ccp,weight);
    end        
    delta_updated=delta_update_func(...
        delta_initial,S_j_data,s_jt_predict,rho,tune_param_BLP);

    output={delta_updated,V_updated};
    other_vars.IV=IV;
    other_vars.s_jt_predict=s_jt_predict;

end
