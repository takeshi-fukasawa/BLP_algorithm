function [output,other_vars]=...
    BLP_Bellman_joint_update_func(...
    delta_initial,V_initial,weight,mu_ijt,rho,...
    S_j_data,weight_V,x_V,beta_C,L,tune_param_BLP)

    [resid_V,other_vars]=Bellman_update_func(...
    V_initial,delta_initial,mu_ijt,beta_C,rho,...
    weight_V,x_V);

    V_updated=V_initial-resid_V{1};
    
    IV=other_vars.IV;
    s_igt_ccp=exp(IV(:,:,:,:,1)-V_updated(:,:,:,:,1));

    s_ijt_given_g_ccp=...
    other_vars.numer_1(:,:,:,:,1)./other_vars.denom_1(:,:,:,:,1);
    s_ijt_ccp=s_ijt_given_g_ccp.*s_igt_ccp;

    s_jt_predict=sum(s_ijt_ccp(:,:,:,:,1).*weight,2);
    

    delta_updated=delta_update_func(...
        delta_initial,S_j_data,s_jt_predict,rho,tune_param_BLP);


resid_delta=delta_initial-delta_updated;
resid_V=V_initial-V_updated;
%%%%%%%%%%%%%%%%%


    output={resid_delta,resid_V};
    other_vars.IV=IV;
    other_vars.s_jt_predict=s_jt_predict;
    other_vars.resid_delta=resid_delta;
    
end
