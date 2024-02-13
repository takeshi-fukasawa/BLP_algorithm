function [output,other_vars]=...
    BLP_Bellman_joint_update_func(...
    delta_initial,V_initial,weight,mu_ij,rho,...
    S_j_data,weight_V,x_V,beta_C,L,tune_param_BLP)

    n_dim_V=size(V_initial,5);
    u_ij_tilde=delta_initial+mu_ij;%J*I*G*T*n_dim_V
    
    IV=(1-rho).*log(sum(exp(u_ij_tilde./(1-rho)),1));%1*ns*G*T*n_dim_V; temporary under IVS
    
    if n_dim_V>=2
    IV=IVS_compute_IV_func(IV(:,:,:,:,1),n_dim_V-1);%1*ns*1*T*n_dim_V
    end
    
    EV=compute_EV_func(V_initial,IV,weight_V,x_V);
    
    
    u_i0_tilde=beta_C*EV(:,:,:,:,1);%J*I*G*T

    [s_jt_predict,ChoiceProb,s_ij_given_g_ccp,s_ig_ccp,...
    numer_1,denom_1,numer_2,denom_2]=...
    share_func(u_ij_tilde,u_i0_tilde,rho,weight);

    IV_state=log(numer_2);
    V_updated=log(denom_2);

    delta_updated=delta_update_func(...
        delta_initial,S_j_data,s_jt_predict,rho,tune_param_BLP);


resid_delta=delta_initial-delta_updated;
resid_V=V_initial-V_updated;
%%%%%%%%%%%%%%%%%


    output={resid_delta,resid_V};
    other_vars.IV=IV_state;
    other_vars.s_jt_predict=s_jt_predict;
    other_vars.resid_delta=resid_delta;
    
end
