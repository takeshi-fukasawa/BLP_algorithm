function [output,other_vars]=...
    BLP_Bellman_joint_update_func(...
    delta_initial,V_initial,weight,mu_ij,rho,...
    S_j_data,weight_V,x_V,beta_C,L,tune_param_BLP)

u_ij_tilde=delta_initial+mu_ij+(beta_C^L).*V_initial;%J*I*G

EV=compute_EV_func(V_initial,[],weight_V,x_V);
u_i0_tilde=beta_C*EV;%J*I*G

    [s_j_predict,ChoiceProb,s_ij_given_g_ccp,s_ig_ccp,...
    numer_1,denom_1,numer_2,denom_2]=...
    share_func(u_ij_tilde,u_i0_tilde,rho,weight);

V_updated=log(denom_2);

speed_param=1-rho; %speed_param=1;
delta_updated=delta_initial+speed_param.*(log(S_j_data)-log(s_j_predict));%J*1*G 

delta_updated=delta_initial+speed_param.*(log(S_j_data)-log(s_j_predict));%J*1*G
  
  if tune_param_BLP>0
      
    s_0_predict=1-sum(s_j_predict,1);
    S_0_data=1-sum(S_j_data);
    delta_updated=delta_updated-tune_param_BLP*(log(S_0_data)-log(s_0_predict));
  end

resid_delta=delta_initial-delta_updated;
resid_V=V_initial-V_updated;
%%%%%%%%%%%%%%%%%


    output={resid_delta,resid_V};
    other_vars=[];
end
