function [output,other_vars]=...
    BLP_update_func(...
    delta_initial,weight,mu_ijt,rho,...
    S_jt_data,tune_param_BLP)
 
  [J,ns,G,T]=size(mu_ijt);
  u_ijt_tilde=delta_initial+mu_ijt;%J*I*G*T
  u_i0t_tilde=zeros(1,ns,1,T);
  [s_jt_predict,ChoiceProb,s_ijt_given_g_ccp,s_igt_ccp]=...
  share_func(u_ijt_tilde,u_i0t_tilde,rho,weight);%J*1*G*T

  % contraction mapping:
delta_updated=delta_update_func(delta_initial,...
S_jt_data,s_jt_predict,rho,tune_param_BLP);

    output={delta_updated};

    other_vars=[];
    
end
