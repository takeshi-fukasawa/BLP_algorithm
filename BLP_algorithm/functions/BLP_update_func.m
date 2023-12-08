function output=...
    BLP_update_func(...
    delta_initial,weight,mu_ijt,rho,...
    S_jt_data,tune_param_BLP)
 
  [J,ns,G,T]=size(mu_ijt);
  u_ijt_tilde=delta_initial+mu_ijt;%J*I*G*T
  u_i0t_tilde=zeros(1,ns,1,T);
  [s_jt_predict,ChoiceProb,s_ijt_given_g_ccp,s_igt_ccp]=...
  share_func(u_ijt_tilde,u_i0t_tilde,rho,weight);%J*1*G*T

  % contraction mapping:
  speed_param=1-rho;
  %speed_param=1;
  delta_updated=delta_initial+speed_param.*(log(S_jt_data)-log(s_jt_predict));%J*1*G
  
      
  if tune_param_BLP>0
    s_0t_predict=1-sum(s_jt_predict(:));
    S_0t_data=1-sum(S_jt_data(:));
    delta_updated=delta_updated-tune_param_BLP*(log(S_0t_data)-log(s_0t_predict));

     if rho>0
        S_gt_data=sum(S_jt_data,1);%1*1*G*T
        [~,ns,G,T]=size(s_igt_ccp);
        s_gt_predict=sum(reshape(weight,1,ns,1).*reshape(s_igt_ccp,1,ns,G,T),2);%1*1*G*T
        delta_updated=delta_updated+tune_param_BLP*rho*(log(S_gt_data)-log(s_gt_predict));

      end% rho>0

  end

  resid=delta_initial-delta_updated;

    output={resid};
end
