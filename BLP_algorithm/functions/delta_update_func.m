function delta_updated=delta_update_func(...
    delta_initial,S_jt_data,s_jt_predict,rho,tune_param_BLP)


    delta_updated=delta_initial+(1-rho).*(log(S_jt_data)-log(s_jt_predict));%J*1*G
  
      
    if tune_param_BLP>0
        s_0t_predict=1-sum(s_jt_predict,[1,3]);%1*1*1*T
        S_0t_data=1-sum(S_jt_data,[1,3]);%1*1*1*T
        delta_updated=delta_updated-tune_param_BLP*(log(S_0t_data)-log(s_0t_predict));

        if rho>0
            S_gt_data=sum(S_jt_data,1);%1*1*G*T
            s_gt_predict=sum(s_jt_predict,1);%1*1*G*T

            delta_updated=delta_updated+tune_param_BLP*rho*(log(S_gt_data)-log(s_gt_predict));

        end% rho>0

  end %tune_param_BLP>0
return
