function [output,other_vars]=...
    IV_update_func(...
    IV_initial,weight,mu_ijt,rho,...
    S_jt_data,S_0t_data,tune_param)
  
    %%% Static case only %%%
    %%% Allow for rho>0,G>=2 case %%%
  
    [J,ns,G,T]=size(mu_ijt);   
    V=log(1+sum(exp(IV_initial),3));%1*ns*1*T
     
    if rho==0
    numer_1=exp(mu_ijt./(1-rho));%J*ns*G*T
    denom_1=sum(reshape(weight,1,ns).*...
        reshape(numer_1,J,ns,G,T).*...
        reshape(exp(-V),1,ns,1,T),2);%J*1*G*T
    end


    %% Update IV
    
    [delta,s_ijt_up_to_scale,s_igt]=...
        compute_delta_from_V_IV_func(...
    mu_ijt,weight,S_jt_data,rho,[],IV_initial);


    if rho==0
        IV_new=log(sum(S_jt_data.*numer_1./denom_1,1));%1*ns*G*T
    else
        v_ijt=reshape(delta,J,1,G,T)+reshape(mu_ijt,J,ns,G,T);%J*ns*G*T
        IV_new=(1-rho).*log(sum(exp(v_ijt/(1-rho)),1));%1*ns*G*T

    end
    
    s_0t_predict=sum(reshape(weight,1,ns).*...
        reshape(exp(-V),1,ns,1,T),2);%1*1*1*T
    s_0_ratio=s_0t_predict./reshape(S_0t_data,1,1,1,T);%1*1*1*T

    if rho==0
        IV_updated=IV_new+tune_param.*log(s_0_ratio);%1*ns*G*T

   else
        S_gt_data=sum(S_jt_data,1);%1*1*G*T
        s_gt_predict=sum(reshape(weight,1,ns,1).*...
            reshape(s_igt,1,ns,G,T),2);%1*1*G*T
        IV_updated=IV_new+...
            rho*(log(S_gt_data)-log(s_gt_predict))+...
            tune_param.*log(s_0_ratio);%1*ns*G*T       
    end

    output={IV_updated};
    other_vars=[];
end
