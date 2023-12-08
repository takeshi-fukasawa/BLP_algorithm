function output=...
    IV_update_func(...
    V_initial,IV_initial,weight,mu_ijt,rho,...
    S_jt_data,S_0t_data,beta_C,L,tune_param)
    
    [J,ns,G,T]=size(mu_ijt);   

    
    if rho==0
    numer_1=exp(mu_ijt./(1-rho));%J*ns*G*T
    denom_1=sum(reshape(weight,1,ns).*...
        reshape(numer_1,J,ns,G,T).*...
        reshape(exp(-V_initial),1,ns,1,T),2);%J*1*G*T
    end


    %% Update IV
    
    [delta,s_ijt_up_to_scale,s_igt]=compute_delta_from_V_IV_func(...
    mu_ijt,weight,S_jt_data,rho,V_initial,IV_initial,beta_C,L);

      weight_V=[];
       EV=compute_EV_func(V_initial,weight_V);

    s_0t_predict=sum(reshape(weight,1,ns).*...
        reshape(exp(beta_C*EV-V_initial),1,ns,1,T),2);%1*1*1*T
    s_0_ratio=s_0t_predict./reshape(S_0t_data,1,1,1,T);%1*1*1*T


    if rho==0

        IV_new=log(sum(S_jt_data.*numer_1./denom_1,1));%1*ns*G*T%%%%%%

        IV_updated=IV_new+tune_param.*log(s_0_ratio);%1*ns*G*T

    else
        v_ijt=reshape(delta,J,1,G,T)+reshape(mu_ijt,J,ns,G,T);%J*ns*G*T
        IV_new=(1-rho).*log(sum(exp(v_ijt/(1-rho)),1));%1*ns*G*T

        %s_ijt=reshape(S_jt_data,J,1,G,T).*...
        %s_ijt_up_to_scale./sum(reshape(weight,1,ns,1).*s_ijt_up_to_scale,2);%J*ns*G*T
        %s_igt=sum(s_ijt,1);%1*ns*G*T

        S_gt_data=sum(S_jt_data,1);%1*1*G*T
        s_gt_predict=sum(reshape(weight,1,ns,1).*...
            reshape(s_igt,1,ns,G,T),2);%1*1*G*T
        IV_updated=IV_new+...
            rho*(log(S_gt_data)-log(s_gt_predict))+...
            tune_param.*log(s_0_ratio);%1*ns*G*T
        
        temp=IV_updated-beta_C*V_initial;
       
    end

    resid_IV=IV_initial-IV_updated;

    %% Update V
    v_i0t=beta_C*V_initial;%1*ns*1*T
    V_updated=log(exp(v_i0t)+sum(exp(IV_updated),3));%% Use IV_updated => faster
    resid_V=V_initial-V_updated;

    output={resid_V,resid_IV};
end
