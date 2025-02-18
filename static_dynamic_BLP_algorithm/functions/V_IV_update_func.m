function [output,other_vars]=...
    V_IV_update_func(...
    V_initial,IV_initial,weight,mu_ijt,rho,...
    S_jt_data,S_0t_data,weight_V,x_V,beta_C,tune_param)
    
    %%% Currently static only?? (Pr0 not introduced??)

    [J,ns,G,T]=size(mu_ijt);   
    n_dim_V=size(V_initial,4);
    
    if T==n_dim_V
        V_initial_obs_pt=V_initial;%1*ns*1*T; V at obs data pts
        IV_initial_obs_pt=IV_initial;%1*ns*G*T; IV at obs data pts
    else
        V_initial_obs_pt=V_initial(:,:,:,1:T);%1*ns*1*T; V at obs data pts
        IV_initial_obs_pt=IV_initial(:,:,:,1:T);%1*ns*G*T; IV at obs data pts
    end
    
    if rho==0
    numer_1=exp(mu_ijt./(1-rho));%J*ns*G*T
    denom_1=sum(reshape(weight,1,ns).*...
        reshape(numer_1,J,ns,G,T).*...
        reshape(exp(-V_initial_obs_pt),1,ns,1,T),2);%J*1*G*T
    end


    %% Update IV
    
    [delta,s_ijt_up_to_scale,s_igt]=...
        compute_delta_from_V_IV_func(...
    mu_ijt,weight,S_jt_data,rho,V_initial_obs_pt,IV_initial_obs_pt);


    if rho==0
        IV_new=log(sum(S_jt_data.*numer_1./denom_1,1));%1*ns*G*T
    else
        v_ijt=reshape(delta,J,1,G,T)+reshape(mu_ijt,J,ns,G,T);%J*ns*G*T
        IV_new=(1-rho).*log(sum(exp(v_ijt/(1-rho)),1));%1*ns*G*T
    end

    if (n_dim_V-T)>=1
        %%% Inclusive value sufficiency
        IV_new=IVS_compute_IV_func(IV_new,n_dim_V-1);%1*ns*G*n_dim_V
    end

    weight_V=[];
    EV=compute_EV_func(V_initial,IV_new,weight_V,x_V,T);%1*ns*1*n_dim_V
    v_i0t=beta_C*EV;%1*ns*1*n_dim_V

    if n_dim_V==T
        v_i0t_obs_pt=v_i0t;%1*ns*1*T; v_i0t at obs data pts
    else
        v_i0t_obs_pt=v_i0t(:,:,:,1:T);%1*ns*1*T; v_i0t at obs data pts
    end

    if tune_param~=0
        s_0t_predict=sum(reshape(weight,1,ns).*...
            reshape(exp(v_i0t_obs_pt-V_initial_obs_pt),1,ns,1,T),2);%1*1*1*T
        s_0_ratio=s_0t_predict./reshape(S_0t_data,1,1,1,T);%1*1*1*T
    end

    if rho==0
        if tune_param==0
            IV_updated=IV_new;%1*ns*G*T
        else
            IV_updated=IV_new+tune_param.*log(s_0_ratio);%1*ns*G*T
        end

   else
        S_gt_data=sum(S_jt_data,1);%1*1*G*T
        s_gt_predict=sum(reshape(weight,1,ns,1).*...
            reshape(s_igt,1,ns,G,T),2);%1*1*G*T

        if tune_param==0
            IV_updated=IV_new+...
                rho*(log(S_gt_data)-log(s_gt_predict));%1*ns*G*T 
        else
            IV_updated=IV_new+...
                rho*(log(S_gt_data)-log(s_gt_predict))+...
                tune_param.*log(s_0_ratio);%1*ns*G*T  
        end     
    end


    %% Update V
%% Use IV_updated => faster
    V_updated=log(exp(v_i0t)+sum(exp(IV_updated),3));
    
    output={V_updated,IV_updated};
    other_vars=[];
end
