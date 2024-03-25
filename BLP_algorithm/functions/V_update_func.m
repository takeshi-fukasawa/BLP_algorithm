function [output,other_vars]=...
    V_update_func(...
    V_initial,weight,mu_ijt,...
    S_jt_data,S_0t_data,...
    weight_V,x_V,beta_C,tune_param,Newton_spec)
 
    %%% G==1 & rho==0 case only
    [J,ns,~,T]=size(mu_ijt);
    n_dim_V=size(V_initial,4);

    if n_dim_V>T %% IVS spec
       V_initial_obs_pt=V_initial(:,:,:,1:T);%1*ns*1*T; V at obs data pts
    else
       n_gird_IV=n_dim_V-T;
       V_initial_obs_pt=V_initial;%1*ns*1*T
    end

    %% Compute delta

    %%%%%%%%%%%%%%%%%%%
    spec_1=1;
    if spec_1==1 & T==n_dim_V
        EV=compute_EV_func(...
                V_initial,[],weight_V,x_V,T);
    
        v_i0t_tilde=beta_C*EV;
        s_i0t_ccp=reshape(exp(v_i0t_tilde-V_initial),1,ns,1,T);
    else
        s_i0t_ccp=[];
    end
    %%%%%%%%%%%%%%%%%%%%%%
 
    %% Compute delta
    [delta_jt,Pr0]=compute_delta_from_V_func(...
       mu_ijt,weight,S_jt_data,V_initial_obs_pt,...
       [],s_i0t_ccp);

    if spec_1==1 & T==n_dim_V
        s_0t_predict=...
            sum(reshape(weight,1,ns,1,1).*(1-Pr0+Pr0.*s_i0t_ccp),[2,5]);%1*1*1*T
        delta_jt=delta_jt-tune_param*(log(S_0t_data)-log(s_0t_predict));
    end 

    %% Update V
    rho=0;
    v_ijt_tilde=delta_jt+mu_ijt;

    [numer_1,denom_1,IV_new,EV]=compute_IV_EV_func(...
            V_initial,v_ijt_tilde,beta_C,rho,weight_V,x_V);

    v_i0t_tilde=beta_C*EV;

    if n_dim_V==T
        s_i0t_ccp=reshape(exp(v_i0t_tilde-V_initial),1,ns,1,T);
    else
        s_i0t_ccp=reshape(exp(v_i0t_tilde(:,:,:,1:T)-V_initial(:,:,:,1:T)),1,ns,1,T);
    end

    s_0t_predict=...
        sum(reshape(weight,1,ns,1,1).*(1-Pr0+Pr0.*s_i0t_ccp),[2,5]);%1*1*1*T 
    S_0_ratio=s_0t_predict./reshape(S_0t_data,1,1,1,T);%1*1*1*T

    if n_dim_V==T & spec_1==1
        V_updated=log(exp(v_i0t_tilde)+exp(IV_new));%1*ns*1*n_dim_V
    elseif n_dim_V==T & spec_1==0
        V_updated=log(exp(v_i0t_tilde)+exp(IV_new).*(S_0_ratio.^tune_param));%1*ns*1*n_dim_V
    else
        V_obs_pt=log(exp(v_i0t_tilde(:,:,:,1:T))+...
            exp(IV_new(:,:,:,1:T).*(S_0_ratio.^tune_param)));%1*ns*1*T
        V_grid=log(exp(v_i0t_tilde(:,:,:,T+1:end))+...
            exp(IV_new(:,:,:,T+1:end)));%1*ns*1*n_grid_IV
        V_updated=cat(4,V_obs_pt,V_grid);%1*ns*1*n_dim_V
    end

%%%%%%%%%%%%%%%%%%%%%%%
    if Newton_spec==1 & n_grid_V==1
        %%% Newton iteration
        prob_0=reshape(exp(v_i0t_tilde),1,ns,1,T)./reshape(exp(V_initial),1,ns,1,T);%1*ns*1*T
        diff=log(exp(v_i0t_tilde)+temp_1.*(temp_2))-V_initial;
        V_updated=V_initial+diff./(1-beta_C*prob_0);%1*ns*1*T
    end
    %%%%%%%%%%%%%%%%%%%%%%%%

    output={V_updated};
    other_vars.s_i0t_ccp=s_i0t_ccp;
    other_vars.v_i0t_tilde=v_i0t_tilde;
    other_vars.IV=IV_new;
    other_vars.delta_jt=delta_jt;
    other_vars.v_ijt_tilde=v_ijt_tilde;
    
end
