
function output=...
    V_update_func(...
    V_initial,weight,mu_ijt,rho,...
    S_jt_data,S_0t_data,numer_1_without_delta,...
    beta_C,L,tune_param,Newton_spec)
 
    %%% G==1 case only
    [J,ns,~,T]=size(mu_ijt);
    n_grid=size(V_initial,5);

   weight_V=[];
   EV=compute_EV_func(V_initial,weight_V);
    v_i0t=beta_C*EV;%1*ns*1*T
    
    numer_1=exp(mu_ijt);%J*ns*1*T
    denom_1=sum(reshape(weight,1,ns,1,1).*reshape(numer_1,J,ns,1,T).*...
        reshape(exp(-V_initial),1,ns,1,T),2);%J*1*1*T
    exp_IV_updated=sum(S_jt_data.*numer_1./denom_1,1);%1*ns*1*T
    
    numer_2=sum(reshape(weight,1,ns).*reshape(exp(v_i0t-V_initial),1,ns,1,T),2);%1*1*1*T
    S_0_ratio=numer_2./reshape(S_0t_data,1,1,1,T);%1*1*1*T


    V_updated=log(exp(v_i0t)+exp_IV_updated.*(S_0_ratio.^tune_param));%1*ns*1*T

    if Newton_spec==1
        %%% Newton iteration
        prob_0=reshape(exp(v_i0t),1,ns,1,T)./reshape(exp(V_initial),1,ns,1,T);%1*ns*1*T
        diff=log(exp(v_i0t)+temp_1.*(temp_2))-V_initial;
        V_updated=V_initial+diff./(1-beta_C*prob_0);%1*ns*1*T
    end
    
    resid=V_initial-V_updated;

    output={resid};
end
