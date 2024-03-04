
function output=...
    joint_update_func(...
    r_initial,Inc_initial,weight,mu_ij,rho,...
    S_j_data,S_0_data,numer_1_without_delta)
 
global denom_1 numer_2 denom_2

    [J,ns]=size(mu_ij);

    %%% r: log(s_i0)
    %log_si0_new=r_normalize_func((r_initial-Inc_initial*rho/(1-rho))+log(weight),S_0_data);
    %%% r: log(s_i0_ccp)
    log_si0_new=r_normalize_func((r_initial)+log(weight),S_0_data);

    log_si0=log_si0_new;
    
    s_i0_initial=exp(log_si0);
    s_i0_ccp_initial=s_i0_initial./weight;

    min_val=0.00001;
    numer_2=exp(log(max(min_val,1-s_i0_ccp_initial))-log(s_i0_ccp_initial));%1*ns
    %numer_2=exp(Inc_initial);

    %%%% Below: Up to scale convergence???
    %numer_2=exp(log(exp(Inc_initial)./(1+exp(Inc_initial)))-log(s_i0_ccp_initial));%1*ns; 
    
    %numer_2=exp(log(max(min_val,1-s_i0_ccp_initial)))-...
    %    log(1./(1+exp(Inc_initial)));%%% Diverge
    
    denom_1=exp(log(numer_2)./(1-rho));%1*ns
    denom_2=1+numer_2;

    % Bayes
    temp=numer_1_without_delta.*s_i0_initial.*numer_2./denom_1;%J*ns
    %%%%temp=numer_1_without_delta.*s_i0_initial.*exp(Inc_initial.*rho/(1-rho)); %% Diverge
    posterior_i_j=temp./sum(temp,2);%J*ns

    s_ig_predict=sum(S_j_data.*posterior_i_j,1);%1*ns

    temp=s_i0_initial;%1*ns
    posterior_i_0=temp./sum(temp,2);%1*ns

    s_i0_predict=S_0_data.*posterior_i_0;%1*ns

    % Either specification OK
    weight_predict=s_ig_predict+s_i0_initial;%1*ns
   
    %weight_predict=s_ig_predict+s_i0_predict;%1*ns
    
    difference=log(weight)-log(weight_predict);

    %%% Not working...
    %difference=log(weight.*(1-s_i0_ccp_initial))-...
        log(s_ig_predict);

    speed_param=1-rho;
    speed_param=(1-rho);
    %speed_param=1;
    
    r_updated=r_initial+speed_param.*difference;%J*1
    
    resid_r=r_initial-r_updated;

    Inc_updated=numer_2;
    Inc_updated=log(1-exp(r_updated))-r_updated;
    %Inc_updated=log(exp(Inc_initial)./(1+exp(Inc_initial)))-r_updated;
    difference_Inc=log(s_i0_ccp_initial)+log(1+exp(Inc_initial));
    %Inc_updated=Inc_initial-0.5*difference_Inc;
    
    resid_Inc=Inc_initial-Inc_updated;

    output={resid_r,resid_Inc};
end
