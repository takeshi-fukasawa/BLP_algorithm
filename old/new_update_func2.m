%%% Slow convergence??
function output=...
    r_update_RCL_func2(...
    V_initial,weight,mu_ij,rho,...
    S_j_data,S_0_data,numer_1_without_delta)
 
global denom_1 numer_2 denom_2

global s_i0_initial s_i0_ccp_initial s_ig_predict weight_predict difference

[J,ns]=size(mu_ij);

    %V_hat_initial=(log(exp(V_initial)-1))./(1-rho);

    r_temp=r_normalize_func(log(weight)-V_initial,S_0_data);%%%%%%%%
    V_initial_new=log(weight)-r_temp;

    denom_2=exp(V_initial_new);
    numer_2=denom_2-1;%%%%
    denom_1=exp(log(numer_2)./(1-rho));%1*ns

    s_i0_ccp_initial=1./denom_2;%1*ns
    s_i0_initial=weight.*s_i0_ccp_initial;%1*ns

    % Bayes
    temp=numer_1_without_delta.*s_i0_initial.*numer_2./denom_1;%J*ns
    posterior_i_j=temp./sum(temp,2);%J*ns

    s_ig_predict=sum(S_j_data.*posterior_i_j,1);%1*ns

    s_i0_ccp_updated=s_i0_initial./(s_i0_initial+s_ig_predict);

    weight_predict=s_ig_predict+s_i0_initial;%1*ns
    
    %difference=log(s_i0_ccp_updated)-log(s_i0_ccp_initial);
    
    difference=log(weight)-log(weight_predict);

    speed_param=1-rho;
    %speed_param=1;

    V_updated=V_initial_new-(1-rho).*difference;%J*1
    %V_hat_updated=(log(exp(V_updated)-1))./(1-rho);
    
    %resid_V_hat=V_hat_initial-V_hat_updated;
    resid_V=V_initial-V_updated;
    
    %output={resid_V_hat,resid_V};
    output={resid_V};

end
