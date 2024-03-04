
function output=...
    V_update_func(...
    V_initial,weight,mu_ij,rho,...
    S_j_data,S_0_data,numer_1_without_delta,...
    beta_C,L,tune_param)
 
global I_ig

    [J,ns]=size(mu_ij);

    s_i0_ccp_temp=exp((beta_C-1).*V_initial);
    s_i0_ccp_initial=S_0_data./sum(weight.*s_i0_ccp_temp,2).*s_i0_ccp_temp;
    s_i0_initial=s_i0_ccp_initial.*weight;
    V_initial=log(s_i0_ccp_initial)./(beta_C-1);

    min_val=1e-14;
    I_ig=log(max(min_val,1-s_i0_ccp_initial))+V_initial;%1*ns

    % Bayes
    temp=numer_1_without_delta.*exp((beta_C.^L).*V_initial./(1-rho)-beta_C.*V_initial).*...
        s_i0_initial.*exp(-rho./(1-rho).*I_ig);%J*ns
    posterior_i_j=temp./sum(temp,2);%J*ns

    s_ig_predict=sum(S_j_data.*posterior_i_j,1);%1*ns

    weight_predict=s_ig_predict+s_i0_initial;%1*ns
       
    difference=log(weight)-log(weight_predict);


    speed_param=(1-rho)*(1-beta_C);%%%%%
    speed_param=1-rho;%% Works well

    %speed_param=1;
    
    V_updated=log(exp(beta_C.*V_initial)+exp(I_ig))-speed_param.*difference;%J*1


    resid=V_initial-V_updated;

    output={resid};
end
