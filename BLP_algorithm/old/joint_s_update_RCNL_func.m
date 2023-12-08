
function output=...
    joint_s_update_RCNL_func(...
    log_s_i0_ccp_initial,log_s_ig_ccp_initial,weight,mu_ij,rho,...
    S_j_data,S_g_data,S_0_data,numer_1_without_delta)
 
global I_ig

 [J,ns]=size(mu_ij);
log_weight=log(weight);

%% Normalization of s_i0, s_ig

%%%%% Avoid the use of s_i0_ccp %%%%%%
%s_i0_ccp_initial=1-sum(exp(log_s_ig_ccp_initial),3);
%log_s_i0_ccp_initial=log(s_i0_ccp_initial);
%%%%%%%%%%%%%%%%%%%%%%%%%

s_i0_temp=exp(log_s_i0_ccp_initial+log_weight);
log_si0=log(s_i0_temp.*S_0_data./sum(s_i0_temp,2));%1*ns*1

s_i0_initial=exp(log_si0);
s_i0_ccp_initial=s_i0_initial./weight;
log_s_i0_ccp_initial=log_si0-log_weight;

%%%%%%
s_ig_temp=exp(log_s_ig_ccp_initial+log_weight);
log_sig=log(S_g_data.*s_ig_temp./sum(s_ig_temp,2));%1*ns*1

s_ig_initial=exp(log_sig);
s_ig_ccp_initial=s_ig_initial./weight;
log_s_ig_ccp_initial=log_sig-log_weight;

%%%%%%%

    min_val=1e-14;

    %%% Old
    %exp_I_ig=exp(log(max(min_val,1-s_i0_ccp_initial))-log(s_i0_ccp_initial));%1*ns*G
    
    %%% New
    %exp_I_ig=exp(log(s_ig_ccp_initial)-log(s_i0_ccp_initial));%1*ns*G
    
    I_ig=log((max(min_val,1-s_i0_ccp_initial)).*s_ig_ccp_initial./sum(s_ig_ccp_initial,3))-...
        log_s_i0_ccp_initial;

    % Bayes
    temp=numer_1_without_delta.*s_i0_initial.*exp(-rho./(1-rho).*I_ig);%J*ns*G

    posterior_i_j=temp./sum(temp,2);%J*ns*G

    s_ig_predict=sum(S_j_data.*posterior_i_j,1);%1*ns*G

    %temp=s_i0_initial;%1*ns
    %posterior_i_0=temp./sum(temp,2);%1*ns
    %s_i0_predict=S_0_data.*posterior_i_0;%1*ns


    % Either specification OK
    sum_s_ig_predict=sum(s_ig_predict,3);
    weight_predict=sum_s_ig_predict+s_i0_initial;%1*ns*1
   
    %weight_predict=sum_s_ig_predict+s_i0_predict;%1*ns
    
    difference=log_weight-log(weight_predict);

    speed_param=1-rho;
    %speed_param=1;
    
    log_s_i0_ccp_updated=log_s_i0_ccp_initial+speed_param.*difference;%J*1
    
    s_i0_ccp_updated=exp(log_s_i0_ccp_updated);%J*1


    s_ig_ccp_updated=(max(min_val,1-s_i0_ccp_updated)).*s_ig_predict./sum_s_ig_predict;%1*ns*G
    %%s_ig_ccp_updated=s_ig_predict./weight_predict;% a little slower than the specification above??
    %%s_ig_ccp_updated=s_ig_predict./(s_i0_ccp_updated.*weight+s_ig_predict);%%% Diverge when G>=2??
    
    log_s_ig_ccp_updated=log(s_ig_ccp_updated);%1*ns*G

    resid_log_s_i0_ccp=log_s_i0_ccp_initial-log_s_i0_ccp_updated;
    resid_log_s_ig_ccp=(log_s_ig_ccp_initial-log_s_ig_ccp_updated)*(1-rho);
    

    output={resid_log_s_i0_ccp,resid_log_s_ig_ccp};
end
