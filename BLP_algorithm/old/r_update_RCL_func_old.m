
function output=...
    r_update_RCL_func(...
    r_initial,weight,mu_ij,rho,...
    S_j_data,S_0_data,numer_1_without_delta,...
    switch_r_spec)
 
global I_ig

    [J,ns]=size(mu_ij);

    if 1==0
    r_initial_new=r_normalize_func(r_initial,S_0_data);%%%%%%

    r_initial=r_initial_new;
    
    s_i0_initial=exp(r_initial);
    s_i0_ccp_initial=s_i0_initial./weight;

    min_val=1e-14;
    I_ig=log(max(min_val,1-s_i0_ccp_initial))-log(s_i0_ccp_initial);%1*ns

    % Bayes
    temp=numer_1_without_delta.*s_i0_initial.*exp(-rho./(1-rho).*I_ig);%J*ns
    posterior_i_j=temp./sum(temp,2);%J*ns

    s_ig_predict=sum(S_j_data.*posterior_i_j,1);%1*ns

    %temp=s_i0_initial;%1*ns
    %posterior_i_0=temp./sum(temp,2);%1*ns

    %s_i0_predict=S_0_data.*posterior_i_0;%1*ns


    % Either specification OK
    weight_predict=s_ig_predict+s_i0_initial;%1*ns
   
    %weight_predict=s_ig_predict+s_i0_predict;%1*ns
    
    difference=log(weight)-log(weight_predict);

    %%% Not working...
    %difference=log(weight.*(1-s_i0_ccp_initial))-...
    %    log(s_ig_predict);

    speed_param=1-rho;
    
    %speed_param=1;
    
    r_updated=r_initial+speed_param.*difference;%J*1
    %r_updated=r_normalize_func(r_updated,S_0_data);%%%%%%%%
    
    else
    %%%%%%%%%%% Kalouptsidi (2012) base spec %%%%%%%
    numer=exp(mu_ij+reshape(r_initial,1,ns));%J*ns
    denom=sum(numer,2);%J*1
    temp1=sum(reshape(S_j_data,J,1).*numer./denom,1);%1*ns

    numer=exp(reshape(r_initial,1,ns));%1*ns
    denom=sum(numer,2);%1*1
    temp2=numer./denom;

    difference=log(weight)-log(temp1+temp2);
    r_updated=r_initial+difference;

    temp=S_0_data-sum(exp(r_initial(1:end-1)));


    if temp>0 & 1==0
        % Kalouptsidi base spec
        r_updated(end)=log(temp);
    elseif temp<=0 & 1==0
            %warning("NaN val")
            %%% Switch to the stable but slow algorithm
            r_tilde_initial=r_initial-r_initial(end);
            r_tilde_updated=r_tilde_initial+difference;
            r_tilde_updated(end)=0;
            r_updated=r_tilde_updated+log(S_0_data)-log(sum(exp(r_tilde_updated)));
            %%%%%

    elseif 1==1
        r_updated(end)=0;% Kalouptsidi conservative spec

    end
    end

    %%%%%%%%%%%%%%%%%%%%
    
    resid=r_initial-r_updated;

    output={resid};
end
