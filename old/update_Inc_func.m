%%% Slow convergence??
function output=...
    update_Inc_func(...
    Inc_initial,weight,mu_ij,rho,...
    S_j_data,S_0_data,numer_1_without_delta)

global denom_1 numer_2 denom_2
global posterior_i_j posterior_i_0

    [J,ns]=size(mu_ij);

    Inc_i0=0;
    %%%%% Normalization %%%%%
    S_g_data=1-S_0_data;
    numer_2=exp(Inc_initial);
    denom_2=numer_2+1;
    S_0_predict=sum(weight.*1./denom_2,2);
    RHS=(1+exp(Inc_initial)).*S_0_predict./S_0_data;
    Inc_initial_new=log(RHS-1);
    
    %%%%%%%%%
    S_g_predict=sum(weight.*numer_2./denom_2,2);
    %S_g_predict=1-S_0_predict;

    RHS=(1+1./exp(Inc_initial)).*S_g_predict./S_g_data;

    %Inc_initial_new=-log(RHS-1);
    %Inc_initial_new=Inc_initial+log(S_g_data)-log(S_0_data)...
    %    -log(S_g_predict)+log(S_0_predict);

    Inc_initial=Inc_initial_new;
    %%%%%%

    numer_2=exp(Inc_initial);
    denom_2=numer_2+1;

    denom_1=exp(Inc_initial./(1-rho));%1*ns

    s_i0_ccp=1./denom_2;%1*ns

    % Bayes
    temp=weight.*numer_1_without_delta.*s_i0_ccp.*numer_2./denom_1;%J*ns
    posterior_i_j=temp./sum(temp,2);%J*ns

    temp=weight.*s_i0_ccp;%1*ns
    posterior_i_0=temp./sum(temp,2);%1*ns

    %%%%%%%%

    s_ig_predict=sum(S_j_data.*posterior_i_j,1);%1*ns

    s_i0_predict=S_0_data.*posterior_i_0;%1*ns

    weight_predict=s_ig_predict+s_i0_predict;%1*ns
    
    difference=-(log(posterior_i_0)+log(denom_2))+...
        (1-rho)*(log(weight)-log(weight_predict));

    difference=-(1-rho)*(log(weight)-log(weight_predict));

    %difference=difference-(log(posterior_i_0)+log(denom_2));

    Inc_updated=Inc_initial+difference;

    %%%%%%%
    %%%%%%%%%
    if 1==0
        denom_temp=sum(S_j_data.*posterior_i_j,1)+S_0_data*posterior_i_0;
        prob_i_given_0=S_0_data*posterior_i_0./denom_temp;

        prob_ig=numer_2./denom_2;
    
        prob_i0_updated=S_0_data.*prob_i_given_0./weight;
        Inc_updated=log(prob_ig)-log(prob_i0_updated);
    end
    %%%%%%%%

    if 1==0
        Inc_updated=log(sum(S_j_data.*posterior_i_j,1))-...
        log(S_0_data.*posterior_i_0);

        speed_param=1-rho; % necessary to stabilize the convergence
        Inc_updated=Inc_initial+speed_param*(Inc_updated-Inc_initial);
    end

    resid_Inc=Inc_initial-Inc_updated;

    output={resid_Inc};
end
