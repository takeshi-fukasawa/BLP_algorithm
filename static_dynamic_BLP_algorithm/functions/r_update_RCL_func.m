
function [output,other_vars]=...
    r_update_RCL_func(...
    r_initial,weight,mu_ij,rho,...
    S_j_data,S_0_data,numer_1_without_delta,...
    switch_r_spec)
 

    [J,ns]=size(mu_ij);

    %%%%%%%%%%% Kalouptsidi (2012) base spec %%%%%%%
    numer=exp(mu_ij+reshape(r_initial,1,ns));%J*ns
    denom=sum(numer,2);%J*1
    temp1=sum(reshape(S_j_data,J,1).*numer./denom,1);%1*ns

    numer=exp(reshape(r_initial,1,ns));%1*ns
    denom=sum(numer,2);%1*1
    temp2=S_0_data.*numer./denom;

    difference=log(weight)-log(temp1+temp2);
    r_updated=r_initial+difference;


     temp=S_0_data-sum(exp(r_initial(1:end-1)));

  
    if switch_r_spec==1 % mixed algorithm
    if temp>0
        % Kalouptsidi base spec
        r_updated(end)=log(temp);
    elseif temp<=0
            %warning("NaN val")
            %%% Switch to the stable but slow algorithm
            r_tilde_initial=r_initial-r_initial(end);
            r_tilde_updated=r_tilde_initial+difference;
            r_tilde_updated(end)=0;
            r_updated=r_tilde_updated+log(S_0_data)-log(sum(exp(r_tilde_updated)));
            %%%%%
    end% mixed algorithm end

    else % conservative algorithm
        r_updated(end)=0;% Kalouptsidi conservative spec

    end

    %%%%%%%%%%%%%%%%%%%%
    
    output={r_updated};
    other_vars.numer=numer;
    other_vars.denom=denom;
    other_vars.temp1=temp1;
    other_vars.temp2=temp2;
    
    
end
