%%%%%
function [delta,s_ijt_up_to_scale,s_igt]=compute_delta_from_V_IV_func(...
    mu_ijt,weight,S_jt_data,rho,V,IV)

%%% rho>0 case allowed
%%% IV,V => delta
%%% If V empty, compute V from IV (static (beta_C==0) case)
%%% Currently, Pr0_spec==0 case only

    [J,ns,G,T]=size(mu_ijt);

    if isempty(V)==1 % static (beta_C==0) case
       V=log(1+sum(exp(IV),3));%1*ns*G*T*n_dim_V
    end

    temp_numer1=exp(mu_ijt./(1-rho));
    temp_denom1=exp(IV./(1-rho));
    temp_numer2=exp(IV);
    temp_denom2=exp(V);

    s_igt=temp_numer2./temp_denom2;
    s_ijt_up_to_scale=temp_numer1./temp_denom1.*s_igt;%J*ns*G*T
    n_state=size(weight,5);
    s_jt_up_to_scale=sum(reshape(weight,1,ns,1,1,n_state).*...
            s_ijt_up_to_scale,[2,5]);%J*1*G*T

    delta=(1-rho).*(log(S_jt_data)-log(s_jt_up_to_scale));%J*1*G*T

end

