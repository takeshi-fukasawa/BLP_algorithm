%%%%%
function [delta,s_ijt_up_to_scale,s_igt]=compute_delta_from_V_IV_func(...
    mu_ijt,weight,S_jt_data,rho,V,IV,beta_C,L)

    [J,ns,G,T]=size(mu_ijt);

    temp_numer1=exp(mu_ijt./(1-rho));
    temp_denom1=exp(IV./(1-rho));
    temp_numer2=exp(IV);
    temp_denom2=exp(V);

    s_igt=temp_numer2./temp_denom2;
    s_ijt_up_to_scale=temp_numer1./temp_denom1.*s_igt;%J*ns*G*T
    s_jt_up_to_scale=sum(reshape(weight,1,ns).*...
            s_ijt_up_to_scale,2);%J*1*G*T

    delta=(1-rho).*(log(S_jt_data)-log(s_jt_up_to_scale));%J*1*G*T


end

