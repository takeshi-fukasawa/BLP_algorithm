%%%%%
function delta=compute_delta_from_V_func(...
    mu_ijt,weight,S_jt_data,rho,V)

    %%% mu_ijt: include continuation values
    [J,ns,G,T]=size(mu_ijt);

    temp=exp(reshape(mu_ijt,J,ns,G,T)-reshape(V,1,ns,1,T));

    %%% rho==0 case 
    delta=log(S_jt_data)-log(sum(reshape(weight,1,ns,1,1).*temp,2));%J*1*G*T

end

