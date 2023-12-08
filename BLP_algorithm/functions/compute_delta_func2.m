%%%%%
function delta=compute_delta_func2(...
    mu_ijt,weight,S_jt_data,rho,V,beta_C,L)

    [J,ns,G,T]=size(mu_ijt);

    temp=exp(reshape(mu_ijt,J,ns,G,T)-reshape(V,1,ns,1,T));

    %%% rho==0 case 
    delta=log(S_jt_data)-log(sum(reshape(weight,1,ns,1,1).*temp,2));%J*1*G*T

end

