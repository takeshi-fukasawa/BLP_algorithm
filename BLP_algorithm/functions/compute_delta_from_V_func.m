%%%%%
function delta=compute_delta_from_V_func(...
    mu_ijt,weight,S_jt_data,rho,V)

    %%% mu_ijt: include continuation values
    [J,ns,G,T]=size(mu_ijt);
    n_dim_V=size(V,5);
    temp=exp(reshape(mu_ijt,J,ns,G,T,1)-reshape(V,1,ns,1,T,n_dim_V));

    %%% rho==0 case 
    n_state=size(weight,5); %% n_dim_V==1 => n_state=1; n_dim_V==n_state>1 under multiple states; IVS => n_dim_V==n_state==1 (only one obs pt)
    delta=log(S_jt_data)-log(sum(reshape(weight,1,ns,1,1,n_state).*temp,[2,5]));%J*1*G*T

end

