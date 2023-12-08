    %%%%%
function delta=compute_delta_func(numer_1_without_delta,I_igt,weight,S_jt_data,rho,outside_util)

    temp=(numer_1_without_delta.*exp(-rho./(1-rho).*I_igt))./(exp(outside_util)+sum(exp(I_igt),3));%J*ns*G*T

    delta=(1-rho).*(log(S_jt_data)-log(sum(weight.*temp,2)));%J*1*G*T

end
