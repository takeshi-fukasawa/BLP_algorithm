function output=...
    V_Pr0_update_func(...
    V_initial,Pr0_initial,weight,mu_ijt,rho,...
    S_jt_data,S_0t_data,...
    weight_V,x_V,beta_C,L,tune_param,Newton_spec)

[out_V,other_vars_V]=V_update_func(V_initial,weight,mu_ijt,rho,...
    S_jt_data,S_0t_data,...
    weight_V,x_V,beta_C,L,tune_param,Newton_spec);

resid_V=out_V{1};
s_i0t_ccp=other_vars_V.s_i0t_ccp;

Pr0_updated=Pr0_initial;
Pr0_updated(:,:,:,:,2:end)=Pr0_updated(:,:,:,:,1:end-1).*s_i0t_ccp;

resid_Pr0=Pr0_initial-Pr0_updated;

output={resid_V,resid_Pr0};

end

