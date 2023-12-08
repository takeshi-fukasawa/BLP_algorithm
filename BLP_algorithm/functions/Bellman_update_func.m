function resid_V=Bellman_update_func(...
    V_initial,delta,mu_ijt,beta_C,L,rho,weight)


u_ijt_tilde=delta+mu_ijt+(beta_C^L).*V_initial;%J*I*G*T

weight_V=[];
EV=compute_EV_func(V_initial,weight_V);
u_i0t_tilde=beta_C*EV;%J*I*1*T


[s_jt,ChoiceProb,s_ijt_given_g_ccp,s_igt_ccp,...
numer_1,denom_1,numer_2,denom_2]=...
share_func(u_ijt_tilde,u_i0t_tilde,rho,weight);

V_updated=log(denom_2);
resid_V={V_initial-V_updated};

end
