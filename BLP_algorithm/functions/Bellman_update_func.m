function [resid_V,other_vars]=Bellman_update_func(...
    V_initial,delta,mu_ijt,beta_C,rho,...
    weight_V,x_V)

    %%% rho>0 allowed; Bellman update

    u_ijt_tilde=delta+mu_ijt;%J*I*G*T*n_dim_V (1 ??)

    [numer_1,denom_1,IV,EV]=compute_IV_EV_func(...
            V_initial,u_ijt_tilde,beta_C,rho,weight_V,x_V);

    v_i0t_tilde=beta_C*EV;%J*I*1*T*n_dim_V

    V_updated=log(exp(v_i0t_tilde)+sum(exp(IV),3));%1*ns*1*T*n_dim_V


resid_V={V_initial-V_updated};
other_vars.IV=IV;
other_vars.numer_1=numer_1;
other_vars.denom_1=denom_1;
other_vars.v_i0t_tilde=v_i0t_tilde;

end
