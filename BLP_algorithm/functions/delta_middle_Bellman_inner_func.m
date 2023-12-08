function out=delta_middle_Bellman_inner_func(...
    delta_initial,V_initial0,mu_ij,S_j_data,beta_C,L,...
    rho,weight,spectral_V_spec,tune_param_BLP)

global count DIST

ITER_MAX=1500;
TOL=1e-12;
V_initial=V_initial0;


if spectral_V_spec==0

for k=1:ITER_MAX

    out_V=Bellman_update_func(V_initial,delta_initial,...
        mu_ij,beta_C,L,rho,weight);

    resid_V=out_V{1};
    V_updated=V_initial-resid_V;
    DIST=max(abs(resid_V(:)));
    if DIST<TOL
      break;
   else
      V_initial=V_updated;
   end
end % for loop

else %spectral_V_spec==1

        output_spectral=...
        spectral_func_cpu(@Bellman_update_func,1,V_initial0,...
        [],[],...
        delta_initial,mu_ij,beta_C,L,rho,weight);

        count;
        DIST_spectral=DIST;
    V_updated=output_spectral{1};

end

resid_delta=BLP_update_func(delta_initial,weight,....
    mu_ij+(beta_C.^L).*V_updated-beta_C.*V_updated,rho,S_j_data,tune_param_BLP);

 out=resid_delta;

end
