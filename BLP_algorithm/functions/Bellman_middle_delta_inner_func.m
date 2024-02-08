function [out,other_vars]=Bellman_middle_delta_inner_func(...
    V_initial,delta_initial0,mu_ij,S_j_data,weight_V,x_V,...
    beta_C,L,...
    rho,weight,spectral_delta_spec,tune_param_BLP,delta_fixed_spec)

global delta_updated count DIST

ITER_MAX=1000;
TOL=1e-12;
delta_initial=delta_initial0;

if delta_fixed_spec==0
if spectral_delta_spec==0

for iter=1:ITER_MAX

    out_delta=BLP_update_func(delta_initial,weight,....
        mu_ij+(beta_C.^L).*V_initial-beta_C.*V_initial,rho,S_j_data,tune_param_BLP);

    resid_delta=out_delta{1};
    delta_updated=delta_initial-resid_delta;
    if max(abs(resid_delta(:)))<TOL
      break;
   else
      delta_initial=delta_updated;
   end
end % for loop


else %spectral_delta_spec==1
        output_spectral=...
        spectral_func(@BLP_update_func,1,[],[],delta_initial0,...
        weight,mu_ij+(beta_C.^L).*V_initial-beta_C.*V_initial,...
        rho,S_j_data,tune_param_BLP);

        count;
        DIST_spectral=DIST;
    delta_updated=output_spectral{1};

end

elseif delta_fixed_spec==1
    delta_updated=delta_initial0;
end

 resid_V=Bellman_update_func(V_initial,delta_updated,mu_ij,beta_C,L,rho,weight,...
 weight_V,x_V);

 out=resid_V;

other_vars=[];

end
