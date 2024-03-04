function [out,other_vars]=delta_middle_Bellman_inner_func(...
    delta_initial,V_initial0,mu_ij,S_j_data,beta_C,L,...
    rho,weight,spectral_V_spec,tune_param_BLP)

global count DIST

spec=[];
if spectral_V_spec==1
    spec.update_spec=4;
else
    spec.update_spec=0;
end

        output_spectral=...
        spectral_func(@Bellman_update_func,spec,{V_initial0},...
        delta_initial,mu_ij,beta_C,rho,...
        weight_V,x_V);

    V_updated=output_spectral{1};

resid_delta=BLP_update_func(delta_initial,weight,....
    mu_ij-beta_C.*V_updated,rho,S_j_data,tune_param_BLP);

 out=resid_delta;

other_vars=[];
end
