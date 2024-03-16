function [out,other_vars]=delta_middle_Bellman_inner_func(...
    delta_initial,V_initial0,mu_ij,S_j_data,beta_C,L,...
    rho,weight,spectral_V_spec,tune_param_BLP)

spec=spec_default;
if spectral_V_spec==1
    spec.update_spec=4;
else
    spec.update_spec=0;
end

        [output_spectral]=...
        spectral_func(@Bellman_update_func,spec,{V_initial0},...
        delta_initial,mu_ij,beta_C,rho,...
        weight_V,x_V);

    V_sol=output_spectral{1};

out=BLP_update_func(delta_initial,weight,....
    mu_ij-beta_C.*V_sol,rho,S_j_data,tune_param_BLP);


other_vars=[];
end
