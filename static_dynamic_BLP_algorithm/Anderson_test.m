spec=spec_default;
spec.SQUAREM_spec=1;

[output_SQUAREM,other_vars_SQUAREM,...
    iter_info_Bellman_SQUAREM]=...
        spectral_func(@Bellman_update_func,spec,{V_initial0},...
        delta_jt_true,mu_ijt,...
    beta_C,rho_est,weight_V,x_V);
V_updated_SQUAREM=output_SQUAREM{1};

%% Anderson acceleration
spec.SQUAREM_spec=2;
[output_Anderson,other_vars_Anderson,...
    iter_info_Bellman_Anderson]=...
        spectral_func(@Bellman_update_func,spec,{V_initial0},...
        delta_jt_true,mu_ijt,...
    beta_C,rho_est,weight_V,x_V);
V_updated_Anderson=output_Anderson{1};
