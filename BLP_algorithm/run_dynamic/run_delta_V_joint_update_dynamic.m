

%% BLP_Bellman_joint_update_func 
dump_param=[];
for method=1:2
if method==1 % fixed point iteration
    vec=0;
elseif method==2 % spectral
   vec=t_dim_id*ones(1,2);
   %vec=[];
   %dump_param=[0.1,0.5];
end

[output_spectral,other_vars,DIST_MAT,iter_info]=...
        spectral_func(@BLP_Bellman_joint_update_func,2,vec,dump_param,...
        delta_initial0,V_initial0,...
        weight,mu_ijt_est,rho_est,...
    S_jt_data,weight_V,x_V,beta_C,L,tune_param_BLP);
    

    delta_sol=output_spectral{1};
    V_sol=output_spectral{2};


results=results_output_func(iter_info,other_vars.s_jt_predict,S_jt_data);
ratio_delta=delta_sol./delta_jt_true;

if method==1

if tune_param_BLP==0
results_V_BLP_0(m,:)=results;
ratio_delta_V_BLP_0=ratio_delta;
elseif tune_param_BLP==1
results_V_BLP_1(m,:)=results;
ratio_delta_V_BLP_1=ratio_delta;
end% tune_param_BLP==0 or 1

elseif method==2

if tune_param_BLP==0
results_V_BLP_0_spectral(m,:)=results;
ratio_delta_V_BLP_0_spectral=ratio_delta;
elseif tune_param_BLP==1
results_V_BLP_1_spectral(m,:)=results;
ratio_delta_V_BLP_1_spectral=ratio_delta;
end% tune_param_BLP==0 or 1
end

end % method
