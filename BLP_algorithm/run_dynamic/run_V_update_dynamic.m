

%% BLP_Bellman_joint_update_func 
dump_param=[];
for method=1:2
if method==1 % fixed point iteration
    vec=0;
elseif method==2 % spectral
   vec=t_dim_id*ones(1,2);
end

[output_spectral,other_vars,DIST_table_spectral,iter_info]=...
        spectral_func(@V_update_func,1,vec,dump_param,...
        V_initial0,...
        weight,mu_ijt_est,...
        S_jt_data,S_0t_data,...
        weight_V,x_V,beta_C,tune_param,Newton_spec);

    V_sol=output_spectral{1};
    delta_sol=log(other_vars.exp_delta_jt);

%%%%%%%%
s_jt_predict=S_jt_data;%%%%%
%%%%%%%%%

results=results_output_func(iter_info,s_jt_predict,S_jt_data);
ratio_delta=delta_sol./delta_jt_true;


if method==1

if tune_param==0
results_V_0(m,:)=results;
ratio_delta_V_0=ratio_delta;
elseif tune_param==1
results_V_1(m,:)=results;
ratio_delta_V_1=ratio_delta;
elseif tune_param>1
results_V_2(m,:)=results;
ratio_delta_V_2=ratio_delta;
end

elseif method==2

if tune_param==0
results_V_0_spectral(m,:)=results;
ratio_delta_V_0_spectral=ratio_delta;
elseif tune_param==1
results_V_1_spectral(m,:)=results;
ratio_delta_V_1_spectral=ratio_delta;
elseif tune_param>1
results_V_2_spectral(m,:)=results;
ratio_delta_V_2_spectral=ratio_delta;

end% tune_param==0 or 1 or others
end

end % method
