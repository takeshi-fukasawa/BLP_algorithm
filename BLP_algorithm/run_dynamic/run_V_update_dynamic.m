
%% BLP_Bellman_joint_update_func 
dump_param=[];
for method=1:1
if method==1 % fixed point iteration
    vec=0;
elseif method==2 % spectral
   vec=t_dim_id*ones(1,2);
   %vec=[];
end

[delta_jt_temp,Pr0_temp,s_ijt_ccp_up_to_scale_temp]=compute_delta_from_V_func(...
       mu_ijt_true,weight,S_jt_data,V_true(:,:,:,:,1));

%%%%%%%%%%%%%%%%%%%%%%
s_ijt_ccp_up_to_scale_temp2=exp(reshape(mu_ijt_true,J,ns,G,T,1)-...
     reshape(V_true(:,:,:,:,1),1,ns,1,T));%J*ns*G*T
s_ijt_ccp_temp=s_ijt_ccp_up_to_scale_temp.*exp(delta_jt_temp);
ratio=s_ijt_ccp_temp./s_ijt_ccp_true-1

%%%%%%%%%%%%%
s_ijt_ccp_temp0=exp(delta_jt_true+...
    reshape(mu_ijt_true,J,ns,G,T,1)-...
        reshape(V_true(:,:,:,:,1),1,ns,1,T,1));%J*ns*G*T*n_dim_V
s_ijt_ccp_temp0./s_ijt_ccp_true-1
%%%%%%%%%%%%%%%%

[output_spectral,other_vars,DIST_table_spectral,iter_info]=...
        spectral_func(@V_update_func,1,vec,dump_param,...
        V_true,...
        weight,mu_ijt_est,...
        S_jt_data,S_0t_data,...
        weight_V,x_V,beta_C,tune_param,Newton_spec);

    V_sol=output_spectral{1};
    delta_sol=other_vars.delta_jt;

    IV=other_vars.IV;
    v_i0t_tilde=other_vars.v_i0t_tilde;
    v_ijt_tilde=other_vars.v_ijt_tilde;
    
    s_jt_predict=...
    share_func(v_ijt_tilde(:,:,:,:,1),...
        v_i0t_tilde(:,:,:,:,1),rho_est,weight);


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
