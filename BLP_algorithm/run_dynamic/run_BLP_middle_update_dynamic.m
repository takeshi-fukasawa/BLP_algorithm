spectral_V_spec=1;
tune_param_BLP=0;

%% Normal iteration
%%%V_initial0=V_true;
V_initial=V_initial0;
delta_initial=delta_initial0;

DIST_MAT=zeros(ITER_MAX,1);
tic
for k=1:ITER_MAX

    out_cell=delta_middle_Bellman_inner_func(delta_initial,...
        V_initial0,mu_ijt_est,S_jt_data,beta_C,L,rho_est,...
        weight,spectral_V_spec,tune_param_BLP);

    delta_updated=out_cell{1};
    DIST=max(abs(delta_updated(:)-delta_initial(:)));%scalar
    DIST_MAT(k,:)=DIST;


      if DIST<TOL
          break;% end the for loop
      elseif isnan(DIST)==1
          k=ITER_MAX;
          break;
      else
          delta_initial=delta_updated;%J by 1
      end

end% for loop
t_delta_middle=toc;
ratio_delta_middle=delta_updated./delta_jt_true;
n_iter_update_delta_middle=k;

spec.update_spec=t_dim_id;

tic
    [output_spectral,other_vars,iter_info]=...
        spectral_func(@delta_middle_Bellman_inner_func,...
        spec,{delta_initial0},V_initial0,...
        mu_ijt_est,S_jt_data,beta_C,L,rho_est,...
        weight,spectral_V_spec,tune_param_BLP);
    
    delta_sol=output_spectral{1};

t_delta_middle_spectral=iter_info.t_cpu;
IV_state=other_vars.IV;

ratio_delta_middle_spectral=delta_sol./delta_jt_true;
DIST_MAT_delta_middle_spectral=DIST_table;
n_iter_update_delta_middle_spectral=iter_info.n_iter;
