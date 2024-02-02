
spectral_delta_spec=1;
tune_param_BLP=1;
%delta_fixed_spec=1;

%% Normal iteration
%%%V_initial0=V_true;
if delta_fixed_spec==0
    delta_initial00=delta_initial0;
else
    delta_initial00=delta_jt_true;
end
V_initial=V_initial0;
%V_initial=V_true;%%%%%%%%%%%%%%

delta_updated=0;

DIST_MAT=zeros(ITER_MAX,1);
tic
for k=1:ITER_MAX

    resid_cell=Bellman_middle_delta_inner_func(...
        V_initial,delta_initial00,mu_ijt_est,S_jt_data,beta_C,L,...
    rho_est,weight,spectral_delta_spec,tune_param_BLP,delta_fixed_spec);

    V_updated=V_initial-resid_cell{1};
    DIST=max(abs(V_updated(:)-V_initial(:)));%scalar
    DIST_MAT(k,:)=DIST;


      if DIST<TOL
          break;% end the for loop
      elseif isnan(DIST)==1
          k=ITER_MAX;
          break;
      else
          V_initial=V_updated;%J by 1
      end

end% for loop
t_V_middle=toc;
ratio_V_middle=delta_updated./delta_jt_true;
n_iter_update_V_middle=k;

%% Spectral
if delta_fixed_spec==0
    delta_initial00=delta_initial0;
else
    delta_initial00=delta_jt_true;
end

tic
output_spectral=...
        spectral_func_cpu(@Bellman_middle_delta_inner_func,...
        1,V_initial0,...
        [],[],...
        delta_initial00,mu_ijt_est,S_jt_data,beta_C,L,...
    rho_est,weight,spectral_delta_spec,tune_param_BLP,delta_fixed_spec);
    
    V_sol=output_spectral{1};
    delta_sol=delta_updated;

t_V_middle_spectral=toc;

ratio_V_middle_spectral=delta_sol./delta_jt_true;

DIST_MAT_V_middle_spectral=DIST_table;
n_iter_update_V_middle_spectral=count;
