%% Joint s update spectral
DIST=100;
DIST_MAT_joint_s=zeros(ITER_MAX,3);

results_joint_s_spectral=zeros(n_market,3);
results_joint_s=zeros(n_market,3);

log_s_i0t_ccp_initial0=log(S_0t_data)+0*weight;
log_s_igt_ccp_initial0=log(S_gt_data)+0*weight;


tic
for kk=1:n_sim
    log_s_i0t_ccp_initial=log_s_i0t_ccp_initial0;
    log_s_igt_ccp_initial=log_s_igt_ccp_initial0;

for k=1:ITER_MAX

    output=...
    joint_s_update_RCNL_func(...
    log_s_i0t_ccp_initial,log_s_igt_ccp_initial,weight,mu_ijt_est,rho_est,...
    S_jt_data,S_gt_data,S_0t_data,numer_1_without_delta);
    resid_log_s_i0t_ccp=output{1};
    resid_log_s_igt_ccp=output{2};

    log_s_i0t_ccp_updated=log_s_i0t_ccp_initial-resid_log_s_i0t_ccp;
    log_s_igt_ccp_updated=log_s_igt_ccp_initial-resid_log_s_igt_ccp;

    delta_updated_joint_s=compute_delta_func(numer_1_without_delta,I_igt,weight,S_jt_data,rho_est,0);%J*1

    DIST1=max(abs(log_s_i0t_ccp_updated(:)-log_s_i0t_ccp_initial(:)));
    DIST2=max(abs(log_s_igt_ccp_updated(:)-log_s_igt_ccp_initial(:)));

    %DIST=max(abs(delta_jt_true(:)-delta_updated_joint_s(:)));%scalar
    DIST_MAT_joint_s(k,1)=DIST1;
    DIST_MAT_joint_s(k,2)=DIST2;
    DIST_MAT_joint_s(k,3)=DIST;

    DIST=max(DIST1,DIST2);
    
    %%%%%
    %k
    %DIST
      if DIST<TOL
          break;% end the for loop
      elseif isnan(DIST)==1
          k=ITER_MAX;
          break;
      else
          log_s_i0t_ccp_initial=log_s_i0t_ccp_updated;%J by 1
          log_s_igt_ccp_initial=log_s_igt_ccp_updated;%J by 1
      end

end% for loop
end %%kk

t_joint_s=toc/n_sim;


n_iter_update_joint_s=k;
ratio_delta_joint_s=delta_updated_joint_s./delta_jt_true;
[min(ratio_delta_joint_s),max(ratio_delta_joint_s)];

log_dist_joint_s_new=log10(DIST_MAT_joint_s);

results_joint_s(m,1)=n_iter_update_joint_s;
results_joint_s(m,2)=t_joint_s;
results_joint_s(m,3)=(results_joint_s(m,1)<ITER_MAX);

if mistake_spec==0
    results_joint_s(m,3)=(results_joint_s(m,1)<ITER_MAX & max(abs(ratio_delta_joint_s(:)-1))<1e-8);
end
    %plot(log_dist_new)
%%%%%

if 1==1

tic;
for kk=1:n_sim
if GPU_spec==1
    output_spectral=...
    spectral_func(@joint_s_update_RCNL_func,2,...
    log_s_i0t_ccp_initial0,log_s_igt_ccp_initial0,...
    [],[],[],[],...
    weight,mu_ijt_est,rho_est,...
    S_jt_data,S_gt_data,S_0t_data,numer_1_without_delta);

else
     output_spectral=...
    spectral_func_cpu(@joint_s_update_RCNL_func,2,...
    log_s_i0t_ccp_initial0,log_s_igt_ccp_initial0,...
    [],[],[],[],...
    weight,mu_ij_est,rho_est,...
    S_jt_data,S_gt_data,S_0t_data,numer_1_without_delta);
   
end
delta_sol_joint_s_spectral=compute_delta_func(numer_1_without_delta,I_igt,weight,S_jt_data,rho_est,0);%J*1

end%%kk
    t_joint_s_spectral=toc/n_sim;

    log_s_i0t_ccp_sol=output_spectral{1};
    log_s_igt_ccp_sol=output_spectral{2};
    


n_iter_joint_s_spectral=count;
DIST_MAT_joint_s_spectral=DIST_table;
ratio_delta_joint_s_spectral=delta_sol_joint_s_spectral./delta_jt_true;
[min(ratio_delta_joint_s_spectral),max(ratio_delta_joint_s_spectral)];

log_dist_joint_s_spectral=log10(DIST_MAT_joint_s_spectral);

results_joint_s_spectral(m,1)=n_iter_joint_s_spectral;
results_joint_s_spectral(m,2)=t_joint_s_spectral;

results_joint_s_spectral(m,3)=(results_joint_s_spectral(m,1)<ITER_MAX);

    if mistake_spec==0
        results_joint_s_spectral(m,3)=(results_joint_s_spectral(m,1)<ITER_MAX & max(abs(ratio_delta_joint_s_spectral(:)-1))<1e-8);
    end
end

%%%%%%%%%%%%%
