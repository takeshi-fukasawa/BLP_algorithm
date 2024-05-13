%% Compute delta, given nonlinear parameters (modified version of Kalouptsidi method)
ITER_MAX=spec.ITER_MAX;
r_initial0=log(S_0t_data.*weight);%1 by ns; Initial value
%%%r_initial0=r_initial0+randn(size(r_initial0));

%r_initial0=r_initial0-r_initial0(end);

%r_initial0=log(s_i0t_true)+0.1*rand(1,ns);%1*ns; true answer + noise

%r_initial0=r_normalize_func(r_initial0,S_0t_data);

tic
for kk=1:n_sim
r_initial=r_initial0;

DIST=100;

for k=1:ITER_MAX

    output=...
    r_update_RCL_func(...
    r_initial,weight,mu_ijt_est,rho_est,...
    S_jt_data,S_0t_data,exp(mu_ijt_est),...
    switch_r_spec);
    
    r_updated=output{1};

    %r_updated=r_normalize_func(r_updated,S_0t_data);%%~=> Faster??


    DIST=max(abs(r_updated(:)-r_initial(:)));%scalar
 
    V_updated_1=log(weight)-r_updated;
    delta_updated_1=compute_delta_from_V_func(...
        mu_ijt_est,weight,S_jt_data,V_updated_1);

    %%delta_updated_1=compute_delta_func(numer_1_without_delta,I_igt,weight,S_jt_data,rho_est,0);%J*1

    %DIST=max(abs(delta_jt_true(:)-delta_updated_1(:)));%scalar
    DIST_MAT(k,2)=DIST;
    
    %%%%%
    %k
    %DIST
      if DIST<TOL
          break;% end the for loop
      elseif isnan(DIST)==1
          k=ITER_MAX;
          break;
      else
          r_initial=r_updated;%J by 1
      end

end% for loop

end%%kk
r_updated=r_normalize_func(r_updated,S_0t_data);%%%%%%%%


V_updated_1=log(weight)-r_updated;
delta_updated_1=compute_delta_from_V_func(...
        mu_ijt_est,weight,S_jt_data,V_updated_1);

t_update_s=toc/n_sim;

feval_update_s=k;
ratio_delta_update_s=delta_updated_1./delta_jt_true;
[min(ratio_delta_update_s),max(ratio_delta_update_s)];


results_s(m,1)=feval_update_s;
results_s(m,2)=t_update_s;

results_s(m,3)=(results_s(m,1)<ITER_MAX);

if mistake_spec==0
    results_s(m,3)=results_s(m,3).*(max(abs(ratio_delta_update_s(:)-1))<1e-8);
end


log_dist_new=log10(DIST_MAT);


%%%%%%%%%%%%%%%%%
if 1==1
%% New spectral (Kalouptsidi method)
spec=spec_default;
for kk=1:n_sim
        [output_spectral,other_vars,iter_info]=...
        spectral_func(@r_update_RCL_func,spec,{r_initial0},...
        weight,mu_ijt_est,rho_est,...
        S_jt_data,S_0t_data,exp(mu_ijt_est),switch_r_spec);

    r_sol=output_spectral{1};

% Necessary????
r_sol_temp=r_normalize_func(r_sol,S_0t_data);

V_updated_1=log(weight)-r_updated;

delta_sol=compute_delta_from_V_func(...
        mu_ijt_est,weight,S_jt_data,V_updated_1);


end%%kk
t_update_s_spectral=iter_info.t_cpu/n_sim;

feval_update_s_spectral=iter_info.feval;
ratio_delta_s_spectral=delta_sol./delta_jt_true;

log_dist_spectral=log10(DIST_MAT_spectral);

results_s_spectral(m,1)=feval_update_s_spectral;
results_s_spectral(m,2)=t_update_s_spectral;
results_s_spectral(m,3)=(results_s_spectral(m,1)<ITER_MAX);

if mistake_spec==0
    results_s_spectral(m,3)=results_s_spectral(m,3).*(max(abs(ratio_delta_s_spectral(:)-1))<1e-8);
end
ratio_r=(log(weight)-log(V_true))./r_sol;

end % spectral


