%% IV_update_func
%%% Jointly update IV and V

IV_initial0=repmat(...
    log(exp(V_initial0)-exp(beta_C*V_initial0)),...
    1,1,G,1);%1*ns*G*T
IV_initial=IV_initial0;
V_initial=V_initial0*1.1;
%IV_initial=IV_true;
%V_initial=V_true;

tune_param=1/(1-beta_C);% fast, even when rho>0.

tune_param=1;

DIST_MAT=zeros(ITER_MAX,2);
tic
%ITER_MAX=1000;
for k=1:ITER_MAX
    output=...
    IV_update_func(...
    V_initial,IV_initial,weight,mu_ijt_est,rho_est,...
    S_jt_data,S_0t_data,beta_C,L,tune_param);
    
    resid_V=output{1};
    resid_IV=output{2};
    V_updated=V_initial-resid_V;
    IV_updated=IV_initial-resid_IV;
    
    DIST_IV=max(abs(IV_updated(:)-IV_initial(:)));%scalar
    DIST_V=max(abs(V_updated(:)-V_initial(:)));%scalar
    DIST=max(DIST_IV,DIST_V);

    DIST_MAT(k,1)=DIST_V;
    DIST_MAT(k,2)=DIST_IV;
    


    %%%%%
    %k
    %DIST
      if DIST<TOL
          break;% end the for loop
      elseif isnan(DIST)==1
          k=ITER_MAX;
          break;
      else
          V_initial=V_updated;
          IV_initial=IV_updated;
          
      end

end% for loop
t_update_IV=toc;

delta_updated=compute_delta_from_V_IV_func(mu_ijt_est,weight,...
        S_jt_data,rho_est,...
        V_updated,IV_updated,beta_C,L);%J*1

ratio_delta_IV=delta_updated./delta_jt_true;

n_iter_update_IV=k;
%%%%%%%%%%%%%%%%%%%%%

%% IV_update_func spectral
tic
for kk=1:n_sim
        output_spectral=...
        spectral_func_cpu(@IV_update_func,2,...
        V_initial0,IV_initial0,...
        [],[],[],[],...
        weight,mu_ijt_est,rho_est,...
        S_jt_data,S_0t_data,...
        beta_C,L,tune_param);

    V_sol=output_spectral{1};
    IV_sol=output_spectral{2};
    
end%%kk

n_iter_update_IV=k;
t_update_IV_spectral=toc/n_sim;

delta_sol=compute_delta_from_V_IV_func(...
        mu_ijt_est,weight,S_jt_data,rho_est,...
        V_sol,IV_sol,beta_C,L);

n_iter_update_IV_spectral=count;
DIST_MAT_IV_spectral=DIST_table;
ratio_delta_IV_spectral=delta_sol./delta_jt_true;
%[min(ratio_delta_IV_spectral),max(ratio_delta_IV_spectral)]
