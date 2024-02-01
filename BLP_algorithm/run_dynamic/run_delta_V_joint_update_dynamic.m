
if tune_param_BLP==0
    results_V_BLP=results_V_BLP_0;
    results_V_BLP_spectral=results_V_BLP_0_spectral;
elseif tune_param_BLP==1
    results_V_BLP=results_V_BLP_1;
    results_V_BLP_spectral=results_V_BLP_1_spectral;
end

%% BLP_Bellman_joint_update_func 

delta_initial=delta_initial0;
V_initial=V_initial0;

%delta_initial=delta_jt_true;
%V_initial=V_true;

DIST_MAT_V_BLP=zeros(ITER_MAX,2);
tic
for k=1:ITER_MAX

 out=BLP_Bellman_joint_update_func(...
    delta_initial,V_initial,weight,mu_ijt_est,rho_est,...
    S_jt_data,beta_C,L,tune_param_BLP);

    resid_delta=out{1};
    resid_V=out{2};

    delta_updated=delta_initial-resid_delta;
    V_updated=V_initial-resid_V;

    DIST_V=max(abs(V_updated(:)-V_initial(:)));%scalar
    DIST_delta=max(abs(delta_updated(:)-delta_initial(:)));%scalar

   DIST=max(DIST_delta,DIST_V);

   DIST_MAT_V_BLP(k,1)=DIST_delta;
   DIST_MAT_V_BLP(k,2)=DIST_V;
   
    %%%%%
    %k
    %DIST
      if DIST<TOL
          break;% end the for loop
      elseif isnan(DIST)==1
          k=ITER_MAX;
          break;
      else
          delta_initial=delta_updated;
          V_initial=V_updated;%J by 1
      end

end% for loop
t_V_BLP=toc;

ratio_delta_V_BLP=delta_updated./delta_jt_true;

n_iter_update_V_BLP=k;

EV=compute_EV_func(V_sol,[],weight_V,x_V);
[s_jt_predict,~]=...
  share_func(delta_sol+mu_ijt_est,beta_C*EV,rho_est,weight);%J*1*G*T
DIST_s_jt=max(abs(log(s_jt_predict(:))-log(S_jt_data(:))));

results_V_BLP(m,1)=n_iter_update_V_BLP;
results_V_BLP(m,2)=t_V_BLP;
results_V_BLP(m,3)=(results_V_BLP(m,1)<ITER_MAX);
results_V_BLP(m,4)=log10(DIST_s_jt);
results_V_BLP(m,5)=(results_V_BLP(m,4)<log10(TOL));

%% BLP_Bellman_joint_update_func Spectral
tic
output_spectral=...
        spectral_func_cpu(@BLP_Bellman_joint_update_func,2,...
        delta_initial0,V_initial0,...
        [],[],[],[],...
        weight,mu_ijt_est,rho_est,...
    S_jt_data,beta_C,L,tune_param_BLP);
    

    delta_sol=output_spectral{1};
    V_sol=output_spectral{2};
t_V_BLP_spectral=toc;

ratio_delta_V_BLP_spectral=delta_sol./delta_jt_true;
DIST_MAT_V_BLP_spectral=DIST_table;
n_iter_update_V_BLP_spectral=count;

EV=compute_EV_func(V_sol,[],weight_V,x_V);
[s_jt_predict,~]=...
  share_func(delta_sol+mu_ijt_est,beta_C*EV,rho_est,weight);%J*1*G*T
DIST_s_jt_spectral=max(abs(log(s_jt_predict(:))-log(S_jt_data(:))));


results_V_BLP_spectral(m,1)=n_iter_update_V_BLP_spectral;
results_V_BLP_spectral(m,2)=t_V_BLP_spectral;
results_V_BLP_spectral(m,3)=(results_V_BLP_spectral(m,1)<ITER_MAX);
results_V_BLP_spectral(m,4)=log10(DIST_s_jt_spectral);
results_V_BLP_spectral(m,5)=(results_V_BLP_spectral(m,4)<log10(TOL));

if tune_param_BLP==0
    results_V_BLP_0=results_V_BLP;
    results_V_BLP_0_spectral=results_V_BLP_spectral;
elseif tune_param_BLP==1
    results_V_BLP_1=results_V_BLP;
    results_V_BLP_1_spectral=results_V_BLP_spectral;
end
