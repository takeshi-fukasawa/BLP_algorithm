%% Compute delta, given nonlinear parameters (sigma)
DIST=100;

    if tune_param_BLP==0
        results_BLP=results_BLP_contraction;
        results_BLP_spectral=results_BLP_contraction_spectral;
    else
        results_BLP=results_BLP_new;
        results_BLP_spectral=results_BLP_new_spectral;
    end


tic
for kk=1:n_sim
    if rho_est>0
        delta_initial=log(S_jt_data)-log(S_0t_data)-rho_est*log(S_jt_given_g_data);%J*1*G*T; Initial value of delta
    else
        delta_initial=log(S_jt_data)-log(S_0t_data);%J*1*G*T; Initial value of delta
    end
%%delta_initial=delta_initial+0.01;
delta_initial=delta_initial+0*randn(size(delta_initial));

%delta_initial=ones(size(delta_initial));
%delta_initial=delta_jt_true*0.999;%%%%

DIST_MAT=zeros(ITER_MAX,3);
for iter_BLP=1:ITER_MAX

  output=BLP_update_func(...
    delta_initial,weight,mu_ijt_est,rho_est,...
    S_jt_data,tune_param_BLP);

   resid=output{1};

   delta_updated=delta_initial-resid;
  
  DIST=max(abs(delta_updated(:)-delta_initial(:)));%scalar
  
  %%%%%
  %DIST=max(abs(delta_jt_true(:)-delta_updated(:)));%scalar
  %%%%%
  DIST_MAT(iter_BLP,1)=DIST;

  if DIST<TOL
      break;% end the for loop
  elseif isnan(DIST)==1
      iter_BLP=ITER_MAX;
      break;
  else
      delta_initial=delta_updated;%J by 1
  end

end% for loop
end
t_BLP=toc/n_sim;
ratio_delta_BLP=delta_updated./delta_jt_true;

n_iter_BLP=iter_BLP;

results_BLP(m,1)=n_iter_BLP;
results_BLP(m,2)=t_BLP;
results_BLP(m,3)=(results_BLP(m,1)<ITER_MAX);

DIST_MAT_BLP=DIST_MAT;

if mistake_spec==0
    results_BLP(m,3)=results_BLP(m,3).*(max(abs(ratio_delta_BLP(:)-1))<1e-8);
end

    %% BLP spectral
    delta_initial0=log(S_jt_data)-log(S_0t_data)-rho_est.*log(S_jt_given_g_data);%J by 1; Initial value of delta

    tic
    for kk=1:n_sim
    if GPU_spec==1
        output_BLP_spectral=...
        spectral_func(@BLP_update_func,1,delta_initial0,...
        [],[],...
        weight,mu_ijt_est,rho_est,S_jt_data,tune_param_BLP);
    else
        output_BLP_spectral=...
    spectral_func_cpu(@BLP_update_func,1,delta_initial0,...
    [],[],...
    weight,mu_ijt_est,rho_est,S_jt_data,tune_param_BLP);

    delta_sol=output_BLP_spectral{1};

    end
    end%% kk

    t_BLP_spectral=toc/n_sim;

    n_iter_BLP_spectral=count;
    DIST_MAT_BLP_spectral=DIST_table;
    

    log_dist_BLP_spectral=log10(DIST_MAT_BLP_spectral);

    ratio_delta_BLP_spectral=delta_sol./delta_jt_true;
    delta_check_BLP_spectral=[min(ratio_delta_BLP_spectral),max(ratio_delta_BLP_spectral)];

    results_BLP_spectral(m,1)=n_iter_BLP_spectral;
    results_BLP_spectral(m,2)=t_BLP_spectral;
    results_BLP_spectral(m,3)=(results_BLP_spectral(m,1)<ITER_MAX);
    
    if mistake_spec==0
        results_BLP_spectral(m,3)=(results_BLP_spectral(m,1)<ITER_MAX & max(abs(ratio_delta_BLP_spectral(:)-1))<1e-8);
    end

    if tune_param_BLP==0
        results_BLP_contraction=results_BLP;
        results_BLP_contraction_spectral=results_BLP_spectral;
    else
        results_BLP_new=results_BLP;
        results_BLP_new_spectral=results_BLP_spectral;
    end
        
%%%%%%%%%%%%%%
