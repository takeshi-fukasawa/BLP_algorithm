%% IV_update_func

if tune_param==0
    results_IV=results_IV_contraction;
    results_IV_spectral=results_IV_contraction_spectral;
else
    results_IV=results_IV_new;
    results_IV_spectral=results_IV_new_spectral;
end

%%%IV_initial0=IV_true;
IV_initial0=repmat(log(S_gt_data)-log(S_0t_data),[1,ns,1,1]);%1*ns*G*T

IV_initial=IV_initial0;
%%IV_initial=IV_true;

DIST_MAT=zeros(ITER_MAX,1);
tic
for iter_IV=1:ITER_MAX
  
    output=...
    IV_update_func(...
    IV_initial,weight,mu_ijt_est,rho_est,...
    S_jt_data,S_0t_data,tune_param);
    resid=output{1};

    IV_updated=IV_initial-resid;

    DIST=max(abs(IV_updated(:)-IV_initial(:)));%scalar
    DIST_MAT(iter_IV,:)=DIST;


    %%%%%
      if DIST<TOL
          break;% end the for loop
      elseif isnan(DIST)==1
          iter_IV=ITER_MAX;
          break;
      else
          IV_initial=IV_updated;%J by 1
      end


end% for loop

delta_updated=compute_delta_from_V_IV_func(mu_ijt_est,weight,...
        S_jt_data,rho_est,...
        [],IV_updated);%J*1
t_update_IV=toc;


ratio_delta_IV=delta_updated./delta_jt_true;

DIST_MAT_IV=DIST_MAT;

n_iter_update_IV=iter_IV;

      results_IV(m,1)=n_iter_update_IV;
      results_IV(m,2)=t_update_IV;
      results_IV(m,3)=(results_IV(m,1)<ITER_MAX);

[s_jt_predict,~]=...
  share_func(delta_updated+mu_ijt_est,zeros(1,ns,1,T),rho_est,weight);%J*1*G*T
DIST_s_jt=max(abs(log(s_jt_predict(:))-log(S_jt_data(:))));
results_IV(m,4)=log10(DIST_s_jt);
results_IV(m,5)=(results_IV(m,4)<log10(TOL_DIST_s_jt));

     
%% IV_update_func spectral
tic
for kk=1:n_sim
        output_spectral=...
        spectral_func(@IV_update_func,1,[],[],IV_initial0,...
        weight,mu_ijt_est,rho_est,...
        S_jt_data,S_0t_data,tune_param);

    IV_sol=output_spectral{1};

end%%kk
t_update_IV_spectral=toc/n_sim;


delta_updated=compute_delta_from_V_IV_func(mu_ijt_est,weight,...
        S_jt_data,rho_est,...
        [],IV_sol);%J*1


n_iter_update_IV_spectral=count;
DIST_MAT_IV_spectral=DIST_table;
ratio_delta_IV_spectral=delta_sol./delta_jt_true;

results_IV_spectral(m,1)=n_iter_update_IV_spectral;
results_IV_spectral(m,2)=t_update_IV_spectral;
results_IV_spectral(m,3)=(results_IV_spectral(m,1)<ITER_MAX);

[s_jt_predict,~]=...
  share_func(delta_sol+mu_ijt_est,zeros(1,ns,1,T),rho_est,weight);%J*1*G*T
DIST_s_jt=max(abs(log(s_jt_predict(:))-log(S_jt_data(:))));
results_IV_spectral(m,4)=log10(DIST_s_jt);
results_IV_spectral(m,5)=(results_IV_spectral(m,4)<log10(TOL_DIST_s_jt));

if mistake_spec==0
    results_IV_spectral(m,3)=(results_IV_spectral(m,1)<ITER_MAX & ...
        max(abs(ratio_delta_IV_spectral(:)-1))<1e-8);
end

if tune_param==0
    results_IV_contraction=results_IV;
    results_IV_contraction_spectral=results_IV_spectral;
else
    results_IV_new=results_IV;
    results_IV_new_spectral=results_IV_spectral;
end
       