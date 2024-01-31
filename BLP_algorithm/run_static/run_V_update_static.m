%% V_update_func

if tune_param==0
    results_V=results_V_contraction;
    results_V_spectral=results_V_contraction_spectral;
else
    results_V=results_V_new;
    results_V_spectral=results_V_new_spectral;
end

%%%V_initial0=V_true;
V_initial0=-log(S_0t_data).*ones(size(weight));
V_initial=V_initial0;

%%%%%%%%%%%
%V_initial=V_initial0*10+randn(size(V_initial0));
%%%%%%%%%%

DIST_MAT=zeros(ITER_MAX,1);
tic
for iter_V=1:ITER_MAX

    output=...
    V_update_func(...
    V_initial,weight,mu_ijt_est,rho_est,...
    S_jt_data,S_0t_data,numer_1_without_delta,beta_C,L,tune_param,Newton_spec);
    resid=output{1};

    V_updated=V_initial-resid;

    %%% The operation is equivalent to multiply (S_0_predict/S_0_data)
    if 1==0
        s_i0t_ccp_temp=exp((beta_C-1).*V_updated);
        s_i0t_ccp=S_0t_data./sum(weight.*s_i0t_ccp_temp,2).*s_i0_ccp_temp;
        s_i0t=s_i0t_ccp.*weight;
        V_updated=log(s_i0t_ccp)./(beta_C-1);
    end

    DIST=max(abs(V_updated(:)-V_initial(:)));%scalar
    DIST_MAT(iter_V,:)=DIST;


    %%%%%
      if DIST<TOL
          break;% end the for loop
      elseif isnan(DIST)==1
          iter_V=ITER_MAX;
          break;
      else
          V_initial=V_updated;%J by 1
      end


end% for loop

delta_updated=compute_delta_from_V_func(mu_ijt_est,weight,...
        S_jt_data,rho_est,...
        V_updated);%J*1
t_update_V=toc;


ratio_delta_V=delta_updated./delta_jt_true;

DIST_MAT_V=DIST_MAT;

n_iter_update_V=iter_V;

      results_V(m,1)=n_iter_update_V;
      results_V(m,2)=t_update_V;
      results_V(m,3)=(results_V(m,1)<ITER_MAX);

[s_jt_predict,~]=...
  share_func(delta_updated+mu_ijt_est,zeros(1,ns,1,T),rho_est,weight);%J*1*G*T
DIST_s_jt=max(abs(log(s_jt_predict(:))-log(S_jt_data(:))));
results_V(m,4)=log10(DIST_s_jt);
results_V(m,5)=(results_V(m,4)<log10(TOL));

     
%% V_update_func spectral
tic
for kk=1:n_sim
    if GPU_spec==1
        output_spectral=...
        spectral_func(@V_update_func,1,V_initial0,...
        [],[],...
        weight,mu_ijt_est,rho_est,...
        S_jt_data,S_0t_data,numer_1_without_delta,...
        beta_C,L,tune_param,Newton_spec);
    else
        output_spectral=...
        spectral_func_cpu(@V_update_func,1,V_initial0,...
        [],[],...
        weight,mu_ijt_est,rho_est,...
        S_jt_data,S_0t_data,numer_1_without_delta,...
        beta_C,L,tune_param,Newton_spec);
    end

    V_sol=output_spectral{1};

end%%kk
t_update_V_spectral=toc/n_sim;


delta_sol=compute_delta_from_V_func(mu_ijt_est,weight,...
        S_jt_data,rho_est,...
        V_sol);%J*1

n_iter_update_V_spectral=count;
DIST_MAT_V_spectral=DIST_table;
ratio_delta_V_spectral=delta_sol./delta_jt_true;
%[min(ratio_delta_V_spectral),max(ratio_delta_V_spectral)]

results_V_spectral(m,1)=n_iter_update_V_spectral;
results_V_spectral(m,2)=t_update_V_spectral;
results_V_spectral(m,3)=(results_V_spectral(m,1)<ITER_MAX);

[s_jt_predict,~]=...
  share_func(delta_sol+mu_ijt_est,zeros(1,ns,1,T),rho_est,weight);%J*1*G*T
DIST_s_jt=max(abs(log(s_jt_predict(:))-log(S_jt_data(:))));
results_V_spectral(m,4)=log10(DIST_s_jt);
results_V_spectral(m,5)=(results_V_spectral(m,4)<log10(TOL));

if mistake_spec==0
    results_V_spectral(m,3)=(results_V_spectral(m,1)<ITER_MAX & ...
        max(abs(ratio_delta_V_spectral(:)-1))<1e-8);
end

if tune_param==0
    results_V_contraction=results_V;
    results_V_contraction_spectral=results_V_spectral;
else
    results_V_new=results_V;
    results_V_new_spectral=results_V_spectral;
end
       