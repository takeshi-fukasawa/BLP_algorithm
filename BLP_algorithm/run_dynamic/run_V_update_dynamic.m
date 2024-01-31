
%% V_update_func
%%%V_initial0=V_true;
V_initial=V_initial0;
%V_initial=V_true;

DIST_MAT=zeros(ITER_MAX,1);
tic
for k=1:ITER_MAX

    output=...
    V_update_func(...
    V_initial,weight,mu_ijt_est,rho_est,...
    S_jt_data,S_0t_data,numer_1_without_delta,beta_C,L,tune_param,Newton_spec);
    resid=output{1};

    V_updated=V_initial-resid;

    if 1==0 % scaling; effective for Newton* iteration??
        s_i0t_ccp_temp=exp((beta_C-1).*V_updated);
        s_i0t_ccp=S_0_data./sum(weight.*s_i0t_ccp_temp,2).*s_i0t_ccp_temp;
        s_i0t=s_i0t_ccp.*weight;
        V_updated=log(s_i0t_ccp)./(beta_C-1);
    end

    DIST=max(abs(V_updated(:)-V_initial(:)));%scalar
    DIST_MAT(k,:)=DIST;


    %%%%%
    %k
    %DIST
      if DIST<TOL
          break;% end the for loop
      elseif isnan(DIST)==1
          k=ITER_MAX;
          break;
      else
          V_initial=V_updated;%J by 1
      end

end% for loop
t_update_V=toc;

delta_updated=compute_delta_from_V_func(mu_ijt_est,weight,...
        S_jt_data,rho_est,...
        V_updated);%J*1

ratio_delta_V=delta_updated./delta_jt_true;

n_iter_update_V=k;

EV=compute_EV_func(V_updated,[],weight_V);
[s_jt_predict,~]=...
  share_func(delta_updated+mu_ijt_est,beta_C*EV,rho_est,weight);%J*1*G*T
DIST_s_jt=max(abs(log(s_jt_predict(:))-log(S_jt_data(:))));

%%%%%%%%%%%%%%%%%%%%%
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


delta_sol=compute_delta_from_V_func(...
        mu_ijt_est,weight,S_jt_data,rho_est,V_sol);

EV=compute_EV_func(V_sol,[],weight_V);
[s_jt_predict,~]=...
  share_func(delta_sol+mu_ijt_est,beta_C*EV,rho_est,weight);%J*1*G*T
DIST_s_jt_spectral=max(abs(log(s_jt_predict(:))-log(S_jt_data(:))));


n_iter_update_V_spectral=count;
DIST_MAT_V_spectral=DIST_table;
ratio_delta_V_spectral=delta_sol./delta_jt_true;
%[min(ratio_delta_V_spectral),max(ratio_delta_V_spectral)]
