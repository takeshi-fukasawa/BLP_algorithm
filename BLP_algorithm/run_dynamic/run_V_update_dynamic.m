

if tune_param==0
    results_V=results_V_0;
    results_V_spectral=results_V_0_spectral;
elseif tune_param==1
    results_V=results_V_1;
    results_V_spectral=results_V_1_spectral;
else
    results_V=results_V_2;
    results_V_spectral=results_V_2_spectral;
end

%% V_update_func
%%%V_initial0=V_true;
V_initial=V_initial0;
%V_initial=V_true;

DIST_MAT=zeros(ITER_MAX,1);
tic
for k=1:ITER_MAX

    [output,other_vars]=...
    V_update_func(...
    V_initial,weight,mu_ijt_est,...
    S_jt_data,S_0t_data,weight_V,x_V,beta_C,tune_param,Newton_spec);
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
        V_updated(:,:,:,:,1));%J*1

IV_state=other_vars.IV;

ratio_delta_V=delta_updated./delta_jt_true;

n_iter_update_V=k;

EV=compute_EV_func(V_updated,IV_state,weight_V,x_V);
[s_jt_predict,~]=...
  share_func(delta_updated+mu_ijt_est,beta_C*EV(:,:,:,:,1),rho_est,weight);%J*1*G*T
DIST_s_jt=max(abs(log(s_jt_predict(:))-log(S_jt_data(:))));

results_V(m,1)=n_iter_update_V;
results_V(m,2)=t_update_V;
results_V(m,3)=(results_V(m,1)<ITER_MAX);
results_V(m,4)=log10(DIST_s_jt);
results_V(m,5)=(results_V(m,4)<log10(TOL_DIST_s_jt));

%%%%%%%%%%%%%%%%%%%%%
%% V_update_func spectral
tic
for kk=1:n_sim
        [output_spectral,other_vars,DIST_table_spectral]=...
        spectral_func(@V_update_func,1,0,t_dim_id,V_initial0,...
        weight,mu_ijt_est,...
        S_jt_data,S_0t_data,weight_V,x_V,...
        beta_C,tune_param,Newton_spec);

    V_sol=output_spectral{1};


end%%kk
t_update_V_spectral=toc/n_sim;

IV_state=other_vars.IV;

delta_sol=compute_delta_from_V_func(...
        mu_ijt_est,weight,S_jt_data,rho_est,V_sol);

EV=compute_EV_func(V_sol,IV_state,weight_V,x_V);
[s_jt_predict,~]=...
  share_func(delta_sol+mu_ijt_est,beta_C*EV(:,:,:,:,1),rho_est,weight);%J*1*G*T
DIST_s_jt_spectral=max(abs(log(s_jt_predict(:))-log(S_jt_data(:))));


n_iter_update_V_spectral=count;
DIST_MAT_V_spectral=DIST_table;
ratio_delta_V_spectral=delta_sol./delta_jt_true;
%[min(ratio_delta_V_spectral),max(ratio_delta_V_spectral)]

results_V_spectral(m,1)=n_iter_update_V_spectral;
results_V_spectral(m,2)=t_update_V_spectral;
results_V_spectral(m,3)=(results_V_spectral(m,1)<ITER_MAX);
results_V_spectral(m,4)=log10(DIST_s_jt_spectral);
results_V_spectral(m,5)=(results_V_spectral(m,4)<log10(TOL_DIST_s_jt));


if tune_param==0
    results_V_0=results_V;
    results_V_0_spectral=results_V_spectral;
elseif tune_param==1
    results_V_1=results_V;
    results_V_1_spectral=results_V_spectral;
else
    results_V_2=results_V;
    results_V_2_spectral=results_V_spectral;
end
