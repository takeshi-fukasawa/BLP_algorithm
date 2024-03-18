
%% IV_update_func 
for method=1:2
    spec=spec_default;
    if method==1 % fixed point iteration
        spec.update_spec=0;
    elseif method==2 % spectral
    if t_dependent_alpha_spec==1
        spec.update_spec=t_dim_id;
    else
        spec.update_spec=[];%%%%%%
    end
end

IV_initial0=repmat(log(S_gt_data)-log(S_0t_data),[1,ns,1,1]);%1*ns*G*T

[output_spectral,other_vars,iter_info]=...
        spectral_func(@IV_update_func,spec,...
        {IV_initial0},...
        weight,mu_ijt_est,rho_est,...
        S_jt_data,S_0t_data,tune_param);


    delta_sol=compute_delta_from_V_IV_func(mu_ijt_est,weight,...
        S_jt_data,rho_est,...
        [],IV_sol);%J*1


    [s_jt_predict,~]=...
  share_func(delta_sol+mu_ijt_est,zeros(1,ns,1,T),rho_est,weight);%J*1*G*T

    results_IV=results_output_func(iter_info,s_jt_predict,S_jt_data);
    ratio_delta=delta_sol./delta_jt_true;


if method==1

if tune_param==0
    results_IV_0(m,:)=results_V;
    ratio_delta_IV_0=ratio_delta;
    iter_info_IV_0=iter_info;
elseif tune_param==1
    results_IV_1(m,:)=results_IV;
    ratio_delta_IV_1=ratio_delta;
    iter_info_IV_1=iter_info;
elseif tune_param>1
    results_IV_2(m,:)=results_IV;
    ratio_delta_IV_2=ratio_delta;
    iter_info_IV_2=iter_info;
end

elseif method==2

if tune_param==0
    results_IV_0_spectral(m,:)=results_IV;
    ratio_delta_IV_0_spectral=ratio_delta;
    iter_info_IV_0_spectral=iter_info;
elseif tune_param==1
    results_IV_1_spectral(m,:)=results_IV;
    ratio_delta_IV_1_spectral=ratio_delta;
    iter_info_IV_1_spectral=iter_info;
elseif tune_param>1
    results_IV_2_spectral(m,:)=results_IV;
    ratio_delta_IV_2_spectral=ratio_delta;
    iter_info_IV_2_spectral=iter_info;
end% tune_param==0 or 1 or others
end

end % method
