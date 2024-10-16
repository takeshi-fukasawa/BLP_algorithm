
global feval_Bellman
hot_start_V_spec=1;

%% V_update_func
for method=2:2
    run spec_settings.m
    feval_Bellman=0;

    tStart=tic;
        [output_spectral,other_vars,iter_info]=...
        spectral_func(@delta_middle_Bellman_inner_func,...
        spec,{delta_initial0},V_initial0,...
        mu_ijt_est,S_jt_data,beta_C,rho_est,...
        weight,tune_param_BLP,spec,...
        weight_V,x_V,hot_start_V_spec);
    iter_info.t_cpu=toc(tStart);

    %%iter_info.t_cpu=tEnd-tStart; % inner and middle loop spectral algorithm => False toc...
    iter_info.feval_Bellman=feval_Bellman;
    delta_sol=output_spectral{1};

    s_jt_predict=other_vars.s_jt_predict;

    results_BLP_middle=results_output_func(iter_info,s_jt_predict,S_jt_data);
    ratio_delta=delta_sol./delta_jt_true;

if method==1

    if tune_param==0
        results_BLP_middle_0(m,:)=results_BLP_middle;
        ratio_delta_BLP_middle_0=ratio_delta;
        iter_info_BLP_middle_0=iter_info;
    elseif tune_param==1
        results_BLP_middle_1(m,:)=results_BLP_middle;
        ratio_delta_BLP_middle_1=ratio_delta;
        iter_info_V_1=iter_info;
    elseif tune_param>1
        results_BLP_middle_2(m,:)=results_BLP_middle;
        ratio_delta_BLP_middle_2=ratio_delta;
        iter_info_BLP_middle_2=iter_info;
    end

elseif method==2
    if tune_param==0
        results_BLP_middle_0_spectral(m,:)=results_BLP_middle;
        ratio_delta_BLP_middle_0_spectral=ratio_delta;
        iter_info_BLP_middle_0_spectral=iter_info;
    elseif tune_param==1
        results_BLP_middle_1_spectral(m,:)=results_BLP_middle;
        ratio_delta_BLP_middle_1_spectral=ratio_delta;
        iter_info_BLP_middle_1_spectral=iter_info;
    elseif tune_param>1
        results_BLP_middle_2_spectral(m,:)=results_BLP_middle;
        ratio_delta_BLP_middle_2_spectral=ratio_delta;
        iter_info_BLP_middle_2_spectral=iter_info;
    end% tune_param==0 or 1 or others

elseif method==3
    if tune_param==0
        results_BLP_middle_0_SQUAREM(m,:)=results_BLP_middle;
        ratio_delta_BLP_middle_0_SQUAREM=ratio_delta;
        iter_info_BLP_middle_0_SQUAREM=iter_info;
    elseif tune_param==1
        results_BLP_middle_1_SQUAREM(m,:)=results_BLP_middle;
        ratio_delta_BLP_middle_1_SQUAREM=ratio_delta;
        iter_info_BLP_middle_1_SQUAREM=iter_info;
    elseif tune_param>1
        results_BLP_middle_2_SQUAREM(m,:)=results_BLP_middle;
        ratio_delta_BLP_middle_2_SQUAREM=ratio_delta;
        iter_info_BLP_middle_2_SQUAREM=iter_info;
    end% tune_param==0 or 1 or others

elseif method==4
    if tune_param==0
        results_BLP_middle_0_Anderson(m,:)=results_BLP_middle;
        ratio_delta_BLP_middle_0_Anderson=ratio_delta;
        iter_info_BLP_middle_0_Anderson=iter_info;
    elseif tune_param==1
        results_BLP_middle_1_Anderson(m,:)=results_BLP_middle;
        ratio_delta_BLP_middle_1_Anderson=ratio_delta;
        iter_info_BLP_middle_1_Anderson=iter_info;
    elseif tune_param>1
        results_BLP_middle_2_Anderson(m,:)=results_BLP_middle;
        ratio_delta_BLP_middle_2_Anderson=ratio_delta;
        iter_info_BLP_middle_2_Anderson=iter_info;
    end% tune_param==0 or 1 or others

end


end % method

results_V_delta_nested=[feval_Bellman,results_BLP_middle];
results_V_delta_joint=[results_V_BLP_1_spectral(end,1),results_V_BLP_1_spectral(end,:)];

results_comparison=[results_V_delta_joint;results_V_delta_nested];

if 1==0 & IVS_spec==1
    filename=append(output_path,"dynamic_BLP_IVS_results_beta_",...
        string(beta_C),"_",string(mistake_spec),"_V_delta_mapping_comparison_IVS.csv");
    writematrix(results_comparison,filename)
end

if 1==0 & IVS_spec==0
    filename=append(output_path,"dynamic_BLP_IVS_results_beta_",...
        string(beta_C),"_",string(mistake_spec),"_V_delta_mapping_comparison_perfect_foresight.csv");
    writematrix(results_comparison,filename)
end
