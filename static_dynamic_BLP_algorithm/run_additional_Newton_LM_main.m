%% Additional results (Try Newton/LM)
rho_true=0.0;
rho_est=rho_true;%%%%%
G=1;

    %% Simulation 1
    J=25;% Number of products per nest
    beta_0=0;
    iniital_scale=2;
    far_initial_delta_spec=0;


    %% Simulation 1 (Base)
    run run_RCNL_iterations_static_additional.m
    results_no_nest_additional_1=table;

    t_mapping_vec_1=zeros(3,1);
    for mapping_spec=1:2
        [t_base,t_additional]=check_mapping_comp_cost...
            (mapping_spec,delta_initial0,mu_ijt_est,...
        weight,S_jt_data,S_0t_data);
        if mapping_spec==1
            t_mapping_vec_1(1:2)=[t_base;t_base+t_additional];
        else% Newton
            t_mapping_vec_1(3)=t_base;
        end
    end

    %% Simulation 1 (Small outside share)
    beta_0=4;
    initial_scale=2;
    far_initial_delta_spec=0;


    run run_RCNL_iterations_static_additional.m
    results_no_nest_additional_1_small_outside_share=table;

    %% Simulation 1 (Large hetero)
    beta_0=0;
    initial_scale=10;
    far_initial_delta_spec=0;

    run run_RCNL_iterations_static_additional.m
    results_no_nest_additional_1_large_hetero=table;

    %% Simulation 1 (Far initial values)
    beta_0=0;
    initial_scale=2;
    far_initial_delta_spec=1;

    run run_RCNL_iterations_static_additional.m
    results_no_nest_additional_1_far_initial_value=table;

    %% Simulation 2
    J=250;% Number of products per nest
    beta_0=0;
    initial_scale=2;
    far_initial_delta_spec=0;

    run run_RCNL_iterations_static_additional.m
    results_no_nest_additional_2=table;

    t_mapping_vec_2=zeros(3,1);
    for mapping_spec=1:2
        [t_base,t_additional]=check_mapping_comp_cost...
            (mapping_spec,delta_initial0,mu_ijt_est,...
        weight,S_jt_data,S_0t_data);
        if mapping_spec==1
            t_mapping_vec_2(1:2)=[t_base;t_base+t_additional];
        else% Newton
            t_mapping_vec_2(3)=t_base;
        end
    end

    t_mapping_vec=round([t_mapping_vec_1,t_mapping_vec_2],5);

    results_additional=[results_no_nest_additional_1;...
        results_no_nest_additional_1_small_outside_share;...
        results_no_nest_additional_1_large_hetero;...
        results_no_nest_additional_1_far_initial_value;...
        results_no_nest_additional_2];

    if n_market>1 & 1==0
        filename=append(output_path,"results_RCL_mapping_comp_cost.csv");
        writematrix(t_mapping_vec,filename)
        
        results_additional(:,end)=results_additional(:,end)*100;%Percentage
        filename=append(output_path,"results_additional_newton_etc.csv");
        writematrix(round(results_additional,5),filename)

    end

    %% Case with small number of consumer types
    I=10;
    ns=I;

    run run_RCNL_iterations_static_additional.m
    results_no_nest_additional_small_number_of_consumer_types=table;

    t_mapping_vec_small_I=zeros(3,1);
    for mapping_spec=1:2
        [t_base,t_additional]=check_mapping_comp_cost...
            (mapping_spec,delta_initial0,mu_ijt_est,...
        weight,S_jt_data,S_0t_data);
        if mapping_spec==1
            t_mapping_vec_small_I(1:2)=[t_base;t_base+t_additional];
        else% Newton
            t_mapping_vec_small_I(3)=t_base;
        end
    end
