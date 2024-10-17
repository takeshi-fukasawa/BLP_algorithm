%%% Run Monte Carlo simulation of static BLP inner-loop algorithms.
%%% Written by Takeshi Fukasawa in May 2024.


clear

run preliminary_path.m

addpath('./functions')
addpath('./run_static')
addpath('./run_static_dynamic')
addpath('./results')

output_path=append(output_path,'static_BLP/')
addpath(spectral_func_path)

%%% In the estimation, researchers do not know true values of nonlinear parameters, and they solve for mean utility delta  by using candidate values of nonlinear parameters.
%%% We replicate the situation by setting mistake_spec==1. If mistake_spec==0, we solve for delta based on the true nonlinear parameters. 
mistake_spec=1; 

large_hetero_spec=0; % See also run_RCNL_main_hetero_test.m
durable_spec=0;% Used in dynamic BLP
t_dependent_alpha_spec=0; % Used in dynamic BLP
skip_contraction_spec=0; % If 1, try only the spectral/SQUAREM iterations. If 0, also try standard fixed point iterations.

spec_default=[];
spec_default.TOL=1e-13;
spec_default.compute_alpha_spec=3;
spec_default.ITER_MAX=1000;

if spec_default.compute_alpha_spec~=3
    skip_contraction_spec=1;
end

n_sim=1;
I=10;
I=1000;
%I=2;
n_dim_V=1;
n_draw=1;

n_market=5;
T=1;

beta_C=0.0;
L=1;

if 1==1
%% Simulation 1
J=25;% Number of products per nest
G=1;
beta_0=0;
rho_true=0.0;
rho_est=rho_true;%%%%%

run run_RCNL_iterations_static.m
results_1=results;

results_no_nest=[results_1];

    %% Simulation 2
    J=250;% Number of products per nest
    G=1;
    beta_0=0;
    rho_true=0.0;
    rho_est=rho_true;%%%%%
    
    run run_RCNL_iterations_static.m
    results_2=results;
    results_no_nest=[results_1;results_2];



if n_market>1 & 1==1
    filename=append(output_path,"no_nest_results_ns_",...
        string(ns),"_",string(mistake_spec),"_S",string(spec.compute_alpha_spec),".csv");

    writematrix(results_no_nest,filename)
end

end % RCL run

if 1==1 & spec.compute_alpha_spec==3 & I>2
%% Simulation 3 (RCNL)
J=25;% Number of products per nest
G=3;
beta_0=0;
rho_true=0.5;
rho_est=rho_true;%%%%%

run run_RCNL_iterations_static.m
results_RCNL=results;

if n_market>1 & 1==1
    filename=append(output_path,"RCNL_results_ns_",...
        string(ns),"_",string(mistake_spec),"_S",string(spec.compute_alpha_spec),".csv");

    writematrix(results_RCNL,filename)

    filename=append(output_path,"RCNL_table_ns_",...
        string(ns),"_",string(mistake_spec),"_S",string(spec.compute_alpha_spec),".csv");
    
end
end % run RCNL

