clear

addpath('./functions')
addpath('./run_static')
addpath('./run_static_dynamic')

addpath('C:/Users/fukas/Dropbox/git/spectral')

mistake_spec=0;
large_hetero_spec=0;
durable_spec=0;
t_dependent_alpha_spec=0;

spec_default=[];
spec_default.compute_alpha_spec=3;

n_sim=1;
I=10;
%I=1000;
n_dim_V=1;
n_draw=1;

n_market=1;
T=1;

beta_C=0.0;
L=1;

%% Simulation 1
J=25;% Number of products per nest
G=1;
beta_0=0;
rho_true=0.0;
rho_est=rho_true;%%%%%

run run_RCNL_iterations_static.m
results_1=results;

%% Simulation 2
J=25;% Number of products per nest
G=1;
beta_0=4;
rho_true=0.0;
rho_est=rho_true;%%%%%

run run_RCNL_iterations_static.m
results_2=results;


results_no_nest=[results_1;results_2];

results_no_nest(:,end)=results_no_nest(:,end)*100;

if 1==0
    filename=append("C:/Users/fukas/Dropbox/light_bulb/simulation_data/Monte_Carlo/no_nest_results_ns_",...
        string(ns),"_",string(mistake_spec),".csv");

    writematrix(results_no_nest,filename)
end

%% Simulation 3 (RCNL)
J=25;% Number of products per nest
G=2;
beta_0=4;
rho_true=0.5;
rho_est=rho_true;%%%%%

run run_RCNL_iterations_static.m
results_3=results;
