clear
global DEBUG FLAG_ERROR DIST DIST_table count ITER_MAX TOL

addpath('./functions')
addpath('./run_static')

large_hetero_spec=0;
mistake_spec=0;
tune_param=1;
tune_param_BLP=0;
durable_spec=0;


n_sim=1;
I=1000;
%I=100;

n_market=1;
T=1;

beta_C=0.0;
L=1;

if 1==0
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

%% Simulation 3
J=25;% Number of products per nest
G=1;
beta_0=0;
rho_true=0.0;
rho_est=rho_true;%%%%%

run run_RCNL_iterations_static.m
results_3=results;

%% Simulation 4
J=25;% Number of products per nest
G=3;
beta_0=0;
rho_true=0.8;
rho_est=rho_true;%%%%%

run run_RCNL_iterations_static.m
results_4=results;

end

%% Simulation 5
J=25;% Number of products per nest
G=3;
beta_0=-4;
rho_true=0.95;
rho_est=rho_true;%%%%%

run run_RCNL_iterations_static.m
results_5=results;

results_no_nest=[results_1;results_2];
results_nest=[results_3;results_4;results_5];

results_no_nest(:,end)=results_no_nest(:,end)*100;
results_nest(:,end)=results_nest(:,end)*100;

if 1==0
writematrix(results_nest,append("C:\Users\fukas\Dropbox\light_bulb\simulation_data\Monte_Carlo\nest_results_",string(mistake_spec),".csv"))

writematrix(results_no_nest,append("C:\Users\fukas\Dropbox\light_bulb\simulation_data\Monte_Carlo\no_nest_results_",string(mistake_spec),".csv"))
end
