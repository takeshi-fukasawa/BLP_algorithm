%%%%%%%%%%%%%%%%%
%% Parameter settings
clear
global DEBUG FLAG_ERROR DIST DIST_table count ITER_MAX TOL

addpath('./functions')
addpath('./run_dynamic')

large_hetero_spec=0;
mistake_spec=0;

rho_true=0;
rho_est=rho_true;
L=100000;

G=1;

GPU_spec=0;

n_sim=1;
I=1000;
I=2;

ns=I;
TOL=1e-14;%% Important??
ITER_MAX=1000;

n_market=1;
%%%%%%%%%%%%%%%%
%% Simulation 1
durable_spec=1;

n_dim_V=6;

J=25;% Number of products per nest
beta_0=4;

beta_C=0.99;

%%% T==1=> stationary;
T=10;

run run_RCNL_iterations_dynamic.m
