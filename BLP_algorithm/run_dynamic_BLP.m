%%%%%%%%%%%%%%%%%
%% Parameter settings
clear
global DEBUG FLAG_ERROR DIST DIST_table count ITER_MAX TOL
global R2 y_mat
global coef_1_true coef_0_true
global Pr0_spec

coef_0_true=[];
coef_1_true=[];

addpath('./functions')
addpath('./run_dynamic')
addpath('C:/Users/fukas/Dropbox/git/spectral')
save_path="C:/Users/fukas/Dropbox/light_bulb/simulation_data/Monte_Carlo/";

large_hetero_spec=0;
mistake_spec=0;

rho_true=0;
rho_est=rho_true;
L=100000;

G=1;

n_sim=1;
I=1000;
I=2;

ns=I;
TOL=1e-13;%% Important??
ITER_MAX=3000;
%ITER_MAX=2000;

n_market=1;
n_draw=1;

%%%%%%%%%%%%%%%%
%% Simulation 1
durable_spec=1;
%%durable_spec=0;

Pr0_spec=1; %%%% Pr0_spec==1=> Introduce endogenous Pr0

t_dependent_alpha_spec=1;

n_dim_V=20;
%n_dim_V=1;

J=25;% Number of products per nest
beta_0=4;

beta_C=0.99;

%%% T==1=> stationary;
T=100;
T=30;
%T=2;
%T=1;

run run_RCNL_iterations_dynamic.m
