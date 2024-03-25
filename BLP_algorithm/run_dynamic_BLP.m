%%%%%%%%%%%%%%%%%
%% Parameter settings
clear

global Pr0_spec inv_multiply_Chebyshev_basis Chebyshev_extrema
rng('default')

addpath('./functions')
addpath('./run_dynamic')
addpath('./run_static_dynamic')
addpath('C:/Users/fukas/Dropbox/git/spectral')
save_path="C:/Users/fukas/Dropbox/light_bulb/simulation_data/Monte_Carlo/";


spec_default=[];
spec_default.compute_alpha_spec=3;
spec_default.TOL=1e-12;
spec_default.common_alpha_spec=0;
spec_default.ITER_MAX=3000;
%spec_default.ITER_MAX=3;

skip_contraction_spec=1;
large_hetero_spec=0;
mistake_spec=0;

rho_true=0;
rho_est=rho_true;
L=100000;

G=1;

n_sim=1;
I=1000;
I=50;
I=2;

ns=I;

n_market=1;
n_draw=1;

Sun_Ishihara_spec=1;

%%%%%%%%%%%%%%%%
%% Simulation 1
durable_spec=1;
%%durable_spec=0;

Pr0_spec=1; %%%% Pr0_spec==1=> Introduce endogenous Pr0

t_dependent_alpha_spec=0;

n_dim_V=1;

J=25;% Number of products per nest
beta_0=4;

beta_C=0.99;
%%beta_C=0;

%%% T==1=> stationary;
T=100;
T=30;
%T=2;
%T=1;

if Sun_Ishihara_spec==1
   T=25;
   n_grid_IV=10;
   beta_0=6;
   Pr0_spec=1;
   durable_spec=1;
   t_dependent_alpha_spec=0;
else
    T=50;
    n_grid_IV=0;
    beta_0=6;
    Pr0_spec=1;
    durable_spec=1;
    t_dependent_alpha_spec=1;
end

run run_RCNL_iterations_dynamic.m
