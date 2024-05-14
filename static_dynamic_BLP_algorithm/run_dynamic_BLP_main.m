
%%% Run Monte Carlo simulation of dynamic BLP inner-loop algorithms.
%%% Written by Takeshi Fukasawa in May 2024.

clear

global Pr0_spec inv_multiply_Chebyshev_basis Chebyshev_extrema
global ratio_s_i0t_ccp_t

rng('default')

run preliminary_path.m

addpath('./functions')
addpath('./run_dynamic')
addpath('./run_static_dynamic')
addpath(spectral_func_path)


spec_default=[];
spec_default.compute_alpha_spec=3;
spec_default.TOL=1e-12;
spec_default.common_alpha_spec=0;
spec_default.ITER_MAX=3000;
%spec_default.ITER_MAX=3;
%%spec_default.stabilize_spec=1;

skip_contraction_spec=0;
large_hetero_spec=0;
mistake_spec=1;%%%%

rho_true=0;
rho_est=rho_true;
L=100000;

G=1;

n_sim=1;
I=1000;
I=50;
%I=2;

ns=I;

n_market=20;
n_draw=1;

IVS_spec=0; % If IVS_spec==1, use inclusive value sufficiency. If IVS_spec==0, assume perfect foresight on the future state transitions of market-level variables. For details, see the paper.

%%%%%%%%%%%%%%%%
%% Simulation 1
durable_spec=1; % Use the model in Sun and Ishihara (2019)
%%durable_spec=0;

Pr0_spec=1; %%%% Pr0_spec==1=> Introduce endogenous Pr0

t_dependent_alpha_spec=0; % Introduce time-dependent step sizes alpha in the spectral / SQUAREM algorithms.

n_dim_V=1;

J=25;% Number of products per nest

beta_C=0.99;
%%beta_C=0;

%%% T==1=> stationary;
T=100;
T=30;
%T=2;
%T=1;

if IVS_spec==1
   T=25;
   n_grid_IV=10;
   beta_0=6;
   Pr0_spec=1;
   t_dependent_alpha_spec=0;%%%
   %%spec_default.alpha_max=5;%%%%
   spec_default.ITER_MAX=3000;
   n_draw=5;
else
    T=50;
    n_grid_IV=0;
    beta_0=6;
    Pr0_spec=1;
    t_dependent_alpha_spec=1;
    spec_default.alpha_max=10;
end

%%T=1;

run run_RCNL_iterations_dynamic.m

if IVS_spec==1
    filename=append(output_path,"dynamic_BLP_IVS_results_beta_",...
        string(beta_C),"_",string(mistake_spec),".csv");
else
    filename=append(output_path,"dynamic_BLP_nonstationary_results_beta_",...
        string(beta_C),"_",string(mistake_spec),".csv");
end

if 1==0
    writematrix(results_table,filename)
end

