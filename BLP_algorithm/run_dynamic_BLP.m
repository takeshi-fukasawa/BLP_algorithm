%%%%%%%%%%%%%%%%%
%% Parameter settings
clear
global DEBUG FLAG_ERROR DIST DIST_table count ITER_MAX TOL
global I_ig
global delta_updated

addpath('./functions')
addpath('./run_dynamic')

large_hetero_spec=0;
mistake_spec=0;

n_dim_V=1;

J=25;% Number of products per nest
G=1;
beta_0=4;

beta_C=0.99;
beta_C=0.99;

rho_true=0;
rho_est=rho_true;
L=100000;


GPU_spec=0;

n_sim=1;
I=1000;
I=10;

T=10;

ns=I;
TOL=1e-12;
ITER_MAX=3000;

n_market=1;

for m=1:n_market

rng(m);
gpurng(m);

%% Generate data and solve equilibrium
%run DGP.m
run DGP_ABLP.m
mu_ij_est=mu_ijt_true*1;

run solve_equil.m

V_initial0=-log(S_0t_data.*weight);
delta_initial0=log(S_jt_data)-log(S_0t_data)-rho_est.*log(S_jt_given_g_data);% Initial value of delta

%% Update V (rho==0 & G==1 case only)
if G==1 & rho_est==0
    %%% tune_param==0
tune_param=0;Newton_spec=0;
run run_V_update_dynamic.m

n_iter_update_V_0(m,1)=n_iter_update_V;
t_update_V_0(m,1)=t_update_V;
n_iter_update_V_spectral_0(m,1)=n_iter_update_V_spectral;
t_update_V_spectral_0=t_update_V_spectral;
ratio_delta_V_0=ratio_delta_V;
ratio_delta_V_spectral_0=ratio_delta_V_spectral;

%%% tune_param==1
tune_param=1;Newton_spec=0;

run run_V_update_dynamic.m

n_iter_update_V_1(m,1)=n_iter_update_V;
t_update_V_1(m,1)=t_update_V;
n_iter_update_V_spectral_1(m,1)=n_iter_update_V_spectral;
t_update_V_spectral_1=t_update_V_spectral;
ratio_delta_V_1=ratio_delta_V;
ratio_delta_V_spectral_1=ratio_delta_V_spectral;

%%% tune_param==1/(1-beta_C)
tune_param=1/(1-beta_C);Newton_spec=0;
run run_V_update_dynamic.m

n_iter_update_V_2(m,1)=n_iter_update_V;
t_update_V_2(m,1)=t_update_V;
n_iter_update_V_spectral_2(m,1)=n_iter_update_V_spectral;
t_update_V_spectral_2=t_update_V_spectral;
ratio_delta_V_2=ratio_delta_V;
ratio_delta_V_spectral_2=ratio_delta_V_spectral;

%%% Newton iteration
if 1==0
tune_param=1;Newton_spec=1;
run run_V_update_dynamic.m

n_iter_update_V_Newton(m,1)=n_iter_update_V;
t_update_V_Newton(m,1)=t_update_V;
n_iter_update_V_spectral_Newton(m,1)=n_iter_update_V_spectral;
t_update_V_spectral_Newton=t_update_V_spectral;
end

end


%% Jointly Update V and IV
if 1==1
run run_V_IV_update_dynamic.m
end

%% Joint update of delta and V
%%% tune_param_BLP==0
tune_param_BLP=0;
run run_delta_V_joint_update_dynamic.m

n_iter_update_V_BLP_0(m,1)=n_iter_update_V_BLP;
t_update_V_BLP_0(m,1)=t_BLP_V;
n_iter_update_V_BLP_spectral_0(m,1)=n_iter_update_V_BLP_spectral; 
t_update_V_BLP_spectral_0(m,1)=t_V_BLP_spectral;

%%% tune_param_BLP==1
tune_param_BLP=1;
run run_delta_V_joint_update_dynamic.m

n_iter_update_V_BLP_1(m,1)=n_iter_update_V_BLP;
t_update_V_BLP_1(m,1)=t_BLP_V;
n_iter_update_V_BLP_spectral_1(m,1)=n_iter_update_V_BLP_spectral; 
t_update_V_BLP_spectral_1(m,1)=t_V_BLP_spectral;

%% delta fixed convergence
delta_fixed_spec=1;
run run_Bellman_middle_update_dynamic.m

n_iter_delta_fixed(m,1)=n_iter_update_V_middle;
t_update_delta_fixed(m,1)=t_V_middle;

n_iter_delta_fixed_spectral(m,1)=n_iter_update_V_middle_spectral;
t_update_delta_fixed_spectral(m,1)=t_V_middle_spectral;

end


if 1==0
    run run_BLP_middle_update_dynamic.m
end

if G==1
results=[...
    mean(n_iter_update_V_BLP_0,1) mean(t_update_V_BLP_0,1);...
    mean(n_iter_update_V_BLP_spectral_0,1) mean(t_update_V_BLP_spectral_0,1);...
    mean(n_iter_update_V_BLP_1,1) mean(t_update_V_BLP_1,1);...
    mean(n_iter_update_V_BLP_spectral_1,1) mean(t_update_V_BLP_spectral_1);...
    mean(n_iter_update_V_0,1) mean(t_update_V_0,1);...
    mean(n_iter_update_V_spectral_0,1) mean(t_update_V_spectral_0,1);...
    mean(n_iter_update_V_1,1) mean(t_update_V_1,1);...
    mean(n_iter_update_V_spectral_1,1) mean(t_update_V_spectral_1,1);...
    mean(n_iter_update_V_2,1) mean(t_update_V_2,1);...
    mean(n_iter_update_V_spectral_2,1) mean(t_update_V_spectral_2,1);...
    %mean(n_iter_update_V_Newton,1) mean(t_update_V_Newton,1);...
    mean(n_iter_delta_fixed,1) mean(t_update_delta_fixed,1);...
    mean(n_iter_delta_fixed_spectral,1) mean(t_update_delta_fixed_spectral,1);...
    ];

results(:,2)=round(results(:,2),3);
results=[zeros(size(results,1),2),results];
results(:,1)=beta_C;
results(:,2)=median(S_0t_data(:));

if 1==0
    filename=append("C:/Users/fukas/Dropbox/light_bulb/simulation_data/Monte_Carlo/no_nest_results_beta_",...
        string(beta_C),"_",string(mistake_spec),".csv");

writematrix(round(results,3),filename)
end

end
