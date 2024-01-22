clear
global DEBUG FLAG_ERROR DIST DIST_table count ITER_MAX TOL

addpath('./functions')
addpath('./run_static')

mistake_spec=0;

GPU_spec=0;

n_sim=1;
I=2;
n_dim_V=1;

n_market=3;
T=1;

beta_C=0.0;
L=1;

%% Simulation 1
J=2;% Number of products per nest
G=1;
f_hetero=10;

beta_0=15;
rho_true=0.0;
rho_est=rho_true;%%%%%

run run_RCNL_iterations.m
results_1=results;

temp=(1-s_i0t_ccp_true).*abs(prob_i_given_j-prob_i_given_0);%J*ns
conv_const_delta_new=max(sum(temp,2))

temp=(1-s_i0t_ccp_true).*abs(prob_i_given_j);%J*ns
conv_const_delta_old=max(sum(temp,2))

%%% G==1 case
temp1=sum(s_ijt_given_g_ccp.*(prob_i_given_j-prob_i_given_0),2);%J*1
modulus=sum(abs(temp1),1)

temp1_temp=sum(abs(s_ijt_given_g_ccp.*(prob_i_given_j-prob_i_given_0)),2);%J*1
modulus_temp=sum(abs(temp1_temp),1);


