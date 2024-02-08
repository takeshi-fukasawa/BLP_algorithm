clear
global DEBUG FLAG_ERROR DIST DIST_table count ITER_MAX TOL

addpath('./functions')
addpath('./run_static')

mistake_spec=0;
large_hetero_spec=1;
durable_spec=0;

n_draw=1;
n_sim=1;

%%I=100;
n_dim_V=1;

n_market=1;
T=1;

beta_C=0.0;
L=1;
G=1;
rho_true=0.0;
rho_est=rho_true;%%%%%

%% Simulation 1
I=2;
J=2;% Number of products per nest
mu_base=0;
f_hetero=10;
beta_0=-1;

run run_RCNL_iterations_static.m
results_1=results;

u_ijt_1=u_ijt_tilde;

BLP_paper_figure_path="C:/Users/fukas/Dropbox/アプリ/Overleaf/BLP/figure/";
temp=(log10(DISTMAT_BLP_1(:,1)));
plot(temp)
saveas(gcf,append(BLP_paper_figure_path,'DIST_large_hetero.png'))


plot(diff(temp))
ylim([-0.001,0.001])
saveas(gcf,append(BLP_paper_figure_path,'DIST_diff_large_hetero.png'))

CCP_table=round([s_ijt_ccp_true;s_i0t_ccp_true],4);
filename=append(BLP_paper_figure_path,'Large_hetero_CCP.csv');
writematrix(CCP_table,filename)

filename=append(BLP_paper_figure_path,'Large_hetero_results.csv');
writematrix(results_1,filename)


if 1==0

%%%%%%%%%%%%%
%% Simulation 1_temp
I=2;
J=2;% Number of products per nest

%%%%%%%%%%%%%
mu_base=-10;
mu_base=0;
f_hetero=10;
beta_0=0;
%%%%%%%%%%%

run run_RCNL_iterations_static.m
results_1_temp=results;

u_ijt_1_temp=u_ijt_tilde;

%%%%%%%%%%%

temp=(1-s_i0t_ccp_true).*abs(prob_i_given_j-prob_i_given_0);%J*ns
conv_const_delta_new=max(sum(temp,2))

temp=(1-s_i0t_ccp_true).*abs(prob_i_given_j);%J*ns
conv_const_delta_old=max(sum(temp,2))

%%% G==1 case
temp1=sum(s_ijt_given_g_ccp.*(prob_i_given_j-prob_i_given_0),2);%J*1
modulus=sum(abs(temp1),1)

temp1_temp=sum(abs(s_ijt_given_g_ccp.*(prob_i_given_j-prob_i_given_0)),2);%J*1
modulus_temp=sum(abs(temp1_temp),1);

%plot(diff(log10(DISTMAT_BLP_1)))

%% Simulation 2
I=2;
J=10;% Number of products per nest
f_hetero=10;
beta_0=-1;

run run_RCNL_iterations_static.m
results_2=results;

%% Simulation 3
I=100;
J=10;% Number of products per nest
f_hetero=10;
beta_0=20;

run run_RCNL_iterations_static.m
results_3=results;

end
