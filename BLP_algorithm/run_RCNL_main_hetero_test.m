clear
addpath('./functions')
addpath('./run_static')
addpath('./run_static_dynamic')

addpath('C:/Users/fukas/Dropbox/git/spectral')
addpath('C:/Users/fukas/Dropbox/git/util/uniformFigureStyle/')

BLP_paper_figure_path="C:/Users/fukas/Dropbox/アプリ/Overleaf/BLP/figure/";


spec_default=[];
spec_default.compute_alpha_spec=3;
spec_default.SQUAREM_spec=0;
spec_default.ITER_MAX=2000;
%%spec_default.norm_spec=2;

mistake_spec=0;
large_hetero_spec=1;%%%%
durable_spec=0;
t_dependent_alpha_spec=0;


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

%% Simulation 1 (spectral)
I=2;
J=2;% Number of products per nest
mu_base=0;
f_hetero=10;
beta_0=-1;

run run_RCNL_iterations_static.m
results_1=results;

iter_info_BLP_contraction_spectral0=iter_info_BLP_contraction_spectral;

if large_hetero_spec==1
    %% Noncontraction figure-1
    temp=(log10(iter_info_BLP_new.DIST_table(:,1)));

    clf
    f = figure;
    set(f,'units','centimeters','position',[10,10,12,9])

    plot(temp)
    xlim([0,2000])
    xlabel('Iterations')

    DIST_label='$\log_{10}\left(\left\Vert \delta_{j}^{(n+1)}-\delta_{j}^{(n)}\right\Vert _{\infty}\right)$';
    ylabel(DIST_label,'Interpreter','latex')
   
    uniformFigureStyle(f);
    saveas(gcf,append(BLP_paper_figure_path,'DIST_large_hetero.png'))
   
    %% Non-contraction figure 2
    clf
    f = figure;
    set(f,'units','centimeters','position',[10,10,12,9])

    plot(diff(temp))
    xlim([0,2000])
    xlabel('Iterations')

    ylim([-0.001,0.001])
    
    Delta_DIST_label='$\Delta \left[ \log_{10}\left(\left\Vert \delta_{j}^{(n+1)}-\delta_{j}^{(n)}\right\Vert _{\infty}\right) \right]$';
    ylabel(Delta_DIST_label,'Interpreter','latex')
   
    uniformFigureStyle(f);

    saveas(gcf,append(BLP_paper_figure_path,'DIST_diff_large_hetero.png'))
    
    CCP_table=round([s_ijt_ccp_true;s_i0t_ccp_true],4);
    filename=append(BLP_paper_figure_path,'Large_hetero_CCP.csv');
    writematrix(CCP_table,filename)
    
    filename=append(BLP_paper_figure_path,'Large_hetero_results.csv');
    writematrix(results_1,filename)
end

if 1==1

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

%%%%%%%%%%%

temp=(1-s_i0t_ccp_true).*abs(prob_i_given_j-prob_i_given_0);%J*ns
conv_const_delta_new=max(sum(temp,2))

temp=(1-s_i0t_ccp_true).*abs(prob_i_given_j);%J*ns
conv_const_delta_old=max(sum(temp,2))

%%% G==1 case
temp1=sum(s_ijt_given_g_ccp_true.*(prob_i_given_j-prob_i_given_0),2);%J*1
modulus=sum(abs(temp1),1)

temp1_temp=sum(abs(s_ijt_given_g_ccp_true.*(prob_i_given_j-prob_i_given_0)),2);%J*1
modulus_temp=sum(abs(temp1_temp),1);

%plot(diff(log10(iter_info_BLP_contraction.DIST_table)))

%% Simulation 2
I=2;
J=10;% Number of products per nest
f_hetero=10;
beta_0=-1;

run run_RCNL_iterations_static.m
results_2=results;

if 1==0
%% Simulation 3
I=100;
J=10;% Number of products per nest
f_hetero=10;
beta_0=20;

run run_RCNL_iterations_static.m
results_3=results;
end% 1==1

end

if 1==1
    %% Simulation 1- SQUAREM
I=2;
J=2;% Number of products per nest

spec_default.SQUAREM_spec=1;

run run_RCNL_iterations_static.m
results_1_SQUAREM=results;

%semilogy(iter_info_V_0_spectral.alpha_table)

comparison_DIST=log10([iter_info_BLP_contraction_spectral0.DIST_table,...
    iter_info_BLP_contraction_spectral.DIST_table]);

    clf
    f = figure;
    set(f,'units','centimeters','position',[10,10,12,9])

plot(comparison_DIST)
xlabel('Iterations')
DIST_label='$\log_{10}\left(\left\Vert \log(S_{j}^{(data)})-\log(s_{j})\right\Vert _{\infty}\right)$';
ylabel(DIST_label,'Interpreter','latex')
legend('Spectral-S3','SQUAREM-S3','Box','off')
uniformFigureStyle(f);

saveas(gcf,append(BLP_paper_figure_path,'spectral_SQUAREM_comparison.png'))

 close all
end
