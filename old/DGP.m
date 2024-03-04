%% Generate variables (x_j,xi_j,z_j,p_j,nu_i)
c=3;
c=0;
theta=1;
alpha=-2;
sigma_true=0.9;
xi_sd=10;
eps_sd=10;

sigma_est=sigma_true;
if 1==0
    sigma_est=0.2;
    rho_est=0.99;
end

gpurng(100);% set seed of random number generation
x_j=randn(J,1,G);%J by 1; obs. product char.
xi_j=xi_sd*randn(J,1,G);%J by 1; unobs. product char.
z_j_1=0.3*randn(J,1,G);%J by 1; IV
z_j_2=0.3*randn(J,1,G);%J by 1; IV
p_j=1+0.1*x_j+0.5*z_j_1+0.5*z_j_2+0.2*xi_j;%J by 1; price
nu_i=randn(1,I,1);%1 by I; consumer heterogeneity
weight=1/ns*ones(1,ns,1);%1*ns*1



mu_ij_est=sigma_est*x_j.*nu_i./eps_sd;%J by I

delta_j_true=(c+theta*x_j+alpha*p_j+xi_j)./eps_sd;%J by 1
mu_ij_true=sigma_true*x_j.*nu_i./eps_sd;%J by I
u_ij_tilde=delta_j_true+mu_ij_true;%J by I


