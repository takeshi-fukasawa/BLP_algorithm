%% Generate variables (x_j,xi_j,z_j,p_j,nu_i)
n_state=1;

gpurng(100);% set seed of random number generation
w_jt=randn(J,1,G,T);

Cov_mat=[1 -0.8 0.3; -0.8 1 0.3; 0.3 0.3 1];

if durable_spec==1
Cov_mat=0.5*eye(3);
end

xi_sd=1;


x_jt=mvnrnd([0,0,0],Cov_mat,J*G*T);%(J*G*T)*3; obs. product char.

if durable_spec==1
    x_jt=repmat(mvnrnd([0,0,0],Cov_mat,J*G),T,1);%(J*G*T)*3
end

xi_jt=1*randn(J,1,G,T);

u_jt=5*rand(J,1,G,T);


beta_x=[1.5;1.5;0.5];
alpha=-3;
sigma_const=0.5;
sigma_x=[0.5 0.5 0.5];
sigma_p=0.2;

if durable_spec==1
beta_x=[1;1;0.5];
alpha=-2;
sigma_const=0;
sigma_x=[0.5 0.5 0];
sigma_p=0.25;

end

if 1==0
    sigma_const=sigma_const*0.0;
    sigma_x=sigma_x*0.0;
    sigma_p=sigma_p*0.0;
end

%sigma_const=sigma_const*5;
%sigma_x=sigma_x*5;
%sigma_p=sigma_p*5;

eps_sd=1;

if eps_sd~=1
    beta_0=6;
end

p_jt=3+1.5*xi_jt+u_jt+...
    reshape(sum(reshape(x_jt,J*G*T,3),2),...
    J,1,G,T);%J*1*G*T

if durable_spec==1
rho_0=0.1;
rho_z=0.95;
z_j0=8*ones(J,1,G,1);
z_jt=zeros(J,1,G,T+1);
z_jt(:,:,:,1)=z_j0;

for t=1:T
eta_jt=0.1*randn(J,1,G,1);
z_jt(:,:,:,t+1)=rho_0+rho_z*z_jt(:,:,:,t)+eta_jt;
end
z_jt=z_jt(:,:,:,2:end);%J*1*G*T

gamma_0=1;
gamma_x=[0.2;0.2;0.1];
gamma_z=1;
gamma_w=0.2;
gamma_xi=0.7;
gamma_p=[0.1;0.1;0.1];

u_jt=0.01*randn(J,1,G,T);
w_jt=randn(J,1,G,T);

%%% Competitors' char ??? %%%%%%
p_jt=gamma_0+reshape(x_jt*(gamma_x+gamma_p),J,1,G,T)+gamma_z*z_jt+...
    gamma_w*w_jt+gamma_xi*xi_jt...
    -sum(reshape(x_jt*gamma_p,J,1,G,T),1)+... %%% reshape(x_jt*gamma_p,J,1,G,T),1)+...
    u_jt;%J*1*G*T
%p_jt=repmat(p_jt(:,:,:,1),[1,1,1,T]);
end


p_mean=reshape(mean(p_jt,[1,3]),T,1);

nu_i=randn(1,I,1);%1 by I; consumer heterogeneity
weight=1/ns*ones(1,ns,1);%1*ns*1


%% Compute product market share, given parameters
delta_jt_true=(beta_0+reshape(x_jt*beta_x,J,1,G,T)+alpha*p_jt+xi_jt)./eps_sd;%J*1*G*T


mu_ijt_true=sigma_const*randn(1,ns,1,1)+...
    sum(reshape(sigma_x,1,1,1,1,3).*reshape(x_jt,J,1,G,T,3).*randn(1,ns,1,1,3),5)+...
    sigma_p*reshape(p_jt,J,1,G,T).*randn(1,ns,1,1);

mu_ijt_true=reshape(mu_ijt_true,J,ns,G,T)./eps_sd;
mu_i0t_true=zeros(1,ns,1,T,n_state)./eps_sd;
mu_i0t_est=mu_i0t_true;

IV_true=log(sum(exp(delta_jt_true+mu_ijt_true)./eps_sd,1));%1*ns*G*T




if 1==0
%%%%% same char over time %%%%%%
delta_jt_true=repmat(delta_jt_true(:,:,:,1),1,1,1,T);
mu_ijt_true=repmat(mu_ijt_true(:,:,:,1),1,1,1,T);
%%%%%%%%%%%
end


if large_hetero_spec==1
   %%% Test large consumer heterogeneity case %%%%
    mu_ijt_true=(mu_base)*ones(size(mu_ijt_true));
    %%%mu_ijt_est(1,2)=2;
    if J==1
        %%% Introduce heterogeneity
        mu_ijt_true(1,1)=0;%-f_hetero;
        mu_ijt_true(1,2)=f_hetero;
    elseif J==2
        mu_ijt_true(1,1)=f_hetero;
        mu_ijt_true(2,2)=f_hetero;
    elseif J==10
        mu_ijt_true=f_hetero*randn(J,ns);
    end

    delta_jt_true=beta_0*ones(J,1);
    %delta_jt_true(1)=delta_jt_true(1)+0.5;
    
    delta_jt_true(1)=0;
    delta_jt_true(2)=-1;
    weight=[0.5,0.5];
    %weight=[0.1,0.9];
    weight=[0.1,0.9];
end

if large_hetero_spec==2
   %%% Test ill-condition problem case %%%%
    delta_jt_true=[-10;10];
    weight=[0.5,0.5];
    J=size(delta_jt_true,1);

    mu_ijt_true=[0,1;0,0];
    mu_ijt_true=1*rand(J,ns,1,1,1);
    mu_ijt_true(1,2)=0;

    mu_ijt_true=zeros(J,ns,1,1,1);
    mu_ijt_true(1,1)=5;
    mu_ijt_true(2,2)=0;
    
    
end


if mistake_spec==0
    mu_ijt_est=mu_ijt_true+randn(J,ns,G,T)*0.0;%J*ns*G*T
else
    mu_ijt_est=mu_ijt_true*2;%J*ns*G*T
    mu_ijt_est=mu_ijt_true*0.8;%J*ns*G*T
    
end

