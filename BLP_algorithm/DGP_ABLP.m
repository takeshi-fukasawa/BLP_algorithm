%% Generate variables (x_j,xi_j,z_j,p_j,nu_i)

gpurng(100);% set seed of random number generation
x_jt=randn(J,1,G,T);% obs. product char.
w_jt=randn(J,1,G,T);

Cov_mat=[1 -0.8 0.3; -0.8 1 0.3; 0.3 0.3 1];
xi_sd=1;

mu=[0 0 0];

x_jt=mvnrnd(mu,Cov_mat,J*G*T);%(J*G*T)*3

xi_jt=1*randn(J,1,G,T);

u_jt=5*rand(J,1,G,T);


beta_x=[1.5;1.5;0.5];
alpha=-3;
sigma_const=0.5;
sigma_x=[0.5 0.5 0.5];
sigma_p=0.2;

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

p_jt=3+1.5*xi_jt+u_jt+reshape(sum(reshape(x_jt,J*G*T,3),2),J,1,G,T);

nu_i=randn(1,I,1);%1 by I; consumer heterogeneity
weight=1/ns*ones(1,ns,1);%1*ns*1


%% Compute product market share, given parameters
delta_jt_true=(beta_0+reshape(x_jt*beta_x,J,1,G,T)+alpha*p_jt+xi_jt)./eps_sd;%J*1*G*T


mu_ijt_true=sigma_const*randn(1,ns,1,1)+...
    sum(reshape(sigma_x,1,1,1,1,3).*reshape(x_jt,J,1,G,T,3).*rand(1,ns,1,1,3),5)+...
    sigma_p*reshape(p_jt,J,1,G,T).*rand(1,ns,1,1);

mu_ijt_true=reshape(mu_ijt_true,J,ns,G,T)./eps_sd;

if 1==0
%%%%% same char over time %%%%%%
delta_jt_true=repmat(delta_jt_true(:,:,:,1),1,1,1,T);
mu_ijt_true=repmat(mu_ijt_true(:,:,:,1),1,1,1,T);
%%%%%%%%%%%
end


if large_hetero_spec==1
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
    weight=[0.5,0.5];
    %weight=[0.1,0.9];
    weight=[0.1,0.9];
end



if mistake_spec==0
    mu_ijt_est=mu_ijt_true+randn(J,ns,G,T)*0.0;%J*ns*G*T
else
    mu_ijt_est=mu_ijt_true*2;%J*ns*G*T
    mu_ijt_est=mu_ijt_true*0.8;%J*ns*G*T
    
end

