clear
%% Parameter settings
c=-3;
c=3;
theta=1;
alpha=-2;
sigma=0.3;
J=100;
I=10;
ns=I;

%% Generate variables (x_j,xi_j,z_j,p_j,nu_i)
rng(100);% set seed of random number generation
x_j=randn(J,1);%J by 1; obs. product char.
xi_j=0.1*randn(J,1);%J by 1; unobs. product char.
z_j_1=0.3*randn(J,1);%J by 1; IV
z_j_2=0.3*randn(J,1);%J by 1; IV
p_j=1+0.1*x_j+0.5*z_j_1+0.5*z_j_2+0.2*xi_j;%J by 1; price
nu_i=randn(1,I);%1 by I; consumer heterogeneity

%% Compute product market share, given parameters
delta_j_true=c+theta*x_j+alpha*p_j+xi_j;%J by 1
mu_ij=sigma*x_j.*nu_i;%J by I
u_ij_tilde=delta_j_true+mu_ij;%J by I
numer=exp(u_ij_tilde);%J by I
denom=1+sum(numer,1);%1 by I
ChoiceProb_true=numer./denom;%J by I
s_j=mean(ChoiceProb_true,2);%J by 1
s_0=1-sum(s_j,1);%1 by 1

S_j_data=s_j;% J by 1
S_0_data=s_0;% J by 1

%% Compute delta, given nonlinear parameters (sigma)
delta_initial=log(S_j_data)-log(S_0_data);%J by 1; Initial value of delta
DIST=100;
TOL=1e-10; % tolerance level to check convergence
ITER_MAX=1000;% maximum number of iterations

u_i0_tilde=zeros(1,ns);

for k=1:ITER_MAX
  u_ij_tilde=delta_initial+sigma*nu_i.*x_j;%J by I
  [s_j_predict,s_ij]=share_func(u_ij_tilde,u_i0_tilde);%J by 1
  %s_i0=1-sum(s_ij,1);% 1 by ns

  %s_j_updated_temp=(mean(exp(u_ij_tilde).*s_i0,2));%J*1
  
  % contraction mapping:
  delta_updated=delta_initial+log(S_j_data)-log(s_j_predict);%J by 1
  DIST=max(abs(delta_updated(:)-delta_initial(:)));%scalar
  
  %%%%%
  DIST=max(abs(delta_j_true(:)-delta_updated(:)));%scalar
  %%%%%

  k
  DIST
  if DIST<TOL
      break;% end the for loop
  else
      delta_initial=delta_updated;%J by 1
  end

end% for loop

n_iteration_BLP=k;

%% Compute delta, given nonlinear parameters (sigma; new)
delta_initial=log(S_j_data)-log(S_0_data);%J by 1; Initial value of delta
s_i0_initial=repmat(S_0_data,1,ns);%1*ns

DIST=100;
TOL=1e-10; % tolerance level to check convergence
ITER_MAX=1000;% maximum number of iterations

for k=1:ITER_MAX
  delta_updated=log(S_j_data)-log(mean(s_i0_initial.*exp(mu_ij),2));%J*1
  u_ij_tilde=delta_updated+sigma*nu_i.*x_j;%J by I
  [s_j_predict,s_ij]=share_func(u_ij_tilde,u_i0_tilde);%J by 1
  s_i0_updated=1-sum(s_ij,1);% 1 by ns

  DIST=max(abs(s_i0_updated(:)-s_i0_initial(:)));%scalar
 
  %%%%%
  DIST=max(abs(delta_j_true(:)-delta_updated(:)));%scalar
  %%%%%
  k
  DIST
  if DIST<TOL
      break;% end the for loop
  else
      s_i0_initial=s_i0_updated;%J by 1
  end

end% for loop

n_iteration_new=k;

%% Compute delta, given nonlinear parameters (Kalouptsidi method)
weight=1/ns*ones(1,ns);%1*ns
r_initial=log(S_0_data.*weight)+10*rand(1,ns);%1 by ns; Initial value
r_initial=r_normalize_func(r_initial,S_0_data);


DIST=100;
TOL=1e-10; % tolerance level to check convergence
ITER_MAX=100;% maximum number of iterations

for k=1:ITER_MAX
    numer=exp(mu_ij+r_initial);%J*ns
    denom=sum(numer,2);%J*1
    frac_j=numer./denom;%J*ns

    numer=exp(r_initial);%1*ns
    denom=sum(numer,2);%1*1
    frac_0=numer./denom;%1*ns
    temp=sum(S_j_data.*frac_j,1)+S_0_data.*frac_0;%1*ns
    
    r_updated=r_initial+log(weight)-log(temp);%J*1
    r_updated=r_normalize_func(r_updated,S_0_data);

    q_updated=exp(r_updated);
    delta_updated=log(S_j_data)-log(sum(q_updated.*exp(mu_ij),2));

    DIST=max(abs(r_updated(:)-r_initial(:)));%scalar
 
    %%%%%
    DIST=max(abs(delta_j_true(:)-delta_updated(:)));%scalar
    %%%%%
    k
    DIST
  if DIST<TOL
      break;% end the for loop
  else
      r_initial=r_updated;%J by 1
  end

end% for loop

n_iteration_new2=k;


% Function to compute market shares s_j given u_ijt_tilde
function [s_j_predict,ChoiceProb]=share_func(u_ij_tilde,u_i0_tilde)
    numer=exp(u_ij_tilde);%J by I
    denom=1+sum(numer,1);%1 by I
    ChoiceProb=numer./denom;%J by I
    s_j_predict=mean(ChoiceProb,2);%J by 1
end

function r_out=r_normalize_func(r,S_0_data)
    r_out=r;
    temp=S_0_data-sum(exp(r(:,1:end-1)),2);%1*1
    if temp>0
        r_out(:,end)=log(temp);
    else
        q_i=exp(r);%1*ns
        q_i=S_0_data*q_i./sum(q_i);%1*ns
        r_out=log(q_i);%1*ns
    end

    %r_out(:,end)=log(temp);

    %S_0_data-sum(exp(r_out),2) %zero??
end