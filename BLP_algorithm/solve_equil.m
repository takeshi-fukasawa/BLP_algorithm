ITER_MAX=3000;

[x_V,weight_V]=gausshermi(n_draw);
x_V=x_V*sqrt(2);%n_draw*1
weight_V=weight_V./sum(weight_V,1);%n_draw*1

V_initial0=zeros(1,ns,1,T,n_dim_V);%1*ns*1*T*n_dim_V

mu_ijt=mu_ijt_true;
run run_Bellman.m

V_true=V_updated;
V_data_true=V_updated(:,:,:,:,1);

IV_true=other_vars.IV;

s_igt_ccp_true=exp(IV_true(:,:,:,:,1)-V_true(:,:,:,:,1));

s_ijt_given_g_ccp_true=...
    other_vars.numer_1(:,:,:,:,1)./other_vars.denom_1(:,:,:,:,1);
s_ijt_ccp_true=s_ijt_given_g_ccp_true.*s_igt_ccp_true;

s_i0t_ccp_true=1-sum(s_igt_ccp_true,3);

if 1==1
    S_jt_data=sum(s_ijt_ccp_true.*weight,2);
    
else
Pr0=ones(1,ns,1,T);
S_jt_data=NaN(J,1,G,T);
s_ijt_ccp_true=NaN(J,ns,G,T);
s_igt_ccp_true=NaN(1,ns,G,T);
for t=1:T
     weight_Pr0_t=reshape(weight,1,ns,1,1).*Pr0(:,:,:,t);%1*ns*1*1
     [S_jt_data_temp,s_ijt_ccp_temp,...
      ~,s_igt_ccp_temp]=...
    share_func(u_ijt_tilde(:,:,:,1,1),u_i0t_tilde(:,:,:,1,1),rho_true,weight_Pr0_t);
     S_jt_data(:,:,:,t)=S_jt_data_temp;
    s_ijt_ccp_true(:,:,:,t)=s_ijt_ccp_temp;
    s_igt_ccp_true(:,:,:,t)=s_igt_ccp_temp;

    if t<=T-1
    Pr0(:,:,:,t+1)=Pr0(:,:,:,t).*(1-sum(s_ijt_ccp_temp,[1,3]));%1*ns*1*1
    end
end % for loop
end % if else statement


S_0t_data=1-sum(S_jt_data,[1,3]);%1*1*1*T

S_gt_data=sum(S_jt_data,1);%1*1*G*T
S_jt_given_g_data=S_jt_data./S_gt_data;%J*1*G*T

IV_state_obs_pt=IV_true(:,:,:,:,1);%1*ns*G*T
IV_state_grid=IV_true(:,:,:,:,2:end);%1*ns*G*T*n_grid_IV
n_grid_IV=size(IV_state_grid,5);

if T>=2
    %%%% Estimate AR1 coefficients of IV state transitions %%%%%%
    X=IV_state_obs_pt(:,:,:,1:end-1);%1*ns*G*(T-1)
    y=IV_state_obs_pt(:,:,:,2:end);%1*ns*G*(T-1)
    X_mean=mean(X,4);%1*ns*G*1
    y_mean=mean(y,4);%1*ns*G*1
    Var_X=sum((X-X_mean).^2,4);%1*ns*G*1; not divided by the number of elements; not normalized
    Cov_X_y=sum((X-X_mean).*(y-y_mean),4);%1*ns*G*1; not divided by the number of elements; not normalized
    coef_1_true=Cov_X_y./Var_X;%1*ns*G*1
    coef_0_true=y_mean-coef_1_true.*X_mean;%1*ns*G*1
    
    y_predict=X.*coef_1_true+coef_0_true;%1*ns*G*(T-1)
    R2=corr(y_predict(:),y(:));
end
