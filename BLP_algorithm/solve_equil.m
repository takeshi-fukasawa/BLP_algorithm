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
u_ijt_tilde_true=other_vars.u_ijt_tilde;
v_i0t_tilde_true=other_vars.v_i0t_tilde;

[s_jt_predict_true,s_ijt_ccp_true,s_ijt_given_g_ccp_true,s_igt_ccp_true,...
    numer_1,denom_1,numer_2,denom_2,Pr0_true]=...
    share_func(u_ijt_tilde_true(:,:,:,:,1),...
    v_i0t_tilde_true(:,:,:,:,1),rho_true,weight);

S_jt_data=s_jt_predict_true;
S_0t_data=1-sum(S_jt_data,[1,3]);%1*1*1*T

S_gt_data=sum(S_jt_data,1);%1*1*G*T
S_jt_given_g_data=S_jt_data./S_gt_data;%J*1*G*T

s_i0t_ccp_true=1-sum(s_ijt_ccp_true,[1,3]);%1*ns*1*T

IV_state_obs_pt=IV_true(:,:,:,:,1);%1*ns*G*T
IV_state_grid=IV_true(:,:,:,:,2:end);%1*ns*G*T*n_grid_IV
n_grid_IV=size(IV_state_grid,5);

if T>=2
    [coef_0_true,coef_1_true,y,y_predict,R2]=estimate_IV_AR1_transition_func(IV_state_obs_pt);    
end
