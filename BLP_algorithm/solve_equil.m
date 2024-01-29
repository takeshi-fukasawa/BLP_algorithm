ITER_MAX=2000;

if GPU_spec==1
    delta_jt_true=gpuArray(delta_jt_true);
    mu_ijt_true=gpuArray(mu_ijt_true);
    mu_ijt_est=gpuArray(mu_ijt_est);
    weight=gpuArray(weight);
end

V_initial=zeros(1,ns,1,T,n_dim_V);%1*ns*1*T*n_dim_V


for iter=1:ITER_MAX

resid_V=Bellman_update_func(...
    V_initial,delta_jt_true,mu_ijt_true,...
    beta_C,L,rho_est,weight);

V_updated=V_initial-resid_V{1};

DIST=max(abs(V_updated(:)-V_initial(:)));

if beta_C==0 | DIST<TOL
    break;
else
    V_initial=V_updated;
end

end %for loop for solving V

u_ijt_tilde=delta_jt_true+mu_ijt_true+(beta_C^L).*V_updated;%J*I*G*T

weight_V=[];
IV=[];
EV=compute_EV_func(V_updated,IV,weight_V);%1*ns*1*T*n_dim_V
u_i0t_tilde=beta_C*EV;%J*I*1*T*n_dim_V

[s_jt_predict,ChoiceProb_true,s_ijt_given_g_ccp,s_igt_ccp_true,...
    numer_1,denom_1,numer_2,denom_2]=...
    share_func(u_ijt_tilde(:,:,:,:,1),u_i0t_tilde(:,:,:,:,1),rho_true,weight);

s_i0t_ccp_true=1-sum(s_igt_ccp_true,3);%1*ns
s_i0t_true=s_i0t_ccp_true.*weight;%1*ns

s_ijt_true=(ChoiceProb_true.*weight);


S_jt_data=sum(weight.*ChoiceProb_true,2);%J*1*G
S_0t_data=1-sum(S_jt_data,[1,3]);

S_gt_data=sum(S_jt_data,1);%1*1*G
S_jt_given_g_data=S_jt_data./S_gt_data;%J*1*G

numer_1_without_delta=exp(mu_ijt_est./(1-rho_est));%J*ns

V_true=V_updated;
