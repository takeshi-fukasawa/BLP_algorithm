global IV_temp0

ITER_MAX=3000;

[x_V,weight_V]=gausshermi(n_draw);
x_V=x_V*sqrt(2);%n_draw*1
weight_V=weight_V./sum(weight_V,1);%n_draw*1

V_initial=zeros(1,ns,1,T,n_dim_V);%1*ns*1*T*n_dim_V

%%%%%%%%%
if 1==1
DIST_MAT=NaN(ITER_MAX,1);

for iter=1:ITER_MAX

resid_V=Bellman_update_func(...
    V_initial,delta_jt_true,mu_ijt_true,...
    beta_C,L,rho_est,weight,weight_V,x_V);

V_updated=V_initial-resid_V{1};

DIST=max(abs(V_updated(:)-V_initial(:)));
DIST_MAT(iter,1)=DIST;

if beta_C==0 | DIST<TOL | isnan(DIST)
    break;
else
    V_initial=V_updated;
end

end %for loop for solving V

if iter==ITER_MAX
    warning("Not converge")
end

end
%%%%%%%%
V_initial=zeros(1,ns,1,T,n_dim_V);
%%V_initial=V_updated;
[output_spectral,other_vars,DIST_table_Bellman,fun_k_cell]=...
        spectral_func(@Bellman_update_func,1,t_dim_id,[],V_initial,...
        delta_jt_true,mu_ijt_true,...
    beta_C,L,rho_est,weight,weight_V,x_V);
V_updated2=output_spectral{1};

u_ijt_tilde=delta_jt_true+mu_ijt_true+(beta_C^L).*V_updated;%J*I*G*T

EV=compute_EV_func(V_updated,IV_temp0,weight_V,x_V);%1*ns*1*T*n_dim_V
u_i0t_tilde=beta_C*EV;%J*I*1*T*n_dim_V


V_true=V_updated;
V_data_true=V_updated(:,:,:,:,1);

if 1==1
[S_jt_data,s_ijt_ccp_true,s_ijt_given_g_ccp,s_igt_ccp_true,...
    numer_1,denom_1,numer_2,denom_2]=...
    share_func(u_ijt_tilde(:,:,:,:,1),u_i0t_tilde(:,:,:,:,1),rho_true,weight);
IV_true=log(numer_2);

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


s_i0t_ccp_true=1-sum(s_igt_ccp_true,3);%1*ns*1*T

S_0t_data=1-sum(S_jt_data,[1,3]);%1*1*1*T

S_gt_data=sum(S_jt_data,1);%1*1*G*T
S_jt_given_g_data=S_jt_data./S_gt_data;%J*1*G*T

