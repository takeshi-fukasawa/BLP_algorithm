%% Compute delta, given nonlinear parameters (Update V)
numer_1_without_delta=exp(mu_ij./(1-rho));%J*ns

gpurng(100);% set seed of random number generation

V_initial0=(-log(S_0_data)).*ones(1,ns,'gpuArray');%1*ns
V_hat_initial0_check=log(exp(V_initial0)-1)./(1-rho);


V_initial=V_initial0;

DIST=100;

global denom_1 numer_2 denom_2;
ITER_MAX=1000;

for k=1:ITER_MAX

    V_initial_temp2=V_initial;
    V_hat_initial_temp2=log(exp(V_initial_temp2)-1)./(1-rho);

    output=...
    r_update_RCL_func2(...
    V_initial,weight,mu_ij,rho,...
    S_j_data,S_0_data,numer_1_without_delta);
    
    resid_V_temp2=output{1};

    V_updated=V_initial-resid_V_temp2;

    V_hat_updated=log(exp(V_updated)-1)./(1-rho);

    DIST=max(abs(V_updated(:)-V_initial(:)));%scalar
 
    delta_updated_2=compute_delta_func(...
        numer_1_without_delta,denom_1,numer_2,denom_2,weight,S_j_data,rho);%J*1

    DIST=max(abs(delta_j_true(:)-delta_updated_2(:)));%scalar
    DIST_MAT(k,3)=DIST;
    
    %%%%%
    %k
    %DIST
      if DIST<TOL
          break;% end the for loop
      else
          V_initial_temp=V_initial;
          V_initial=V_updated;%J by 1
      end

end% for loop

denom_1_temp2=denom_1;
numer_2_temp2=numer_2;
denom_2_temp2=denom_2;

V_updated_temp2=V_updated;

s_i0_initial_temp2=s_i0_initial;
s_i0_ccp_initial_temp2=s_i0_ccp_initial;
s_ig_predict_temp2=s_ig_predict;
weight_predict_temp2=weight_predict;
difference_temp2=difference;
V_updated_temp2=V_updated;
V_hat_updated_temp2=V_hat_updated;

n_iteration_update_V=k;
ratio_delta_V=delta_updated_2./delta_j_true;
[min(ratio_delta_V),max(ratio_delta_V)]

log_dist_new=log10(DIST_MAT);

%plot(log_dist_new)
