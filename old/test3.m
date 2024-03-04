%% Compute delta, given nonlinear parameters (Update V spec2)
numer_1_without_delta=exp(mu_ij./(1-rho));%J*ns

gpurng(100);% set seed of random number generation

V_initial0=(-log(S_0_data)).*ones(1,ns,'gpuArray');%1*ns
V_hat_initial0=log(exp(V_initial0)-1)./(1-rho);
V_initial0_test=log(1+exp((1-rho)*V_hat_initial0));

ITER_MAX=800;

V_hat_initial=V_hat_initial0;

DIST=100;

global denom_1 numer_2 denom_2;


for k=1:ITER_MAX

    V_hat_initial_temp3=V_hat_initial;
    V_initial_temp3=log(1+exp((1-rho)*V_hat_initial_temp3));

    output=...
    r_update_RCL_func3(...
    V_hat_initial,weight,mu_ij,rho,...
    S_j_data,S_0_data,numer_1_without_delta);
    resid_V_hat_temp3=output{1};

    V_hat_updated=V_hat_initial-resid_V_hat_temp3;
    V_updated=log(1+exp((1-rho)*V_hat_updated));

    DIST=max(abs(V_hat_updated(:)-V_hat_initial(:)));%scalar
 
    delta_updated_3=compute_delta_func(...
        numer_1_without_delta,denom_1,numer_2,denom_2,weight,S_j_data,rho);%J*1

    DIST=max(abs(delta_j_true(:)-delta_updated_3(:)));%scalar
    DIST_MAT(k,4)=DIST;
    
    %%%%%
    %k
    %DIST
      if DIST<TOL
          break;% end the for loop
      else

          V_hat_initial_temp=V_hat_initial;
          V_hat_initial=V_hat_updated;%J by 1
      end

end% for loop

denom_1_temp3=denom_1;
numer_2_temp3=numer_2;
denom_2_temp3=denom_2;


s_i0_initial_temp3=s_i0_initial;
s_i0_ccp_initial_temp3=s_i0_ccp_initial;
s_ig_predict_temp3=s_ig_predict;
weight_predict_temp3=weight_predict;
difference_temp3=difference;
V_hat_updated_temp3=V_hat_updated;

V_updated_temp3=V_updated;

n_iter_update_V_2=k;
ratio_delta_V_hat=delta_updated_3./delta_j_true;
[min(ratio_delta_V_hat),max(ratio_delta_V_hat)]

log_dist_new=log10(DIST_MAT);

%plot(log_dist_new)
