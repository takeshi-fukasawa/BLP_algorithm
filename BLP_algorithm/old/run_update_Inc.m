

%% Compute delta, given nonlinear parameters (Update Inclusive value)
numer_1_without_delta=exp(mu_ij_est./(1-rho_est));%J*ns
DIST_MAT_Inc=zeros(ITER_MAX,1);

global posterior_i_j posterior_i_0

gpurng(100);% set seed of random number generation

Inc_initial0=log(sum(S_j_data,1))-log(S_0_data)+0.0*rand(1,ns);

Inc_initial=Inc_initial0;

DIST=100;

global denom_1 numer_2 denom_2;


for k=1:ITER_MAX

    output=...
    update_Inc_func(...
    Inc_initial,weight,mu_ij_est,rho_est,...
    S_j_data,S_0_data,numer_1_without_delta);
    resid=output{1};

    Inc_updated=Inc_initial-resid;

    DIST=max(abs(Inc_updated(:)-Inc_initial(:)));%scalar

    delta_updated_Inc=compute_delta_func(numer_1_without_delta,denom_1,numer_2,denom_2,weight,S_j_data,rho_est);%J*1

    %DIST=max(abs(delta_j_true(:)-delta_updated_Inc(:)));%scalar
    DIST_MAT_Inc(k,2)=DIST;
    
    %%%%%
    %k
    %DIST
      if DIST<TOL
          break;% end the for loop
      else
          Inc_initial=Inc_updated;%J by 1
      end

end% for loop


denom_1_temp1=denom_1;
numer_2_temp1=numer_2;
denom_2_temp1=denom_2;

n_iter_update_Inc=k;
ratio_delta_Inc=delta_updated_Inc./delta_j_true;
[min(ratio_delta_Inc),max(ratio_delta_Inc)]

log_dist_Inc=log10(DIST_MAT_Inc);

%plot(log_dist_new)

%% Update Inc val spectral
output_spectral=...
    spectral_func(@update_Inc_func,1,Inc_initial0,...
    [],[],...
    weight,mu_ij_est,rho_est,...
    S_j_data,S_0_data,numer_1_without_delta);
    
    Inc_sol=output_spectral{1};


delta_sol_Inc_spectral=compute_delta_func(numer_1_without_delta,denom_1,numer_2,denom_2,weight,S_j_data,rho_est);%J*1

n_iter_update_Inc_spectral=count
DIST_MAT_Inc_spectral=DIST_table;
ratio_delta_Inc_spectral=delta_sol_Inc_spectral./delta_j_true;
[min(ratio_delta_Inc_spectral),max(ratio_delta_Inc_spectral)]

log_dist_Inc_spectral=log10(DIST_MAT_Inc_spectral);

