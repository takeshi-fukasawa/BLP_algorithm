function EV=compute_EV_func(V,IV_state,weight_V,x_V)

[~,ns,~,T,n_dim_V]=size(V);


if T==1 & n_dim_V==1 %% stationary expec
   EV=V;
elseif T>=2 & n_dim_V==1 %% Perfect foresight
   %EV=cat(4,V(:,:,:,2:T),zeros(1,ns,1,1));% After the terminal period, no market ...?
   EV=cat(4,V(:,:,:,2:T),V(:,:,:,T)); % After the terminal period, stationary environment
   
else % Inclusive value sufficiency (IVS); Currently, G==1 case only
    IV_state_obs_pt=IV_state(:,:,:,:,1);%1*ns*G*T
    IV_state_grid=IV_state(:,:,:,:,2:end);%1*ns*G*T*n_grid_IV
    n_grid_IV=size(IV_state_grid,5);

    %%%% Estimate AR1 coefficients of IV state transitions %%%%%%
    X=IV_state_obs_pt(:,:,:,1:end-1);%1*ns*G*(T-1)
    y=IV_state_obs_pt(:,:,:,2:end);%1*ns*G*(T-1)
    X_mean=mean(X,4);%1*ns*G*1
    y_mean=mean(y,4);%1*ns*G*1
    Var_X=sum((X-X_mean).^2,4);%1*ns*G*1; not divided by the number of elements; not normalized
    Cov_X_y=sum((X-X_mean).*(y-y_mean),4);%1*ns*G*1; not divided by the number of elements; not normalized
    coef_1=Cov_X_y./Var_X;%1*ns*G*1
    coef_0=y_mean-coef_1.*X_mean;%1*ns*G*1
    
    y_predict=X.*coef_1+coef_0;%1*ns*G*(T-1)
    R2=corr(y_predict(:),y(:));

    sigma=sum((y-y_predict).^2,4)./(T-1);%1*ns*G*1; std of y
    
    n_draw=size(weight_V,1);

    sigma=sigma.*0;%%%%%%% test %%%%%

    IV_t1_state_draw=...
        IV_state.*coef_1+coef_0+...
        sigma.*reshape(x_V,1,1,1,1,1,n_draw);%1*ns*G*T*n_grid_IV*n_draw

    %[min(IV_t1_state_draw(:)),max(IV_t1_state_draw(:))]

    %[min(IV_state_grid(:)),max(IV_state_grid(:))]


    %%% Construct Chebyshev basis
    n_dim_Chebyshev=n_grid_IV;

    Chebyshev_extrema=cos([0:n_dim_Chebyshev-1]*pi/(n_dim_Chebyshev-1));% fixed in the iteration

    basis_t_grid=construct_Chebyshev_basis_func(...
        reshape(Chebyshev_extrema,n_grid_IV,1),...
        n_dim_Chebyshev);%n_grid_IV*n_dim_Chebyshev; fixed in the iteration

    
    y=(reshape(V(1,:,1,1,2:end),ns,n_grid_IV))';%n_grid_IV*ns

    coef=inv(basis_t_grid'*basis_t_grid)*basis_t_grid'*y;%n_dim_Chebyshev*ns

    IV_max=IV_state_grid(:,:,:,:,2);%1*ns*1*1
    IV_min=IV_state_grid(:,:,:,:,end);%1*ns*1*1

    IV_t1_state_draw_scaled=(IV_t1_state_draw-IV_min)./(IV_max-IV_min);%1*ns*1*T*n_draw; in [0,1]
    
    IV_t1_state_draw_scaled=2*...
        reshape(IV_t1_state_draw_scaled,...
        ns*1*T*n_dim_V*n_draw,1)-1;%(ns*1*T*n_dim_V*n_draw)*1; in [-1,1]

    IV_t1_state_draw_scaled_true=IV_t1_state_draw_scaled;

    %%%%%%%%%%%%%%%%%%
    %%% Use grid points as IV_t1 (for validation)
   
    IV_t1_state_draw_scaled_temp=zeros(ns,1,T,n_dim_V,n_draw);
    IV_t1_state_draw_scaled_temp(:,:,:,2:end,:)=...
        repmat(reshape(Chebyshev_extrema,1,1,1,n_grid_IV,1),[ns,1,T,1,1,1]);
    IV_t1_state_draw_scaled_temp=IV_t1_state_draw_scaled(:);
   
    IV_t1_state_draw_scaled=IV_t1_state_draw_scaled_true;%%%%%
   
    %[max(IV_t1_state_draw_scale_temp(:)),...
    %   min(IV_t1_state_draw_scaled_temp(:))]
    %IV_t1_state_draw_scaled=IV_t1_state_draw_scaled_true;%%%%%

    %[max(IV_t1_state_draw_scaled_true(:)),...
    %    median(IV_t1_state_draw_scaled_true(:)),...
    %    min(IV_t1_state_draw_scaled_true(:))]

    %%%%%%%%%%%%%%

    basis_t1=construct_Chebyshev_basis_func(...
        IV_t1_state_draw_scaled,n_dim_Chebyshev);%(ns*T*n_dim_V*n_draw)*n_dim_Chebyshev

    V_t1_draw=sum(reshape(basis_t1,1,ns,1,T,n_dim_V,n_draw,n_dim_Chebyshev).*...
        reshape(coef',1,ns,1,1,1,1,n_dim_Chebyshev),7);%1*ns*G*T*n_dim_V*n_draw

    %%%%%%
    %V_t1_draw=repmat(V,[1,1,1,1,1,n_draw]);%1*ns*G*T*n_dim_V*n_draw
    %%%%%%%

    EV=sum(V_t1_draw.*reshape(weight_V,1,1,1,1,1,n_draw),6);%1*ns*G*(T-1)*n_dim_V
    
    %%%EV=cat(4,V(:,:,:,2:T,1),V(:,:,:,T,1));% Stationary after the terminal period

end

end