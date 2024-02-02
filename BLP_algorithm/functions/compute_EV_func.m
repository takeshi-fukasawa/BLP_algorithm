function EV=compute_EV_func(V,IV_state,weight_V,x_V)

[~,ns,~,T,n_dim_V]=size(V);

EV=V;

if T==1 & n_dim_V==1 %% stationary expec
   EV=V;
elseif T>=2 & n_dim_V==1 %% Perfect foresight
   %EV=cat(4,V(:,:,:,2:T),zeros(1,ns,1,1));
   EV=cat(4,V(:,:,:,2:T),V(:,:,:,T));
   %EV=cat(4,V(:,:,:,2:T),V(:,:,:,1));
   %EV=cat(4,V(:,:,:,1),V(:,:,:,3:T),V(:,:,:,1));
   
   %EV=repmat(V(:,:,:,1),1,1,1,T);
   
   EV=repmat(V(:,:,:,:,1),1,1,1,1,n_dim_V);%%%%%
   %EV=V*0.9+0.5;
   %EV=V*0.99+mean(V,4)/100;%%%%
   %EV=mean(V,4);

else % Inclusive value sufficiency (IVS); Currently, G==1 case only
    IV_state_obs_pt=IV_state(:,:,:,:,1);%1*ns*G*T
    IV_state_grid=IV_state(:,:,:,:,2:end);%1*ns*G*T*n_grid_IV
    n_grid_IV=size(IV_state_grid,5);

    %%%% Compute Ar1 coefficients %%%%%%
    X=IV_state_obs_pt(:,:,:,1:end-1);%1*ns*G*(T-1)
    y=IV_state_obs_pt(:,:,:,2:end);%1*ns*G*(T-1)
    X_mean=mean(X,4);%1*ns*G*1
    y_mean=mean(y,4);%1*ns*G*1
    Var_X=sum((X-X_mean).^2,4);%1*ns*G*1
    Cov_X_y=sum((X-X_mean).*(y-y_mean),4);%1*ns*G*1
    coef_1=Cov_X_y./Var_X;%1*ns*G*1
    coef_0=y_mean-coef_1.*X_mean;%1*ns*G*1
    
    y_predict=X.*coef_1+coef_0;%1*ns*G*(T-1)
    R2=corr(y_predict(:),y(:))
    
    n_draw=size(weight_V,1);

    x_V=x_V*0;%%%%%test %%%%%

    IV_t1_state_draw=...
        IV_state.*coef_1+coef_0+...
        reshape(x_V,1,1,1,1,1,n_draw);%1*ns*G*T*n_grid_IV*n_draw


    %%% Construct Chebyshev basis
    n_dim_Chebyshev=n_grid_IV;

    Chebyshev_extrema=cos([0:n_dim_Chebyshev-1]*pi/(n_dim_Chebyshev-1));% fixed in the iteration

    basis_t_grid=construct_Chebyshev_basis_func(...
        reshape(Chebyshev_extrema,n_grid_IV,1),...
        n_dim_Chebyshev);%n_grid_IV*n_dim_Chebyshev; fixed in the iteration

    
    y=(reshape(V(1,:,1,1,2:end),ns,n_grid_IV))';%n_grid_IV*ns

    coef=inv(basis_t_grid'*basis_t_grid)*basis_t_grid'*y;%n_dim_Chebyshev*ns

    IV_max=IV_state_grid(1,1,1,1,2);
    IV_min=IV_state_grid(1,1,1,1,end);
    IV_t1_state_draw_scaled=reshape(IV_t1_state_draw-IV_min,ns*1*T*n_dim_V*n_draw,1)/(IV_max-IV_min);% in [0,1]
    IV_t1_state_draw_scaled=2*IV_t1_state_draw_scaled-1;%(ns*1*T*n_dim_V*n_draw)*1; in [-1,1]

    basis_t1=construct_Chebyshev_basis_func(...
        IV_t1_state_draw_scaled,n_dim_Chebyshev);%(ns*T*n_grid_IV*n_draw)*n_dim_Chebyshev

    V_t1_draw=sum(reshape(basis_t1,1,ns,1,T,n_dim_V,n_draw,n_dim_Chebyshev).*...
        reshape(coef',1,ns,1,1,1,1,n_dim_Chebyshev),7);%1*ns*G*T*n_grid_IV*n_draw

    EV=sum(V_t1_draw.*reshape(weight_V,1,1,1,1,1,n_draw),6);%1*ns*G*(T-1)*n_grid
    
    %EV=cat(4,V(:,:,:,2:T,1),V(:,:,:,T,1));
    %EV=repmat(EV(:,:,:,:,1),1,1,1,1,n_dim_V);%%%%%

end

end