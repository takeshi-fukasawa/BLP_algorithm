function EV=compute_EV_func(V,IV_state,weight_V,x_V,T)

global inv_multiply_Chebyshev_basis Chebyshev_extrema

[~,ns,~,n_dim_V]=size(V);
n_grid_IV=n_dim_V-T;
n_dim_Chebyshev=n_grid_IV;

if T==1 %% stationary expec
   EV=V;
elseif T>=2 & n_grid_IV==0 %% Perfect foresight
   %EV=cat(4,V(:,:,:,2:T),zeros(1,ns,1,1));% After the terminal period, no market ...?
   EV=cat(4,V(:,:,:,2:T),V(:,:,:,T)); % After the terminal period, stationary environment
   
else % Inclusive value sufficiency (IVS); Currently, G==1 case only
    IV_state_obs_pt=IV_state(:,:,:,1:T);%1*ns*G*T
    IV_state_grid=IV_state(:,:,:,T+1:end);%1*ns*G*n_grid_IV
    [coef_0,coef_1,y,y_predict,R2]=estimate_IV_AR1_transition_func(IV_state_obs_pt);

    y_mat=[y_predict(:),y(:)];

    sigma=sum((y-y_predict).^2,4)./(T-1);%1*ns*G*1; std of y
    
    n_draw=size(weight_V,1);

    if n_draw==1
        sigma=sigma.*0;%%%%%%% test %%%%%
    end

    IV_t1_state_draw=...
        IV_state.*coef_1+coef_0+...
        sigma.*reshape(x_V,1,1,1,1,n_draw);%1*ns*G*n_dim_V*1*n_draw

    %%% Construct Chebyshev basis
    y=(reshape(V(1,:,1,T+1:end),ns,n_grid_IV))';%n_grid_IV*ns
    coef=inv_multiply_Chebyshev_basis*y;%n_dim_Chebyshev*ns

    IV_max=IV_state_grid(:,:,:,1);%1*ns*1*1
    IV_min=IV_state_grid(:,:,:,end);%1*ns*1*1

    IV_t1_state_draw_scaled=(IV_t1_state_draw-IV_min)./(IV_max-IV_min);%1*ns*1*n_dim_V*n_draw; in [0,1]
    
    IV_t1_state_draw_scaled=2*...
        reshape(IV_t1_state_draw_scaled,...
        ns*1*n_dim_V*n_draw,1)-1;%(ns*1*n_dim_V*n_draw)*1; in [-1,1]

    IV_t1_state_draw_scaled_true=IV_t1_state_draw_scaled;

    %%%%%%%%%%%%%%%%%%
   
    if 1==0
    %%% Use grid points as IV_t1 (for validation)
    IV_t1_state_draw_scaled_temp=zeros(ns,1,1,n_dim_V,n_draw);
    IV_t1_state_draw_scaled_temp(:,:,:,T+1:end,:)=...
        repmat(reshape(Chebyshev_extrema,1,1,1,n_grid_IV,1),[ns,1,1,1,1]);
    IV_t1_state_draw_scaled_temp=IV_t1_state_draw_scaled(:);
    end

    IV_t1_state_draw_scaled=IV_t1_state_draw_scaled_true;%%%%%
   

    %%%%%%%%%%%%%%

    basis_t1=construct_Chebyshev_basis_func(...
        IV_t1_state_draw_scaled,n_dim_Chebyshev);%(ns*n_dim_V*n_draw)*n_dim_Chebyshev

    V_t1_draw=sum(reshape(basis_t1,1,ns,1,n_dim_V,n_draw,n_dim_Chebyshev).*...
        reshape(coef',1,ns,1,1,1,n_dim_Chebyshev),6);%1*ns*G*n_dim_V*n_draw

    %%%%%%
    %V_t1_draw=repmat(V,[1,1,1,1,n_draw]);%1*ns*G*n_dim_V*n_draw
    %%%%%%%

    EV=sum(V_t1_draw.*reshape(weight_V,1,1,1,1,n_draw),5);%1*ns*G*n_dim_V
    
end

end