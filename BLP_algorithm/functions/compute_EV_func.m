function EV=compute_EV_func(V,IV,weight_V)

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

else
    IV_obs_pt=IV(:,:,:,:,1);%1*ns*G*T
    IV_grid=IV(:,:,:,:,2:end);%1*ns*G*T*n_grid_IV
    
    EV=cat(4,V(:,:,:,2:T,1),V(:,:,:,T,1));
    EV=repmat(EV(:,:,:,:,1),1,1,1,1,n_dim_V);%%%%%

%%%%% Insert 
    %%%%%%%%%%
    %%% Compute coefficients
    coef_1=0.9;
    coef_0=0;
    %%%%%%%%%%%%%%

end

end