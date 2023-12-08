function EV=compute_EV_func(V,weight_V)

[~,ns,~,T,n_grid]=size(V);

EV=V;

if T==1 & n_grid==1 %% stationary expec
   EV=V;
elseif T>=2 & n_grid==1 %% Perfect foresight
   %EV=cat(4,V(:,:,:,2:T),zeros(1,ns,1,1));
   %EV=cat(4,V(:,:,:,2:T),V(:,:,:,T));
   %EV=repmat(V(:,:,:,1),1,1,1,T);
   
else
%%%%%%%%
end

end