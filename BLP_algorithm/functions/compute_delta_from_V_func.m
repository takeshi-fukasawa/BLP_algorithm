%%%%%
function delta=compute_delta_from_V_func(...
    mu_ijt,weight,S_jt_data,V)

    global Pr0_spec

    %%% rho==0 case 
    %%% mu_ijt: include continuation values
   %%% Sum of weight is not necesasarily equal to 1

    [J,ns,G,T]=size(mu_ijt);
    n_dim_V=size(V,5);


    n_state=size(weight,5); %% n_dim_V==1 => n_state=1; n_dim_V==n_state>1 under multiple states; IVS => n_dim_V==n_state==1 (only one obs pt)

    s_ijt_ccp_up_to_scale=exp(reshape(mu_ijt,J,ns,G,T,1)-reshape(V,1,ns,1,T,n_dim_V));

   if Pr0_spec==0
    s_jt_up_to_scale=sum(reshape(weight,1,ns,1,1,1).*...
        s_ijt_ccp_up_to_scale,[2,5]);%J*1*G*T


    delta=log(S_jt_data)-log(s_jt_up_to_scale);%J*1*G*T

    else % sequentially compute Pr0 and delta
        Pr0=ones(1,ns,1,T+1,1);
        delta=NaN(J,1,G,T);

       for t=1:T
           s_jt_up_to_scale_t=...
               sum(reshape(weight,1,ns,1,1,n_state).*Pr0(:,:,:,t,:).*...
               s_ijt_ccp_up_to_scale(:,:,:,t,:),[2,5]);%J*1*G*T
    
        exp_delta_t=S_jt_data(:,:,:,t)./s_jt_up_to_scale_t;%J*1*G*T
        delta(:,:,:,t)=log(exp_delta_t);%J*1*G*T
    
         s_ijt_ccp_t=s_ijt_ccp_up_to_scale_t.*exp_delta_t;%J*ns*G*T
         s_i0t_ccp_t=1-sum(s_ijt_ccp_t,[1,3]);%1*ns*1*1*n_state

        Pr0(:,:,:,t+1,:)=Pr0(:,:,:,t,:).*s_i0t_ccp_t;
    end % for loop

    
  end % if statement

end
