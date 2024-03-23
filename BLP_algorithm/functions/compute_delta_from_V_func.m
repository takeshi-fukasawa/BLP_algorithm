function [delta,Pr0,s_ijt_ccp_up_to_scale]=compute_delta_from_V_func(...
    mu_ijt,weight,S_jt_data,V,Pr0)

    global Pr0_spec

    %%% Used in V_update_func etc.
    %%% rho==0 case 
    %%% mu_ijt: include continuation values
    %%% Sum of weight is not necesasarily equal to 1

    [J,ns,G,T]=size(mu_ijt);

    s_ijt_ccp_up_to_scale=exp(reshape(mu_ijt,J,ns,G,T,1)-...
        reshape(V,1,ns,1,T));%J*ns*G*T


   if Pr0_spec==0
        s_jt_up_to_scale=sum(reshape(weight,1,ns,1,1,1).*...
            s_ijt_ccp_up_to_scale,[2,5]);%J*1*G*T
    
        delta=log(S_jt_data)-log(s_jt_up_to_scale);%J*1*G*T
        Pr0=ones(1,ns,1,T);
        
   elseif Pr0_spec==1 & isempty(Pr0)==1 %%% sequentially compute Pr0 and delta
        delta=NaN(J,1,G,T);
        Pr0=ones(1,ns,1,T);

       for t=1:T
           Pr0_t=Pr0(:,:,:,t,:);
           s_jt_up_to_scale_t=...
               sum(reshape(weight,1,ns,1,1).*Pr0_t.*...
               s_ijt_ccp_up_to_scale(:,:,:,t),2);%J*1*G*1
    
          exp_delta_t=S_jt_data(:,:,:,t)./s_jt_up_to_scale_t;%J*1*G*1
          delta(:,:,:,t)=log(exp_delta_t);%J*1*G*1
      
          s_ijt_ccp_t=s_ijt_ccp_up_to_scale(:,:,:,t,:).*...
              exp_delta_t;%J*ns*G*1
          s_i0t_ccp_t=1-sum(s_ijt_ccp_t,[1,3]);%1*ns*1*1

          if t<=T-1
             Pr0(:,:,:,t+1)=Pr0_t.*s_i0t_ccp_t;%1*ns*1*1
          end
       end % for loop

   elseif Pr0_spec==1 & isempty(Pr0)==0 %%% Recover delta by matrix opetation

           s_jt_up_to_scale=...
               sum(reshape(weight,1,ns,1,1).*Pr0.*...
               s_ijt_ccp_up_to_scale,2);%J*1*G*T

          exp_delta=S_jt_data./s_jt_up_to_scale;%J*1*G*T
          delta=log(exp_delta);%J*1*G*T

          s_ijt_ccp=s_ijt_ccp_up_to_scale.*...
              exp_delta;%J*ns*G*T
          s_i0t_ccp=1-sum(s_ijt_ccp,[1,3]);%1*ns*1*T
          Pr0=cat(4,ones(1,ns,1,1),Pr0(:,:,:,1:T-1).*s_i0t_ccp(:,:,:,1:T-1));%1*ns*1*T

    end % if statement

end
