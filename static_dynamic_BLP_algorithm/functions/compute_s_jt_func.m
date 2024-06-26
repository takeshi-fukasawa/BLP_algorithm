function [s_jt_predict,Pr0]=...
             compute_s_jt_func(s_ijt_ccp,s_i0t_ccp,weight)
    
    %%% Used in share_func etc.
    global Pr0_spec
    
    [J,ns,G,T]=size(s_ijt_ccp);
    
    if Pr0_spec==0
        s_jt_predict=sum(s_ijt_ccp.*weight,2);%J*1*G*T
        Pr0=ones(1,ns,1,T);
    
    else
        [J,ns,G,T,n_state]=size(s_ijt_ccp);

        %%%%%%%%%%%%
        %%% Due to minute numerical errors, s_i0t_ccp<0??
        %%s_i0t_ccp=1-sum(s_ijt_ccp,[1,3]);%1*ns*1*T*n_state;
        %%%%%%%%%%%
        
        %%% For loop version (to be consistent with V_update_func spec)
        s_jt_predict=NaN(J,1,G,T);
        Pr0=ones(1,ns,1,T);
        for t=1:T
            s_jt_predict(:,:,:,t)=sum(...
                Pr0(:,:,:,t).*weight.*s_ijt_ccp(:,:,:,t),2);%J*1*G*1
            Pr0(:,:,:,t+1)=Pr0(:,:,:,t).*s_i0t_ccp(:,:,:,t);
        end

        %%% Vectorized version (Fast, but not comparable with V_update_func
        %%% (use for loop...)
        %cumprod_s_i0t_ccp=cumprod(s_i0t_ccp,4);%1*ns*1*T*n_state
        %Pr0=cat(4,ones(1,ns,1,1,n_state),...
        %    cumprod_s_i0t_ccp(:,:,:,1:end-1,:));%1*ns*1*T*n_state
        %s_jt_predict=sum(Pr0.*weight.*s_ijt_ccp,2);%J*1*G*T

    end

end
     