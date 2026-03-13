function [t_mapping,t_mapping_additional]=check_mapping_comp_cost(...
    mapping_spec,delta,mu_ijt_est,...
    weight,S_jt_data,S_0t_data)
    %% Check comp cost (Newton update vs BLP contraction mapping)
    n_trial=100;
    t_mapping_additional=0;
    tic
    for trial=1:n_trial
        [J,ns,~,T]=size(mu_ijt_est);
        u_ijt=delta+mu_ijt_est;%J*I*1*T
        
        numer=exp(u_ijt);%J*ns*1*T
        denom=sum(numer,1)+1;%1*ns*1*T
        s_ijt=numer./denom;%J*ns*1*T
        
        s_jt_predict=sum(s_ijt.*weight,2);%J*1*1*T
        resid=log(S_jt_data)-log(s_jt_predict);%J*1*1*T
    
        if mapping_spec==1
            delta_updated=delta+resid;
        else % Newton update
            Jac=...
                -eye(J)+...
                (sum(reshape(weight,1,1,ns).*reshape(s_ijt,J,1,ns).*...
                reshape(s_ijt,1,J,ns),3))./reshape(S_jt_data,J,1,1,1);%J*J
                delta_updated_Newton=delta-Jac\resid;
        end
    end

     t_mapping=toc/n_trial;

     if mapping_spec==1
         tic
         for trial=1:n_trial
             s_0t_predict=1-sum(s_jt_predict);%J*1*1*T
             delta_updated=delta_updated-(log(S_0t_data)-log(s_0t_predict));
         end
         t_mapping_additional=toc/n_trial;
     end
    



end
