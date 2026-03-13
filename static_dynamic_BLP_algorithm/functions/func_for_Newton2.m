function [resid,Jac]=func_for_Newton2(delta,weight,mu_ijt,...
    S_jt_data)

    global s_jt_predict

    %%% resid=S_jt_data-s_jt_predict;
    
    %%% T==1 only !! %%%%
      [J,ns,~,T]=size(mu_ijt);
      u_ijt=delta+mu_ijt;%J*I*1*T

    numer=exp(u_ijt);%J*ns*1*T
    denom=sum(numer,1)+1;%1*ns*1*T
    s_ijt=numer./denom;%J*ns*1*T

    s_jt_predict=sum(s_ijt.*weight,2);%J*1*1*T
    resid=S_jt_data-s_jt_predict;%J*1*1*T


    if nargout>=2
        Jac=...
            -diag(s_jt_predict)+...
            (sum(reshape(weight,1,1,ns).*reshape(s_ijt,J,1,ns).*...
            reshape(s_ijt,1,J,ns),3));%J*J
    end
    
end
