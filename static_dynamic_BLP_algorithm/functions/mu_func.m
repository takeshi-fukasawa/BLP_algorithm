function mu_ijt=mu_func(sigma_vec,v,X_jt)
   [J,~,G,T,K]=size(X_jt);
   ns=size(v,2);
   mu_ijt=sum(reshape(sigma_vec,1,1,1,1,K).*...
       reshape(v,1,ns,1,1,K).*...
       reshape(X_jt,J,1,G,T,K),5);%J*ns*G*T

end
