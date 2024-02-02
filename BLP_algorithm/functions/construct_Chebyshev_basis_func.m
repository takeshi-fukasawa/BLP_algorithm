 function basis=...
     construct_Chebyshev_basis_func(...
     x,n_dim_Chebyshev)

    n_obs=size(x,1);

    basis=zeros(n_obs,n_dim_Chebyshev);
    basis(:,1)=1;
    basis(:,2)=reshape(x,n_obs,1);
    for i=3:n_dim_Chebyshev
       basis(:,i)=2*x.*basis(:,i-1)-basis(:,i-2);
    end
   
end
