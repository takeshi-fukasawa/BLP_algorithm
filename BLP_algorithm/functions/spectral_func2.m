function [x_sol_cell,other_output_k,DIST_table,iter_info,fun_k_cell]=...
    spectral_func2(fun,n_var,vec,dampening_param,varargin)
%%%function [x_sol_cell,DIST_table,fun_k_cell]=...
%%%        spectral_func(fun,n_var,vec,dampening_param,varargin)

%%% Allow output of additional vars %%%%
%%% Input %%%%
% x1_0,x2_0,...: initial value
% vec: dimension of a variable whose updating tune parameters  should be independently chosen (e.g. time)
%%% vec==0 => Implement standard fixed point iteration
% dampening_param

global DEBUG FLAG_ERROR DIST count ITER_MAX TOL
global k alpha_table

tic

alpha_0=1;
alpha_0=1e-1; %% large alpha_0 lead to divergence or slow convergence...

if isempty(vec)==0
    if sum(vec(:))==0
         alpha_0=0.1;
    end
end

%%alpha_0=1e-2; %% large alpha_0 lead to divergence or slow convergence...


ITER_MAX=3000;
%%ITER_MAX=118;

line_search_spec=0;


ITER_MAX_LINE_SEARCH=1;
rho=0.5;
M=10;gamma=10^(-4);

if line_search_spec==0
    ITER_MAX_LINE_SEARCH=1;
end

% varargin:1*XXX

fun_k_cell={};


%% Read inputs
other_input_cell=varargin(n_var+1:length(varargin));

for i=1:n_var
   x_0_cell{1,i}=varargin{1,i};
end

DIST_table=NaN(ITER_MAX,n_var);
alpha_table=NaN(ITER_MAX,n_var);

   [fun_0_cell,other_output_0]=fun(x_0_cell{:},other_input_cell{:});
    
    %%% DIST: sup norm of F(x)=x-Phi(x). 
    DIST_vec=ones(1,n_var);
    for i=1:n_var
      DIST_vec(1,i)=max(abs(fun_0_cell{1,i}),[],'all','omitnan');
    end

    DIST=nanmax(DIST_vec);
    DIST_table(1,:)=DIST_vec;
    alpha_table(1,:)=alpha_0;

x_k_cell=x_0_cell;
fun_k_cell=fun_0_cell;
other_output_k=other_output_0;

%%%%%%%% Loop %%%%%%%%%%%

eps_val=1e-14;%%%%%%
eps_val=0;

if DIST>TOL
for k=0:ITER_MAX-1

   %%% In the current iteration, fun_k_cell is assumed to be far from zero


    for i=1:n_var
     if k>=1
     Delta_x_i=x_k_cell{i}-x_k_minus_1_cell{i};
     Delta_fun_i=fun_k_cell{i}-fun_k_minus_1_cell{i};

    %%vec=4*ones(n_var,1);%%%%%

    %%%%%%%%%%%%%%%
    sum_Delta_x_fun=sum(Delta_x_i.*Delta_fun_i,'all','omitnan');%vector
      sum_Delta_x_x=sum(Delta_x_i.^2,'all','omitnan');
      sum_Delta_fun_fun=sum(Delta_fun_i.^2,'all','omitnan');
%%%%%%%%%%%%%%%%%

    if isempty(vec)==1 

       sign_i=sign(sum(Delta_x_i.*Delta_fun_i,'all','omitnan'));%scalar
       numer_i=sqrt(sum(Delta_x_i.^2,'all','omitnan'))+eps_val;
       denom_i=sqrt(sum(Delta_fun_i.^2,'all','omitnan'))+eps_val;

     elseif isempty(vec)==0 & sum(vec(:))>0 %%% XXX-dependent tune parameters
      sum_dim_ids=1:size(vec(:),1);
      sum_dim_ids=sum_dim_ids(sum_dim_ids~=vec(i));

      sum_dim_ids=[1:3,5];%%%%%

      sum_Delta_x_fun=sum(Delta_x_i.*Delta_fun_i,sum_dim_ids,'omitnan');%vector
      sum_Delta_x_x=sum(Delta_x_i.^2,sum_dim_ids,'omitnan');
      sum_Delta_fun_fun=sum(Delta_fun_i.^2,sum_dim_ids,'omitnan');

       sign_i=sign(sum_Delta_x_fun);%1*1*n_dim etc.
       numer_i=sqrt(sum_Delta_x_x)+eps_val;
       denom_i=sqrt(sum_Delta_fun_fun)+eps_val;

     end

     %%%%%%%%%%%%%%
     %%% BB (otherwise, Varadnahn Roland 2008)
     %%% Worse performance???
     %sign_i=1;%%%%%

     %numer_i=sum_Delta_x_x;
     %denom_i=sum_Delta_x_fun;

     %%% BB second
     %%% Better than the first, but worse than the VR??
     %numer_i=sum_Delta_x_fun;
     %denom_i=sum_Delta_fun_fun;

     
    %%%%%%%%%%%%%%%%%%%%

     if vec==0
        sign_i=1;
        numer_i=1;
        denom_i=1;
     end

    sign_i=1;%%%%%% S3 spec (const sign) %%%%
    
    alpha_k_i=sign_i.*numer_i./denom_i;%scalar or vector (wrt the dimension specified in "vec")

    alpha_k_i((isnan(alpha_k_i)==1))=1;%%%
    alpha_k_i((isinf(alpha_k_i)==1))=1;%%%
    alpha_k_i(((alpha_k_i==0)))=1;%%%

    %%%%%%%%%%%%%%%%%%%
    if (abs(alpha_k_i)<1e-2 |abs(alpha_k_i)>100) & vec~=0
        %warning("too small")
    end
    %%%%%%%%%%%%%%%

   else%% k==1
     alpha_k_i=alpha_0;  
   end

    %alpha_k_i=1;
    
    if isempty(dampening_param)==0
        alpha_k_i=alpha_k_i*dampening_param(i);
    end

      alpha_k{1,i}=alpha_k_i;

      alpha_table(k+1,i)=mean(alpha_k_i(:));

   end % for loop wrt i


   %%% Update variables %%%%%%%%%%%%%%%
   for n=1:ITER_MAX_LINE_SEARCH
    for i=1:n_var
        x_k_plus_1_cell{1,i}=x_k_cell{1,i}-alpha_k{1,i}.*fun_k_cell{1,i};
    end
       [fun_k_plus_1_cell,other_output_k_plus_1]=fun(x_k_plus_1_cell{:},other_input_cell{:});
  
    %%% DIST: sup norm of F(x)=x-Phi(x). 
    DIST_vec=ones(1,n_var);
    for i=1:n_var
      DIST_vec(1,i)=max(abs(fun_k_plus_1_cell{1,i}),[],'all','omitnan');
    end

    
    if line_search_spec==1
       DIST_PAST_MAX_vec=max(DIST_table(max(1,k+1-M+1):k+1,:),[],1);

        id=1;
        eta_k=DIST_table(1,id)/((1+k)^2);
        LHS=DIST_vec(id);
        RHS=DIST_PAST_MAX_vec(id)+eta_k-gamma*(sum(alpha_k{1,id}(:)).^2)*DIST_table(k+1,id);

        if LHS>RHS
            alpha_k_i=alpha_k_i.*rho;
        else
            break;
        end
    end% line search spec==1

    end % end n=1:ITER_MAX_LINE_SEARCH loop
       
    DIST_table(k+2,:)=DIST_vec;
    DIST=nanmax(DIST_vec);

    if isnan(DIST)==1|isinf(sum(DIST))==1|isnan(sum(DIST))==1
       %warning("Error ?? ")
       x_k_plus_1_cell=x_k_cell;
       FLAG_ERROR=1;
       break;
    end

       
   interval=100;
    if k-floor(k/interval)*interval==0&DEBUG==1
        DIST_vec
    end


    if DIST<TOL
        %DIST
      FLAG_ERROR=0;
        break;
    end



    
    %%% Replace variables for the next iteration
	x_k_minus_1_cell=x_k_cell;
	x_k_cell=x_k_plus_1_cell;
	fun_k_minus_1_cell=fun_k_cell;
    fun_k_cell=fun_k_plus_1_cell;
    other_output_k=other_output_k_plus_1;

   
    
end %% end of for loop wrt k=0:ITER_MAX-1
count=k;
%DIST_vec

else % no iteration
    count=1;
    k=0;
    x_k_plus_1_cell=x_k_cell;
            x_sol_cell=x_0_cell;
        other_output_k=other_output_0;
        DIST_table=[];
        fun_k_cell=fun_0_cell;

end

%% Output
x_sol_cell=x_k_plus_1_cell;

t_cpu=toc ;
iter_info.t_cpu=t_cpu;
iter_info.n_iter=k;
iter_info.ITER_MAX=ITER_MAX;

return


