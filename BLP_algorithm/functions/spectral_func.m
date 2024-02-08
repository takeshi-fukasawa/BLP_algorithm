function [x_sol_cell,other_output_k,DIST_table,fun_k_cell]=...
    spectral_func(fun,n_var,vec,dampening_param,varargin)
%%%function [x_sol_cell,DIST_table,fun_k_cell]=...
%%%        spectral_func(fun,n_var,vec,dampening_param,varargin)

%%% Allow output of additional vars %%%%
%%% Input %%%%
% x1_0,x2_0,...: initial value
% vec: dimension of a variable whose updating tune parameters  should be independently chosen (e.g. time)
%%% vec==0 => Implement standard fixed point iteration
% dampening_param

global DEBUG FLAG_ERROR DIST count ITER_MAX TOL
global k

alpha_0=1;
alpha_0=1e-1; %% large alpha_0 lead to divergence or slow convergence...
%alpha_0=1e-2; %% large alpha_0 lead to divergence or slow convergence...

if isempty(vec)==0
    if sum(vec(:))==0
         alpha_0=0.1;
    end
end


ITER_MAX=3000;
%%ITER_MAX=118;

% varargin:1*XXX

fun_k_cell={};

%% Read inputs
other_input_cell=varargin(n_var+1:length(varargin));

for i=1:n_var
   x_0_cell{1,i}=varargin{1,i};
end

DIST_table=NaN(ITER_MAX,n_var);

%% Compute fun_0,fun_1

[fun_0_cell,other_output_0]=fun(x_0_cell{:},other_input_cell{:});
%%%fun_0_cell=fun(x_0_cell{:},other_input_cell{:});

for i=1:n_var
    x_1_cell{1,i}=x_0_cell{1,i}-alpha_0.*fun_0_cell{1,i};
    DIST_table(1,i)=max(abs(fun_0_cell{1,i}(:)),[],'all','omitnan');
end
    DIST=max(DIST_table(1,:));

if ITER_MAX>1 & DIST>TOL

[fun_1_cell,other_output_1]=fun(x_1_cell{:},other_input_cell{:});
%%%%fun_1_cell=fun(x_1_cell{:},other_input_cell{:});

%% Update x_k
x_k_cell=x_1_cell;
x_k_minus_1_cell=x_0_cell;
fun_k_cell=fun_1_cell;
fun_k_minus_1_cell=fun_0_cell;

other_output_k=[];

%% Loop

eps_val=1e-14;%%%%%%
eps_val=0;

for k=1:ITER_MAX
   DIST_vec=ones(1,n_var);
  for i=1:n_var
     Delta_x_i=x_k_cell{i}-x_k_minus_1_cell{i};
     Delta_fun_i=fun_k_cell{i}-fun_k_minus_1_cell{i};

    %%vec=4*ones(n_var,1);%%%%%
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
     else % vec==0
        sign_i=1;
        numer_i=1;
        denom_i=1;
        %%%%warning("XXX")
     end

    alpha_k_i=sign_i.*numer_i./denom_i;%scalar or vector (wrt the dimension specified in "vec")

    alpha_k_i((isnan(alpha_k_i)==1))=1;%%%
    alpha_k_i((isinf(alpha_k_i)==1))=1;%%%
    alpha_k_i(((alpha_k_i==0)))=1;%%%

    
    %alpha_k_i=1;
    
    if isempty(dampening_param)==0
        alpha_k_i=alpha_k_i*dampening_param(i);
    end

    
    x_k_plus_1_cell{1,i}=x_k_cell{1,i}-alpha_k_i.*fun_k_cell{1,i};
        
   DIST_vec(1,i)=max(abs(fun_k_cell{1,i}),[],'all','omitnan');

    end % for loop wrt i 

    DIST=nanmax(DIST_vec);
    DIST_PAST=DIST;
    DIST_table(k,:)=DIST_vec;
    
    if isnan(DIST)==1|isinf(sum(DIST))==1|isnan(sum(DIST))==1
        break;
    end


    interval=100;
    if k-floor(k/interval)*interval==0&DEBUG==1
        DIST_vec
    end


    if DIST<TOL
        if DIST==0
            warning("zero convergence")
        end
        %DIST
      FLAG_ERROR=0;
        break;
    end
    
    %%% UPDATE VARIABLES
	x_k_minus_1_cell=x_k_cell;
	x_k_cell=x_k_plus_1_cell;
	fun_k_minus_1_cell=fun_k_cell;

   [fun_k_cell,other_output_k]=fun(x_k_cell{:},other_input_cell{:});
   %%%fun_k_cell=fun(x_k_cell{:},other_input_cell{:});
   
    
end %% end of for loop wrt iter=1:ITER_MAX
count=k;
%DIST_vec

%% Output
x_sol_cell=x_k_plus_1_cell;

if isempty(other_output_k)==1
    other_output_k=other_output_1;
end

else % no further iteration case
    x_sol_cell=x_1_cell;
    other_output_k=other_output_0;
    DIST_MAT=[];
    fun_k_cell=fun_0_cell;
end

return


