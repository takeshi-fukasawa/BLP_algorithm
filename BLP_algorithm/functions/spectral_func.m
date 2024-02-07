%%%function [x_sol_cell,other_output_k,DIST_table,fun_k_cell]=...
%%%    spectral_func(fun,n_var,vec,dampening_param,varargin)
function [x_sol_cell,DIST_table,fun_k_cell]=...
        spectral_func(fun,n_var,vec,dampening_param,varargin)

%%% Allow output of additional vars %%%%
%%% Input %%%%
% x1_0,x2_0,...: initial value
% vec: dimension of a variable that should be independently updated (e.g. time)
% dampening_param

global DEBUG FLAG_ERROR DIST count ITER_MAX TOL
global k

alpha_0=1;
alpha_0=1e-1; %% large alpha_0 lead to divergence or slow convergence...
%alpha_0=1e-2; %% large alpha_0 lead to divergence or slow convergence...

ITER_MAX=3000;

% varargin:1*XXX

fun_k_cell={};

%% Read inputs
other_input_cell=varargin(n_var+1:length(varargin));

for i=1:n_var
   x_0_cell{1,i}=varargin{1,i};
end

%% Derive fun_0,fun_1

%%%[fun_0_cell,other_output_0]=fun(x_0_cell{:},other_input_cell{:});
fun_0_cell=fun(x_0_cell{:},other_input_cell{:});

for i=1:n_var
    x_1_cell{1,i}=x_0_cell{1,i}-alpha_0.*fun_0_cell{1,i};
end

if ITER_MAX>1

%%%[fun_1_cell,other_output_1]=fun(x_1_cell{:},other_input_cell{:});
fun_1_cell=fun(x_1_cell{:},other_input_cell{:});

%% Update x_k
x_k_cell=x_1_cell;
x_k_minus_1_cell=x_0_cell;
fun_k_cell=fun_1_cell;
fun_k_minus_1_cell=fun_0_cell;

other_output_k=[];

%% Loop
DIST_table=zeros(ITER_MAX,n_var);
eps_val=0;
for k=1:ITER_MAX

   DIST_vec=ones(1,n_var);
  for i=1:n_var
     delta_x_i=x_k_cell{i}-x_k_minus_1_cell{i};
     delta_fun_i=fun_k_cell{i}-fun_k_minus_1_cell{i};

    eps_val=1e-14;%%%%%%
    eps_val=0;%%%%%%%

    vec=4*ones(n_var,1);
     if isempty(vec)==1
       sign_i=sign(sum(delta_x_i.*delta_fun_i,'all','omitnan'));%scalar
       numer_i=sqrt(sum(delta_x_i.^2,'all','omitnan'))+eps_val;
       denom_i=sqrt(sum(delta_fun_i.^2,'all','omitnan'))+eps_val;

     elseif isempty(vec)==0
      sum_delta_x_fun=sum(delta_x_i.*delta_fun_i,1:vec(i)-1,'omitnan');
      sum_delta_x_x=sum(delta_x_i.^2,1:vec(i)-1,'omitnan');
      sum_delta_fun_fun=sum(delta_fun_i.^2,1:vec(i)-1,'omitnan');

       sign_i=sign(sum_delta_x_fun);%1*1*n_dim etc.
       numer_i=sqrt(sum_delta_x_x)+eps_val;
       denom_i=sqrt(sum_delta_fun_fun)+eps_val;

       %%% Nonstationary: Not implementing the following code => Faster
       %%% convergence (Dynamic BLP nonstationary code)
       %sign_i(2:end)=sign(sum(sum_delta_x_fun(2:end),'all','omitnan'));%scalar
       %numer_i(2:end)=sqrt(sum(sum_delta_x_x(2:end),'all','omitnan'))+eps_val;
       %denom_i(2:end)=sqrt(sum(sum_delta_fun_fun(2:end),'all','omitnan'))+eps_val;

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
    if k-floor(k/interval)*interval==0&DEBUG==1&n_var>=2
        DIST_vec
    end


    if DIST<TOL
       if DIST_PAST>TOL*10
           warning("Divergence")
           FLAG_ERROR=1;
       else
           FLAG_ERROR=0;
       end   
        break;
    end
    
	x_k_minus_1_cell=x_k_cell;
	x_k_cell=x_k_plus_1_cell;
	fun_k_minus_1_cell=fun_k_cell;

   %%%[fun_k_cell,other_output_k]=fun(x_k_cell{:},other_input_cell{:});
   fun_k_cell=fun(x_k_cell{:},other_input_cell{:});
   
    
end %% end of for loop wrt iter:1:ITER_MAX
count=k;
%DIST_vec

%% Output
x_sol_cell=x_k_plus_1_cell;

%%%if isempty(other_output_k)==1
%%%    other_output_k=other_output_1;
%%%end

elseif ITER_MAX==1
    x_sol_cell=x_1_cell;
    %%%other_output_k=other_output_0;
    DIST_MAT=[];
    fun_k_cell=fun_0_cell;
end

return


