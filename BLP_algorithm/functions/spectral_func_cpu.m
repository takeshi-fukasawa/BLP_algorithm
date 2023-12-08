function [x_sol_cell]=spectral_func(fun,n_var,varargin)

%%% Input %%%%
% x1_0,x2_0,...: initial value
% x1_max,x2_max,...: max of x ([] implies empty)
% x1_min, x2_min,...: min of x ([] implies empty)

global DEBUG FLAG_ERROR DIST DIST_table count ITER_MAX TOL

alpha_0=1e-1; %% large alpha_0 lead to divergence or slow convergence...
%TOL=1e-11;
%ITER_MAX=100;%%300;

% varargin:1*XXX
%size(varargin)

%% Read inputs
x_elem=NaN(1,n_var); %Number of elements
for i=1:n_var
	x_size{1,i}=size(varargin{1,i});%%%
	x_elem(1,i)=prod(x_size{1,i},2);
end
x_0_mat=zeros(max(x_elem,[],2),n_var);
x_max_mat=zeros(max(x_elem,[],2),n_var);
x_min_mat=zeros(max(x_elem,[],2),n_var);

for i=1:n_var
	temp=varargin{1,i};
	x_0_mat(1:x_elem(1,i),i)=temp(:);
    if isnan(varargin{1,n_var+i})==0
        temp=varargin{1,n_var+i};
        x_max_mat(1:x_elem(1,i),i)=temp(:);
    else
        x_max_mat(:,i)=abs(x_0_mat(:,i))*100+100;
    end
    
    if isnan(varargin{1,2*n_var+i})==0
        temp=varargin{1,2*n_var+i};
        x_min_mat(1:x_elem(1,i),i)=temp(:);
    else
        x_min_mat(:,i)=-abs(x_0_mat(:,i))*100-100;
    end
end


x_0_cell=mat_to_cell_func(x_0_mat,n_var,x_elem,x_size);

%% Derive fun_0,fun_1
[fun_0_cell,fun_0_mat]=...
    fun_out(fun,x_0_cell,n_var,x_elem,varargin{1:length(varargin)});

x_1_mat=x_0_mat-alpha_0.*fun_0_mat;
x_1_cell=mat_to_cell_func(x_1_mat,n_var,x_elem,x_size);

[fun_1_cell,fun_1_mat]=fun_out(fun,x_1_cell,n_var,x_elem,varargin{1:length(varargin)});


%% Update x_k
x_k_mat=x_1_mat;
x_k_minus_1_mat=x_0_mat;
fun_k_mat=fun_1_mat;
fun_k_minus_1_mat=fun_0_mat;

DIST_table=zeros(ITER_MAX,n_var);

%% Loop
for k=1:ITER_MAX
    
	delta_x=x_k_mat-x_k_minus_1_mat;%x_elem*n_var
	delta_fun=fun_k_mat-fun_k_minus_1_mat;%x_elem*n_var
	%% L2 norm
	alpha_k=sqrt(sum(delta_x.^2,1,'omitnan'))./sqrt(sum(delta_fun.^2,1,'omitnan'));%1*n_var
    %alpha_k
    
    %sum(abs(fun_k_mat))
    alpha_k(isnan(alpha_k))=0;
    alpha_k(isinf(alpha_k))=0;
    %alpha_k=1;
    alpha_k=max(alpha_k,1e-10);% If not, alpha_k==0no update... (fun_0==0 case)

    %update_dummy=(isinf(alpha_k)==0|isnan(alpha_k)==0);
    x_k_plus_1_mat=x_k_mat-alpha_k.*fun_k_mat;
    
    
    %x_k_plus_1_mat=...
    %    (x_k_plus_1_mat<=x_max_mat & x_k_plus_1_mat>=x_min_mat & isinf(alpha_k)==0).*x_k_plus_1_mat+...
    % 	(x_k_plus_1_mat>x_max_mat&isinf(alpha_k)==0).*x_max_mat+...
    %		(x_k_plus_1_mat>x_max_mat&isinf(alpha_k)==0).*x_min_mat+...
    %   	(isinf(alpha_k)==1).*x_k_mat;
    

	DIST_MAT=max(abs(x_k_plus_1_mat-x_k_mat));%1*n_var
	DIST_table(k,:)=DIST_MAT;

    DIST=nanmax(DIST_MAT);
    DIST_PAST=DIST;
    
    if isnan(DIST)==1
        break;
    end
    if DIST<TOL
       if DIST_PAST>TOL*10 | isnan(DIST)==1
           warning("Divergence")
           FLAG_ERROR=1;
       else
           FLAG_ERROR=0;
       end   
        break;
    end

	x_k_minus_1_mat=x_k_mat;
	x_k_mat=x_k_plus_1_mat;
	fun_k_minus_1_mat=fun_k_mat;

	x_k_cell=mat_to_cell_func(x_k_mat,n_var,x_elem,x_size);
	[fun_k_cell,fun_k_mat]=...
        fun_out(fun,x_k_cell,n_var,x_elem,varargin{1:length(varargin)});

    if DEBUG==1
        DIST_MAT %%%%%
    end
    
end %% of for loop
count=k+2;
DIST;

%% Output
x_sol_cell=mat_to_cell_func(x_k_plus_1_mat,n_var,x_elem,x_size);

return


%% Nested functions
function [fun_out_cell,fun_out_mat]=fun_out(fun,x_cell,n_var,x_elem,varargin)
temp=varargin(n_var*3+1:length(varargin));

fun_out_cell=fun(x_cell{:},temp{:});

fun_out_mat=zeros(max(x_elem),n_var);
for i=1:n_var
	temp=fun_out_cell{i};
	fun_out_mat(1:x_elem(i),i)=temp(:);
end

return

function x_cell=mat_to_cell_func(x_mat,n_var,x_elem,x_size)
for i=1:n_var
	x_cell{1,i}=reshape(x_mat(1:x_elem(i),i),x_size{i});
end
return
