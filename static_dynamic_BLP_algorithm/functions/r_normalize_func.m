
%%%%%%%%%%%%%%%%%%%%%%%%%

function r_out=r_normalize_func(r,S_0_data)
    %r_out=r;
    %temp=S_0_data-sum(exp(r(:,1:end-1)),2);%1*1
    %if temp>0
    %    r_out(:,end)=log(temp);
    %else

    %r_out=log(S_0_data*exp(r)./sum(exp(r),2));%1*ns
    %end

    %r_out(:,end)=log(temp);

    %%% Always scale; faster??
    r_out=log(S_0_data.*exp(r)./sum(exp(r),2));%1*ns

    %S_0_data-sum(exp(r_out),2) %zero??
end