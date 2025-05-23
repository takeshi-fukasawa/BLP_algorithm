run preliminary_path.m

%spec_validation="BLP"
spec_validation="Nevo"

if spec_validation=="BLP"
    load(append(pyblp_results_path,'data_BLP.mat'))
    filename="BLP_c_val_results.csv";
elseif spec_validation=="Nevo"
    load(append(pyblp_results_path,'data_Nevo.mat'))
    filename="Nevo_c_val_results.csv";
    %load('data_Nevo.mat')
elseif spec_validation=="MonteCarlo"
    filename="MonteCarlo_c_val_results.csv";
end


market_ids=findgroups(market_ids);
min_market_id=min(market_ids);
max_market_id=max(market_ids);

T=max_market_id-min_market_id+1;

    c0_tilde=NaN(T,1);
    c1_tilde=NaN(T,1);
    c0_tilde2=NaN(T,1);

for t=min_market_id:max_market_id

    %if spec_validation=="BLP"| spec_validation=="Nevo"
    
        weights=weights(:);
        I=size(choice_prob,2);
        tt=t-min_market_id+1;
        ids=find(market_ids==t);
        J_t=size(ids(:),1);
        s_jt=reshape(shares(ids),J_t,1);
        weight_t=weights(I*(tt-1)+1:I*tt);
        weight_t=weight_t./sum(weight_t);%I*1
        s_ijt=choice_prob(ids,:);%J_t*I
    %else

    %end


    %%% General code
    s_i0t=1-sum(s_ijt,1);%1*I
    s_0t=1-sum(s_jt(:),1);%1*1
    term1=reshape(weight_t,1,I).*reshape(s_ijt,J_t,I)./reshape(s_jt,J_t,1);%J_t*I
    term2=reshape(weight_t,1,I).*reshape(s_i0t,1,I)./reshape(s_0t,1,1);%J_t*I
    
    
    c0_tilde(tt)=max(1-s_i0t(:));
    c1_tilde(tt)=max(s_i0t(:))-min(s_i0t(:));
    c0_tilde2(tt)=1-sum(weight_t(:).*s_i0t(:));
end

c_val_table=[mean(c0_tilde),std(c0_tilde);...
    mean(c1_tilde),std(c1_tilde)]

writematrix(c_val_table,append('results/static_BLP/',filename))




