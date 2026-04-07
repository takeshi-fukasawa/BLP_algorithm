result_hybrid_1=NaN(n_market,6);
result_hybrid_2=NaN(n_market,6);
result_LM=NaN(n_market,6);

for m=1:n_market
    market_id=m
    
    rng(m);

    run DGP.m
    

    %% Compute product market share, given parameters
    run solve_equil.m

    if J<100 & 1==1
        %% BLP contraction mapping (G>=2 case allowed)
        tune_param_BLP=0;
        run run_BLP_contraction.m
        
        tune_param_BLP=1;
        run run_BLP_contraction.m
    end

    %% Try Iaria & Wang 2025's Hybrid method (FPI => Newton)
    run run_Newton.m

    %% Levenberg Marquart
    run run_LM.m
    
end

table_hybrid_1=mean(result_hybrid_1);


table_LM=mean(result_LM);
table=[table_hybrid_1;table_LM];

if J<100
    table_FPI=[...
    mean(results_BLP_0);...
    mean(results_BLP_1);...
    mean(results_BLP_1_Anderson)];
    table_FPI=[table_FPI(:,1),NaN(3,1),table_FPI(:,2:end)];
    table=[table_FPI;table];
end