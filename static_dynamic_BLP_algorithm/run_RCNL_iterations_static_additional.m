result_hybrid_1=NaN(n_market,6);
result_hybrid_2=NaN(n_market,6);
result_LM=NaN(n_market,6);

for m=1:n_market
    market_id=m
    
    rng(m);

    run DGP.m
    

    %% Compute product market share, given parameters
    run solve_equil.m

    if J<100
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


if J==25
    table_hybrid_1=mean(result_hybrid_1);
    table_LM=mean(result_LM);
    table=[table_hybrid_1;table_LM];
else
    table_LM=mean(result_LM);
    table=[table_hybrid_1;table_LM];
end


