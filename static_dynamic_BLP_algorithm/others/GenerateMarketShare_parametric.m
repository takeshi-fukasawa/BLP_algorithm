function [diff_ms,dev_diff_ms] = GenerateMarketShare_parametric(...
    mean_utility,NbsSample,random_coef,X2,theta2,market_share)
    %% parametric implementation of GMM
    %% Replication code of Iaria and Wang (2025)
    %% Modified by Takeshi Fukasawa in April 2026

    % three parameters, (th1 th2; 0 th3)*
    
    random_coef_para = theta2*random_coef;%NbsX2*NbsSample
    weight=ones(NbsSample,1)/NbsSample;
    
    for iter=1:100 % To precisely measure the comp time, repeat the computation 100 times
        exp_utility_all=exp(repmat(mean_utility,1,NbsSample)+X2*random_coef_para);%J*NbsSample
        
        ms_all=exp_utility_all*diag(1./(sum(exp_utility_all,1)+1));%J*Sample

        %% More efficient code
        ms_all2=exp_utility_all./(sum(exp_utility_all,1)+1);

        ms=ms_all*weight;
        diff_ms = ms-market_share;%J*1
        
        if nargout >1
            %% Compute a term corresponding to the Jacobian
            dev_diff_ms=diag(ms)-(ms_all*diag(weight)*ms_all');

            %% More efficient code
            J=size(ms,1);
            dev_diff_ms2=diag(ms)-(sum(reshape(weight,1,1,NbsSample).*reshape(ms_all,J,1,NbsSample).*...
                    reshape(ms_all,1,J,NbsSample),3));
        end
    end
end

