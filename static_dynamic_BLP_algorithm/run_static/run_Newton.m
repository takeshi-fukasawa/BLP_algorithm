global s_jt_predict

%% Try Iaria & Wang (2025)'s hybrid algorithm
%%% Fixed-point iteration => Newton

TOL=1e-12;
alpha0=0.130707;
ITER_MAX=1000;

%%% hybrid_spec==1: Use Iaria & Wang 2025's threshold
%%% hybrid_spec==2: Use threshold 1e-2 for switching from FP to Newton


for hybrid_spec=1:1
    iter_BLP_contraction=0;
    iter_Newton=0;

    delta=delta_initial0;

    tic
    for iter=1:ITER_MAX
        [resid,Jac]=func_for_Newton(delta,weight,mu_ijt_est,...
        S_jt_data);
        DIST=max(abs(resid));
        if DIST<TOL
            break;
        else
            
            if hybrid_spec==1
                s_0t_predict=1-sum(s_jt_predict(:));
                threshold=alpha0.*s_0t_predict.*min(s_jt_predict(:))./...
                    (2.64*J)...
                    .*s_0t_predict;
                threshold=threshold.^(1/3);
                %threshold=max(threshold,1e-1);%%%%

            elseif hybrid_spec==2
                threshold=1e-2;
            end

            if DIST>threshold
                iter_BLP_contraction=iter_BLP_contraction+1;
                delta=delta+resid;
            else
                iter_Newton=iter_Newton+1;
                delta=delta-Jac\resid;
            end
        end

    end
    t_hybrid=toc;

    if DIST<TOL
        conv=1;
    else
        conv=0;
    end

    [s_jt_predict_hybrid,~]=...
        share_func(delta+mu_ijt_est,zeros(1,ns,1,T),rho_est,weight);%J*1*G*T

    DIST=max(abs(log(S_jt_data)-log(s_jt_predict_hybrid)));

    small_DIST=(DIST<TOL);

    if hybrid_spec==1
        result_hybrid_1(m,:)=[iter_BLP_contraction,iter_Newton,t_hybrid,...
            conv,log10(DIST),small_DIST];
    elseif hybrid_spec==2
        result_hybrid_2(m,:)=[iter_BLP_contraction,iter_Newton,t_hybrid,...
            conv,log10(DIST),small_DIST];
    end
end

