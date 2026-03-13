options = optimoptions('fsolve',...
    'StepTolerance', TOL, ... 
    'FunctionTolerance',(1e-12)^2,...
    'Algorithm','levenberg-marquardt',...
    'Display', 'off',...
    'SpecifyObjectiveGradient',true);


func_for_Newton_anony=@(delta) func_for_Newton(delta, weight, mu_ijt_est, S_jt_data);

tic
[delta_sol_fsolve, fval, exitflag_fsolve, output_fsolve] = fsolve(...
    func_for_Newton_anony, ...
    delta_initial0, ...
    options);
t_LM=toc;

if exitflag_fsolve>=1
    conv=1;
else
    conv=0;
end

[s_jt_predict_fsolve,~]=...
      share_func(delta_sol_fsolve+mu_ijt_est,zeros(1,ns,1,T),rho_est,weight);%J*1*G*T

    DIST=max(abs(log(S_jt_data)-log(s_jt_predict_fsolve)));

    small_DIST=(DIST<TOL);


    result_LM(m,:)=[NaN,NaN,t_LM,...
            conv,log10(DIST),small_DIST];
    
