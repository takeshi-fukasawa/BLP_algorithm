%%% If initial values far from true values,
%%% slow conv even with Newton method...

V_initial=3*ones(2,1);
P=eye(2);
P=[0.5,0.5;0.5,0.5];
%P=[0.7,0.3;0.2,0.8];
%%P=[0.7,0.3;0.3,0.7];

delta=0;
beta_C=0.99;

v_0=[0;0.0];

ITER_MAX=1000;
for iter=1:ITER_MAX
    prob_0=exp(v_0+beta_C*P*V_initial)./exp(V_initial);
    difference=log(exp(beta_C*P*V_initial)+exp(delta))-V_initial;
    diff_V=beta_C*P.*prob_0;

    V_updated=log(exp(beta_C*P*V_initial)+exp(delta));
    
    
    %%% Derivative based method
    alpha_vec=1./(1-beta_C.*prob_0);
    alpha_vec=1./(1-beta_C*P(1,1).*prob_0(1));
    %V_updated=V_initial+difference.*alpha_vec;

    %% Newton method (not correct??)
    %%V_updated=V_initial+difference;
    %V_updated=inv(eye(2)-diff_V).*(V_initial-V_updated);

    DIST=max(abs(V_updated-V_initial));
    if DIST<1e-10
        DIST
        iter
        break;
    else
        V_initial=V_updated;
    end
end
V_updated