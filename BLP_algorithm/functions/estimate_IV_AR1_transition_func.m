function [coef_0,coef_1,y,y_predict,R2]=estimate_IV_AR1_transition_func(IV_state_obs_pt)

    %%%% Estimate AR1 coefficients of IV state transitions %%%%%%
    X=IV_state_obs_pt(:,:,:,1:end-1);%1*ns*G*(T-1)
    y=IV_state_obs_pt(:,:,:,2:end);%1*ns*G*(T-1)
    X_mean=mean(X,4);%1*ns*G*1
    y_mean=mean(y,4);%1*ns*G*1
    Var_X=sum((X-X_mean).^2,4);%1*ns*G*1; not divided by the number of elements; not normalized
    Cov_X_y=sum((X-X_mean).*(y-y_mean),4);%1*ns*G*1; not divided by the number of elements; not normalized
    coef_1=Cov_X_y./Var_X;%1*ns*G*1
    coef_0=y_mean-coef_1.*X_mean;%1*ns*G*1

    y_predict=X.*coef_1+coef_0;%1*ns*G*(T-1)
    
    %R2=corr(y_predict(:),y(:));
    R2=[];
    
end