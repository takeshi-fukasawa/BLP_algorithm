% Function to compute market shares s_j given u_ijt_tilde
function [s_jt_predict,ChoiceProb,s_ijt_given_g_ccp,s_igt_ccp,...
    numer_1,denom_1,numer_2,denom_2]=...
    share_func(u_ijt_tilde,u_i0t_tilde,rho,weight)

    numer_1=exp(u_ijt_tilde./(1-rho));%J*ns*G*T
    denom_1=sum(numer_1,1);%1*ns*G*T
    numer_2=exp((1-rho).*log(denom_1));%1*ns*G*T
    denom_2=exp(u_i0t_tilde)+sum(numer_2,3);%1*ns*1*T
    s_ijt_given_g_ccp=numer_1./denom_1;%J*ns*G*T
    s_igt_ccp=numer_2./denom_2;%J*ns*G*T; prob of purchasing anything
    ChoiceProb=s_ijt_given_g_ccp.*s_igt_ccp;%J*ns*G*T
    s_jt_predict=sum(ChoiceProb.*weight,2);%J*ns*G*T

    s_0t_predict=1-sum(s_jt_predict,[1,3]);%1*ns*1*T

end
