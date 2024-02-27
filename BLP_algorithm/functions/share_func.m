% Function to compute market shares s_j given u_ijt_tilde
function [s_jt_predict,s_ijt_ccp,s_ijt_given_g_ccp,s_igt_ccp,...
    numer_1,denom_1,numer_2,denom_2]=...
    share_func(v_ijt_tilde,v_i0t_tilde,rho,weight)

    global Pr0_spec

    %%% weight:1*ns*1*1*n_state; sum of weight need not be equal to 1
    numer_1=exp(v_ijt_tilde./(1-rho));%J*ns*G*T*(1 or n_state)
    denom_1=sum(numer_1,1);%1*ns*G*T*n_state
    numer_2=exp((1-rho).*log(denom_1));%1*ns*G*T*n_state; exp of IV(inclusive value)
    denom_2=exp(v_i0t_tilde)+sum(numer_2,3);%1*ns*1*T*n_state
    s_ijt_given_g_ccp=numer_1./denom_1;%J*ns*G*T*n_state
    s_igt_ccp=numer_2./denom_2;%1*ns*G*T*n_state; prob of purchasing anything
    s_ijt_ccp=s_ijt_given_g_ccp.*s_igt_ccp;%J*ns*G*T*n_state

    [s_jt_predict]=...
             compute_s_jt_func(s_ijt_ccp,weight);

end
