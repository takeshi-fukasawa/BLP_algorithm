%%% Static case
numer=reshape(weight,1,ns).*reshape(exp(mu_ijt_est),J,ns,G,T).*reshape(exp(-V_true),1,ns,1,T);%J*ns*G*T
denom=sum(numer,2);%J*1*G*T
prob_i_given_j=numer./denom;%J*ns*G*T

numer=reshape(weight,1,ns).*reshape(exp(-V_true),1,ns,1,T);%J*ns*1*T
denom=sum(numer,2);%J*1*1*T
prob_i_given_0=numer./denom;%J*ns*1*T

diff_prob=sum(abs(prob_i_given_j-prob_i_given_0),2);
diff_prob_max=max(diff_prob);
