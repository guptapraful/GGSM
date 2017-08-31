function [sigma beta] = estimate_params_mv_ggd_mm(X)

% Function to compute the parameters of a MVGGD
% X is a Nxn matrix, where each row is a sample and the N rows are assumed
% iid
[N n] = size(X);
 
S = cov(X);
% estimate beta
gamma2_cap =  sum(sum((X*pinv(S)).*X,2).^2)/N ;  
fun = @(beta) n^2.*gamma(n./(2.*beta)).*gamma((n+4)./(2.*beta))./(gamma((n+2)./(2.*beta) )).^2 - gamma2_cap;
bb = 0.05:0.01:10;
tt = abs(fun(bb));
[minval ind] = min(tt);
beta = bb(ind);

sigma = S*n*gamma(n/(2*beta))/(2^(1/beta)*gamma((n+2)/(2*beta)));
