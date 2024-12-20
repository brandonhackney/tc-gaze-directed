function output = BIC(n,k,SSE)
% Calculates a Bayesian Information Criterion from the given values:
% n = sample size, k = # of model parameters, SSE = sum of squared errors
%
% Null models have no parameters(everything is held constant)
% The model with the lowest BIC and number of parameters is the best
% A marginal reduction in BIC at the cost of an extra parameter is bad

output = (n * log(SSE / n)) + (k * log(n));
% NOTE: the function log() returns the NATURAL log (i.e. ln)
% the standard base-10 log is log10(), and there is no ln()
end