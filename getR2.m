function [R2, SSEout] = getR2(actual, predicted)
% Get a "signed-squared correlation", or R^2, based on predicted data
% Instead of using the Matlab regression functions like a normal person,
% we're manually calculating regression based on dataset 1,
% then using the betas to generate a predicted dataset,
% and finally comparing that prediction to a dataset 2.
% I'm not sure the built-ins allow you to do that;
% they just seem to measure how well your model fits dataset 1.
%
% Either way, SSE is sum of squared ERROR (from prediction),
% and SST is sum of squares TOTAL (i.e. from the sample mean)
% I'm assuming the mean here is within each voxel? So just down each row.
% The result should be a unique R2 per voxel,
% signifying how well your prediction explains variance in actual data.

assert(all(size(actual) == size(predicted)), 'Inputs are not the same size! Matrix 1 is %s while matrix 2 is %s', num2str(size(actual)), num2str(size(predicted))); 

% SSE = sum((actual - predicted) .^2, 1);
% SST = sum((actual - mean(actual, 1)).^2, 1);
% R2 = 1 - (SSE ./ SST);

% A second definition is based on variance explained:
% SSR is the "regression" SSE, i.e. deviance of prediction from avg
% Then SSR / SST is the proportion of variance explained by the model
% But each must be scaled by the number of observations
n = height(actual);
SSR = sum((predicted - mean(actual, 1)) .^2, 1) / n;
% SST = sum((actual - mean(actual, 1)).^2, 1) / n;
% R2 = SSR ./ SST;
SSEout = sum(SSR, 'all');

% The alternative is to literally square the r value, i.e. Pearson corr
% A "signed-squared" value preserves the original sign.
% This is the method used by McMahon Bonner & Isik.
% A strong negative R2 is considered worse than a weak positive R2.
R = columncorr(actual, predicted);
s = sign(R);
R2 = s .* (R.^2);

% OPTIONAL OUTPUT: return the SSE so you can calculate a BIC
end