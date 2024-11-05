function r2 = columncorr(actual, predicted)
% Given two matrices, correlate the paired columns (ie NOT a full pairwise)
% e.g. correlate A1 to B1, A2 to B2, etc. but NOT A1 to B2
A_centered = actual - mean(actual, 1);
B_centered = predicted - mean(predicted, 1);

numerator = sum(A_centered .* B_centered, 1); % sum of error products
denominator = sqrt(sum(A_centered .^2, 1)) .* sqrt(sum(B_centered .^2, 1)); % product of root SSEs

r2 = numerator ./ denominator;