function output = averageRunResiduals(residuals, numTRs)
% The residuals you get are going to cover many runs,
% but you want to predict data for a single run,
% so... average all those other runs together??

numT = height(residuals);
numRuns = numT / numTRs; % ought to be an integer
assert(~rem(numRuns, 1), ' The math aint mathin!');
s = 1;
e = numTRs;
for i = 1:numRuns
    tmp(:,:,i) = residuals(s:e,:);
    s = e + 1;
    e = e + numTRs;
end
output = mean(tmp, 3);
