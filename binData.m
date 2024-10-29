function data2 = binData(time1, data1, time2)
% With the expectation that data1 was sampled at a higher rate than data2,
% dimension-reduce data1 to the timescale of time2 by averaging.
% For example, if time1 is at 60Hz and time2 is at 3Hz,
% each element of the data2 represents 20 samples of the data1.
% Do this by averaging data1(1:20), data1(21:40), etc.

% This is not the same as a histogram because we're not tallying, we're
% averaging.

numItems = length(time2);
data2 = zeros(size(time2));

for i = 1:numItems
    % Find everything b/w current time and next time
    % (with an exception for the final element, when there is no next time)
    starttime = time2(i);
    if i == numItems
        endtime = max(time1) + 1;
    else
        endtime = time2(i+1);
    end
    % time2 is in sec, but time1 is in msec, so rescale
%     starttime = starttime * 1000;
%     endtime = endtime * 1000;
    % Get the data
    subset = time1 >= starttime & time1 < endtime;
    data2(i) = max(data1(subset));
end