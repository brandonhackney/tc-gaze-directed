function output = stackRuns(dataCell)
% Given a cell array with data matrices of the same width,
% stack them all on top of each other. Ugly due to preallocation.

% Since this is going to be a HUGE matrix, we'll want to preallocate,
% so do some recon to figure out how big the final result will be.
runLengths = [];
widths = [];
for i = 1:length(dataCell)
    runLengths(i) = height(dataCell{i});
    widths(i) = width(dataCell{i}); % ought to all be equal
end
numRows = sum(runLengths);
assert(range(widths) == 0, 'Some of these runs have a different number of vertices! Aborting')

% Now go through and drop everything in.
output = zeros(numRows, widths(1));
start = 1; % initialize
for i = 1:length(dataCell)
    % Figure out what the current range is
    endpoint = start + runLengths(i) - 1;
    inds = start:endpoint;
    % Drop the data in
    output(inds, :) = dataCell{i};
    % Increment the start point
    start = endpoint + 1;
end