function output = expandROIs(data, roiStruct)
% Given an i*j matrix of i timepoints by j ROIs,
% expand to an i*k matrix with k vertices.
% Requires a list of vertices per parcel, which we take as an input.
% Expected format is a struct with j rows and 2 fields: Label and Verts

% Parse inputs
dw = width(data);
numROIs = height(roiStruct);
assert(width(data) == length(roiStruct), 'Unequal number of ROIs! Input 1 is %i columns wide, but Input 2 is %i rows tall.', dw, numROIs)

% Init output
numTRs = height(data);
numVerts = numel(roiStruct(1).Verts);
output = zeros([numTRs, numVerts]);

for r = 1:numROIs
    avgTS = data(:,r);
    vertList = roiStruct(r).Verts;
    if sum(vertList) > 0
        output(:,vertList) = repmat(avgTS, 1, sum(vertList));
    end
end