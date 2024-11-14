function [dataSmall, atlas] = splitByROI(dataIn, hem)
% Given a whole-brain timeseries matrix, split into freesurfer ROIs.
% Reads in the HCP (or "Glasser") atlas, subsets the data into parcels,
% then returns an average timeseries per parcel.
% Also returns an optional struct with ROI labels and vertex indices.

pths = specifyPaths;
assert(any(strcmp(hem, {'lh', 'rh'})), '2nd input must be ''lh'' or ''rh''');
% Load ROI data
fname = [hem, '.HCPMMP1.annot'];
fpath = fullfile(pths.atlas, fname);
[~, labels, ctable] = read_annotation(fpath);
numROIs = ctable.numEntries;
[numTRs, numVerts] = size(dataIn);
% Initialize outputs
dataOut = zeros([numTRs, numVerts]);
dataSmall = zeros([numTRs, numROIs]);
atlas = struct();
for r = 1:numROIs
    % Figure out which vertices to use
    roi = ctable.struct_names{r};
    roiID = ctable.table(r,5);
    vertList = labels == roiID;
    % Get the average timecourse within this parcel
    % Replace all timeseries in this parcel with the average
    if sum(vertList) > 0
        avgTS = mean(dataIn(:,vertList), 2);
        dataOut(:,vertList) = repmat(avgTS, 1, sum(vertList));
        % Also spit out something with one column per ROI
        dataSmall(:,r) = avgTS;
    end
    atlas(r).Label = roi;
    atlas(r).Verts = vertList; % the first one will be empty.
end