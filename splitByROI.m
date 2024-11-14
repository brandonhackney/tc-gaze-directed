function [dataOut, dataLabels] = splitByROI(dataIn, hem)
% Given a whole-brain timeseries matrix, split into freesurfer ROIs.
% Reads in the HCP (or "Glasser") atlas, subsets the data into parcels,
% then returns an average timeseries per parcel.
% (replace each vertex with its parcel's average)

pths = specifyPaths;
assert(any(strcmp(hem, {'lh', 'rh'})), '2nd input must be ''lh'' or ''rh''');
% Load ROI data
fname = [hem, '.HCPMMP1.annot'];
fpath = fullfile(pths.atlas, fname);
[~, labels, ctable] = read_annotation(fpath);
numROIs = ctable.numEntries;
dataOut = zeros(size(dataIn));
for r = 1:numROIs
    % Figure out which vertices to use
    roi = ctable.struct_names{r};
    roiID = ctable.table(r,5);
    vertList = labels == roiID;
    % Get the average timecourse within this parcel
    if sum(vertList) > 0
        avgTS = mean(dataIn(:,vertList), 2);
        dataOut(:,vertList) = repmat(avgTS, 1, sum(vertList));
    end
    dataLabels{r} = roi;
end