function numTRs = findRunLength(subNum, runNum)
% Get the number of TRs for a specific scan based on subID and run num
% Read from a specific fMRIPrep output instead of a larger .nii
    % First make sure the file we want actually exists
    subID = sprintf('sub-%02.f', subNum);
    p = specifyPaths;
    testName = fullfile(p.derivatives, subID, 'ses-02', 'func', [subID, '_ses-02_task-tricopa_run-', num2str(runNum), '_desc-confounds_timeseries.tsv']);
    assert(exist(testName, 'file'), 'Cannot find file %s', testName);
    % Then if the data does exist, we just want the # of rows in that file
    data = readtable(testName, 'FileType', 'text', 'Delimiter', '\t');
    numTRs = height(data);
end