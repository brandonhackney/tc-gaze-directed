function makeDummySSM(sub, hem)
% BV requires a coregistration file of some sort to do group analyses,
% but the data for this project is already in fsaverage space,
% so we don't need to further coregister everything.
% Trick BV into cooperating by generating a dummy coregistration file
% that maps each subject's fsaverage 1:1 with the 'group' fsaverage

% Generate filename
pths = specifyPaths;
fname = strjoin({sub, hem, 'smoothwm', 'fsaverage.ssm'}, '_');
fpath = pths.bv; % ...?
fout = fullfile(fpath, fname);

% Only run this script if the file doesn't yet exist
if ~exist(fout, 'file')
    
    % Determine the number of vertices used
    % Load in a related gifti and measure it
    ftest = findSubData(pths, sub, 1, ['hemi-', hem, '_space-fsaverage_bold.func.gii']);
    data = gifti(ftest);
    numVerts = size(data.cdata, 1);
    clear data;
    
    % Now do the thing we set out to do
    ssm = xff('new:ssm'); % init an empty SSM

    % Write in the data
    ssm.NrOfTargetVertices = numVerts;
    ssm.NrOfSourceVertices = numVerts;
    ssm.SourceOfTarget = [1:numVerts]'; % make vertical

    ssm.SaveAs(fout);

    % XFF requires manual garbage collection
    ssm.clearobj;
    xff(0, 'unwindstack');
end