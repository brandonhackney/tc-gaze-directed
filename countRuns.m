function numRuns = countRuns(subNum, taskName)
% Checks all 'ses-xx' folders in BIDS/sub-xx for a bold.json
% Reports the number that match 'task-yyy'
    assert(isstring(taskName) || ischar(taskName), '2nd input must be a string/char!');
    p = specifyPaths;
    subID = sprintf('sub-%02.f', subNum);
    numSes = height(dir(fullfile(p.bids, subID, 'ses-*')));
    tmp = zeros(numSes, 1);
    for i = 1:numSes
        sesID = sprintf('ses-%02.f', i);
        fdir = fullfile(p.bids, subID, sesID, 'func');
        fspec = strjoin({subID, sesID, ['task-', taskName], 'run-*', 'bold.json'}, '_');
        flist = dir(fullfile(fdir,fspec));
        tmp(i) = height(flist);
    end
    numRuns = sum(tmp);
    if numRuns == 0
        warning('No data found for "task-%s". Verify spelling?', taskName);
    end
end