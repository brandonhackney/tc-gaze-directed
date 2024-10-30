function out = findTxt(subNum, runNum)
subID = sprintf('sub-%02.f', subNum);
p = specifyPaths;
fname = fullfile(p.beh, subID, [subID, '_task-tricopa_date-*_run-', num2str(runNum), '.txt']);
% Above uses a wildcard for data - see if you can find a matching file
x = dir(fname);
assert(~isempty(x), 'Could not find file %s', fname);
out = fullfile(x(1).folder, x(1).name);
end