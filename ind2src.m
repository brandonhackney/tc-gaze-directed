function srcNum = ind2src(i)
% The subject IDs in source do not match the subject IDs in derivatives
% e.g. if you have source [1 2 4 6], derivatives changes it to [1 2 3 4]
% This function returns the source subject number based on the deriv num.

% Scan the source folder for a subject list
pths = specifyPaths;
dlist = dir(fullfile(pths.beh, 'sub-*'));

% Convert text to numbers
slist = {dlist.name};
slist = str2num(cell2mat(erase(slist, 'sub-')')); %#ok<ST2NM>

% Index subID from the source list based on the derivatives subID
srcNum = slist(i);