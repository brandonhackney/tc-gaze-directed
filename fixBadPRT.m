function fixBadPRT(prt)
% The PRT files in source have the wrong data:
% They use onset & duration instead of onset & offset
% This simply adds onset to duration to get offset, then overwrites.

timing = prt.Cond.OnOffsets;
% Make sure you don't double-dip on accident
if timing(2,1) > timing(2,2)
    % If onset is larger than offset,
    % then 'offset' is likely actually duration.
    % So sum the onset and duration to get offset.
    prt.Cond.OnOffsets(:,2) = timing(:,1) + timing(:,2);
end
% No else statement - if it looks good, skip it.
end