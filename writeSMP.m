function writeSMP(fullMap, subj, label)
% Find parameters
numVert = length(fullMap);
numTRs = 177; % how to find??
df = numTRs - 1;
LowerThreshold = 30; % R^2 of 0.3 * 100
UpperThreshold = 100; % R^2 of 1 * 100
hem = 'rh'; % how to determine??
% [cluster,LowerThreshold] = fdrCluster(fullMap,df); % expects T values
pths = specifyPaths;
% dataDir = pths.bv;
dataDir = '/data/home/brandon/Documents/surf/'; % test
% Write to variable
SMP = xff('new:smp');
SMP.NrOfVertices = numVert;
SMP.NameOfOriginalSRF = fullfile(dataDir, subj, [subj, '-Surf2BV'], [subj '_' hem '_smoothwm.srf']);
% SMP.NameOfOriginalSRF = [dataDir filesep subj filesep subj '-Freesurfer' filesep subj '-Surf2BV' filesep subj '_' hem '_smoothwm.srf'];
SMP.Map.Type = 16; % 1 = t-stat, 2 = correlation, 4 = F-stat, 5 = Z-stat, 15 = beta, 16 = probability, there are others.
SMP.Map.SMPData = fullMap;
SMP.Map.Name = label;
SMP.Map.DF1 = df;
SMP.Map.BonferroniValue = numVert;
% SMP.Map.LowerThreshold = LowerThreshold;
% SMP.Map.UpperThreshold = UpperThreshold;

% Write to file
% SMP.SaveAs(fullfile(dataDir, subj, [subj '_' atlasName '_' contName '_' hem '.smp']));
SMP.SaveAs(fullfile(dataDir, subj, [label '.smp']));
SMP.clearobj;