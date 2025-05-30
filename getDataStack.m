function [dataStack,  roiLabels] = getDataStack(subNum, hem, varargin)
% Extract data for this subject, subset to ROIs, and get design matrix

if nargin > 2
    dspace = varargin{1};
else
    dspace = [];
end

numRuns = countRuns(subNum, 'tricopa');   
dataStack = cell(numRuns, 1);

hem1 = {'L', 'R'}; % different functions want different formats
hem2 = {'lh', 'rh'};

% Now go:
fprintf(1, 'Subject %i: ', subNum);
fprintf(1, 'Getting data for %i runs\n', numRuns);
tic;
for r = 1:numRuns
      % % PREDICTORS:
        % Need to account for nuisance regressors, like head motion
        % Except the full files are not the same width between runs!!
        % So only grab a select few predictors
        % And keep these separate, since we will predict in steps
        nuisance = loadNuisance(subNum, r);
        
        % z-score your predictors
        nuisance = zscore(nuisance);
        
        % Insert an intercept column
        nuisance = [ones(height(nuisance), 1), nuisance];

      % % DATA:
        % Load the actual MRI data for all runs
        dataStack{r} = loadData(subNum, r, hem1{hem}, dspace);
        
%%%%%%%%% TEST
        % Replace the real data with random data
%         dataStack{r} = rand(size(dataStack{r})) .* mean(dataStack{r}, 'all');
%%%%%%%%% TEST
        

        % PREPROCESS DATA:
        % Immediately z-score the data, to remove any run-specific effect
%         dataStack{r} = zscore(dataStack{r});
        % Regress the nuisance out of the data, like head motion
        [~, dataStack{r}, r2(r,:)] = simpleGLM(dataStack{r}, nuisance);
        % Outlier rejection: 
        % Regress a stimulus-timing predictor (unused in main analysis)
        % Keep the top 10% of vertices in each parcel, based on the R^2
        % This helps reduce noise using a data-driven approach,
        % while avoiding circularity issues by using a unique predictor.
        
        % feed this into splitByROI to use as your test values
        [timing, ~] = getSDM(subNum, r, {'Timing'});
        [~,~,testr2(1,:)] = simpleGLM(dataStack{r}, timing);

        % Then avg within ROIs to further reduce computational load
%         [dataStack{r}, roiLabels] = splitByROI(dataStack{r}, hem2{hem});
        [dataStack{r}, roiLabels] = splitByROI(dataStack{r}, hem2{hem}, testr2);
        
        % z-score those residuals independently within each vertex?
%         dataStack{r} = zscore(dataStack{r});
        % OR z-score those residuals relative to whole-brain variance??
%         dataStack{r} = zscore(dataStack{r}, [], 'all');

end % for each run
fprintf(1, 'Done processing data for subject %i.\n', subNum);
toc

% Report the nuisance regression R2s for each run from this subject
% plotR2s(r2, subNum);

end % function