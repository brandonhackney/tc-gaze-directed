function [dataStack,  roiLabels, fullPred, predList] = getDataStack(subNum, hem)
% Extract data for this subject, subset to ROIs, and get design matrix

numRuns = countRuns(subNum, 'tricopa');
fullPred = [];    
dataStack = cell(numRuns, 1);

hem1 = {'L', 'R'}; % different functions want different formats
hem2 = {'lh', 'rh'};

% Now go:
fprintf(1, 'Subject %i: ', subNum);
fprintf(1, 'Getting data for %i runs\n', numRuns);
tic;
for r = 1:numRuns
      % % PREDICTORS:
        % Get the design matrix for each run
        [tmpPred, predList] = getSDM(subNum, r);
%         numPredictors = width(tmpPred) + 1;
        numPredictors = length(predList);
        % Also need to account for nuisance regressors, like head motion
        % Except the full files are not the same width between runs!!
        % So only grab a select few predictors
        % And keep these separate, since we will predict in steps
        nuisance = loadNuisance(subNum, r);
%         tmpPred = [tmpPred, nuisance, timing];
%         tmpPred = [tmpPred, timing, nuisance];
        
        % z-score your predictors
        tmpPred = zscore(tmpPred);
        nuisance = zscore(nuisance);
        
        % Insert an intercept column for both
        tmpPred = [ones(height(tmpPred), 1), tmpPred];
        nuisance = [ones(height(nuisance), 1), nuisance];
        

        % Insert run num as a trailing column, for later subsetting
        tmpPred(:, end+1) = r;
        fullPred = [fullPred; tmpPred];

        
      % % DATA:
        % Load the actual MRI data for all runs
        dataStack{r} = loadData(subNum, r, hem1{hem});
        
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