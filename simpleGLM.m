function [varargout] = simpleGLM(timeseries,designmatrix,contrast)
% [betas, residuals, R2, tMap] = simpleGLM(timeseries,designmatrix,(contrast))
% Intended to take a whole-brain timeseries and compute a GLM contrast
% timeseries is an n * m matrix of fMRI values at n timepoints for m voxels
% designmatrix is an n * p matrix of n timepoints for p predictor variables
% contrast is a 1 x p vector of weights on each predictor, incl constant
% contrast can be generated by getContrastVector first
%
% betas is a p * m matrix of estimates
% residuals is an n * m matrix of error terms
% R2 is a 1 * m matrix of model fits
% tMap will return a 1 * m vector of t statistics for each voxel

    % Verify designmatrix includes an intercept term
    % Check the front and the back. If neither works, insert in the front.
    if ~all(designmatrix(:,1) == 1) && ~all(designmatrix(:,end) == 1)
        designmatrix = [ones(height(designmatrix), 1), designmatrix];
    end
    
    % Calculate regression via matrix math.
    % Old method was "unstable", as X'*X can become near-singular
%     betas = (designmatrix' * designmatrix) \ designmatrix' * timeseries;
    % Instead, try QR decomposition to improve "stability"
    [Q,R] = qr(designmatrix, 0); % The 0 ensures an economy-sized decomp
    betas = R \ (Q' * timeseries); % Equivalent to R^(-1) * (Q' * Y)
    residuals = timeseries - designmatrix * betas;
    
    varargout = cell(1,nargout);
    for i = 1:length(varargout)
        switch i
            case 1
                varargout{i} = betas;
            case 2
                varargout{i} = residuals;
            case 3
                % Calculate R2
                SS_total = sum((timeseries - mean(timeseries)).^2, 1);
                SS_residual = sum(residuals.^2, 1);
                varargout{i} = 1 - (SS_residual ./ SS_total);
            case 4
                % Calculate tMap
                assert(nargin > 2, 'No contrast vector provided! Cannot calculate t-statistics')
                varargout{i} = (contrast * betas) ./ sqrt((std(residuals).^2) .* (contrast / (designmatrix' * designmatrix) * contrast'));

            otherwise
                warning('Too many output arguments. Options are [tMap,betas,residuals]. Extra outputs ignored.')
                break
        end
    end
end