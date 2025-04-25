function [varargout] = simpleGLM(timeseries,designmatrix,lambda, contrast)
% [betas, residuals, R2] = simpleGLM(timeseries,designmatrix,(lambda))
% Intended to take a whole-brain timeseries and compute a GLM contrast
% INPUTS:
% timeseries is an n * m matrix of fMRI values at n timepoints for m voxels
% designmatrix is an n * p matrix of n timepoints for p predictor variables
% lambda is a scalar weight between 0 and 1 used for ridge regression
% lambda == 0 is the same as normal OLS regression (it cancels the penalty)
% lambda == 1 is the most aggressive regularization possible
% 
% OUTPUTS:
% betas is a p * m matrix of estimates
% residuals is an n * m matrix of error terms
% R2 is a 1 * m matrix of model fits

    % Verify designmatrix includes an intercept term
    % Check the front and the back. If neither works, insert in the front.
    if ~all(designmatrix(:,1) == 1) && ~all(designmatrix(:,end) == 1)
        designmatrix = [ones(height(designmatrix), 1), designmatrix];
    end
    
    if nargin < 3
        % Cancel out the ridge regression unless specifically asked for.
        % Lambda must be a value between 0 and 1, where 1 is aggressive.
        lambda = 0;
    end
    
    % Calculate regression via matrix math.
    % Old method was "unstable", as X'*X can become near-singular
%     betas = (designmatrix' * designmatrix) \ designmatrix' * timeseries;
    % Instead, try QR decomposition to improve "stability"
    [Q,R] = qr(designmatrix, 0); % The 0 ensures an economy-sized decomp
    % Add L2 regularization to R (R + lambda * identity matrix)
    % Helps prevent overfitting, e.g. if removing a var IMPROVES model fit
    d = size(R, 2); % Number of predictors (columns in design matrix)
    R_reg = R + lambda * eye(d); % Add lambda * I to the diagonal of R
    % Calculate regression
    betas = R_reg \ (Q' * timeseries); % Equivalent to R^(-1) * (Q' * Y)
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
                if nargin < 4
                    contrast = ones(1,width(designmatrix));
                end
                varargout{i} = (contrast * betas) ./ sqrt((std(residuals).^2) .* (contrast / (designmatrix' * designmatrix) * contrast'));
            otherwise
                warning('Too many output arguments. Options are [betas,residuals, R2]. Extra outputs ignored.')
                break
        end
    end
end