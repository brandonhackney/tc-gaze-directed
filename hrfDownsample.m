function output = hrfDownsample(sdm, frameCol, TRvec)
% SDM is a matrix of n timepoints by m parameters.
% frameCol is a vector describing the timing of those values.
% TRvec is another vector describing the timepoints to downsample to.
% Take the SDM, convolve it with an HRF to get an expected BOLD signal,
% then interpolate from the sample rate in frameCol to the rate in TRvec.
% Convolve with the HRF first to preserve as much data as possible,
% relative to downsampling first.

% Get kernel
SR = diff(frameCol); SR = SR(1);
hrf = spm_hrf(SR);
% Init output
output = zeros(length(TRvec), width(sdm));
% Operate over each column individually
for i = 1:width(sdm)
    col = sdm(:,i);
    col = conv(col, hrf, 'full');
    col = col(1:height(sdm)); % chop off the trailing portion
    output(:,i) = interp1(frameCol, col, TRvec);
end

end