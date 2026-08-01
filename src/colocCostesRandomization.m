function [pValue, robs, rNull] = colocCostesRandomization(ch1, ch2, mask, nIterations, blockSize, rngSeed)
%COLOCCOSTESRANDOMIZATION  Costes' randomisation significance test for
% pixel-based colocalisation.
%
%   [pValue, robs, rNull] = colocCostesRandomization(ch1, ch2, mask, ...
%       nIterations, blockSize, rngSeed)
%
% A raw Pearson's r (or Manders' M1/M2) computed on pixel data has no
% valid parametric significance test -- neighbouring pixels are
% spatially autocorrelated (via the imaging point-spread function), so a
% naive IID-based p-value is badly overoptimistic. Costes' randomisation
% instead builds an empirical null distribution: channel 2 is cut into
% square blocks of BLOCKSIZE pixels (chosen to be roughly the imaging
% PSF, so each block's own internal texture/noise structure survives
% intact), the blocks are randomly permuted across the image
% NITERATIONS times, and Pearson's r is recomputed against the
% unshuffled channel 1 each time. The observed r is then compared
% against this null distribution.
%
% Only the largest sub-rectangle that tiles exactly into BLOCKSIZE x
% BLOCKSIZE blocks is shuffled; any leftover border strip narrower than
% BLOCKSIZE is left unshuffled (identical to the original channel 2 in
% every iteration). This avoids padding artefacts and is a negligible
% approximation whenever the image is many blocks wide/tall.
%
% INPUTS
%   ch1, ch2    - equal-size 2D numeric intensity images.
%   mask        - (optional) logical, same size, restricting the pixel
%                 population. Default: all finite pixels in both
%                 channels.
%   nIterations - number of block-randomisations (default 100, matching
%                 Coloc2's default).
%   blockSize   - block edge length in pixels for the spatial scramble
%                 (default 3; should be set to roughly the PSF FWHM in
%                 pixels -- too small destroys the very autocorrelation
%                 structure the test needs to preserve, too large leaves
%                 too few blocks to randomise meaningfully). Must be
%                 <= min(size(ch1)).
%   rngSeed     - (optional) seed for reproducible randomisation. The
%                 global RNG state is saved and restored afterwards.
%
% OUTPUTS
%   pValue - fraction of randomised iterations whose r met or exceeded
%            the observed r, i.e. P(r_random >= r_observed). Small
%            values (e.g. <0.05, or 0/nIterations) support genuine
%            positive colocalisation rather than coincidental pixel
%            overlap.
%   robs   - the observed (unrandomised) Pearson's r, for reference.
%   rNull  - [nIterations x 1] the null distribution of randomised r
%            values.
%
% REFERENCES
%   Costes SV, Daelemans D, Cho EH, Dobbin Z, Pavlakis G, Lockett S (2004)
%   Automatic and quantitative measurement of protein-protein
%   colocalization in live cells. Biophys J 86:3993-4003.

if nargin < 3 || isempty(mask)
    mask = true(size(ch1));
end
if nargin < 4 || isempty(nIterations)
    nIterations = 100;
end
if nargin < 5 || isempty(blockSize)
    blockSize = 3;
end

[nY, nX] = size(ch1);
if blockSize > min(nY, nX)
    error('colocCostesRandomization:blockTooLarge', ...
        'blockSize (%d) must be <= min(size(ch1)) (%d).', blockSize, min(nY,nX));
end

restoreRng = [];
if nargin >= 6 && ~isempty(rngSeed)
    prevState = rng(rngSeed);
    restoreRng = onCleanup(@() rng(prevState)); %#ok<NASGU>
end

mask = mask & isfinite(ch1) & isfinite(ch2);
robs = colocPearson(ch1, ch2, mask);

nBy = floor(nY / blockSize);
nBx = floor(nX / blockSize);
yLim = nBy * blockSize;
xLim = nBx * blockSize;
ch2Core = ch2(1:yLim, 1:xLim);

rNull = nan(nIterations, 1);
for iIter = 1:nIterations
    ch2Shuf = ch2;
    ch2Shuf(1:yLim, 1:xLim) = blockScramble(ch2Core, blockSize, nBy, nBx);
    rNull(iIter) = colocPearson(ch1, ch2Shuf, mask);
end

pValue = mean(rNull >= robs);

end % colocCostesRandomization


% =========================================================================
function scrambled = blockScramble(img, blockSize, nBy, nBx)
%BLOCKSCRAMBLE  Randomly permute blockSize x blockSize tiles of img.

blocks = mat2cell(img, repmat(blockSize, 1, nBy), repmat(blockSize, 1, nBx));
order  = randperm(nBy * nBx);
scrambled = cell2mat(reshape(blocks(order), nBy, nBx));

end % blockScramble
