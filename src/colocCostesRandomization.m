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
% intact), the blocks are randomly permuted within the mask
% NITERATIONS times, and Pearson's r is recomputed against the
% unshuffled channel 1 each time. The observed r is then compared
% against this null distribution.
%
% Shuffling stays WITHIN the mask (Costes 2004 scrambles within the region
% being analysed): the BLOCKSIZE x BLOCKSIZE tile grid is anchored at the
% mask's bounding box, tiles lying wholly inside the mask are permuted
% among each other, and the remaining in-mask pixels (the mask's ragged
% edge) are permuted individually among themselves. Nothing from outside
% the mask is ever swapped in -- that would dilute the region with
% background and make p too small. pValue (and rNull) are NaN when fewer
% than two tiles fit inside the mask.
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

rNull = nan(nIterations, 1);
pValue = NaN;
if ~any(mask, 'all')
    return
end

% Tile grid anchored at the mask's bounding box, so a rectangular mask
% tiles exactly and an irregular one loses as few edge pixels as possible.
[maskRows, maskCols] = find(mask);
r0 = min(maskRows);
c0 = min(maskCols);
nBy = floor((max(maskRows) - r0 + 1) / blockSize);
nBx = floor((max(maskCols) - c0 + 1) / blockSize);

insideTiles = [];
pixIdx = zeros(blockSize^2, 0);
if nBy > 0 && nBx > 0
    yLim = nBy * blockSize;
    xLim = nBx * blockSize;
    [colIdx, rowIdx] = meshgrid(c0:c0+xLim-1, r0:r0+yLim-1);
    tileIdx = (ceil((colIdx - c0 + 1) / blockSize) - 1) * nBy + ceil((rowIdx - r0 + 1) / blockSize);
    linCore = sub2ind([nY nX], rowIdx(:), colIdx(:));

    % Only tiles lying wholly inside the mask are block-shuffled, and only
    % among each other -- Costes scrambles within the region being
    % analysed, so background from outside the mask (e.g. outside the
    % cell) must not be swapped in, which would pull the null r towards 0
    % and make p too small.
    tileInside = accumarray(tileIdx(:), double(mask(linCore)), [nBy*nBx 1]) == blockSize^2;
    insideTiles = find(tileInside);

    % pixIdx(:, j) = linear indices (into ch2) of tile insideTiles(j)'s
    % pixels, in the same within-tile order for every tile (stable sort).
    [~, sortOrder] = sort(tileIdx(:));
    pixByTile = reshape(linCore(sortOrder), blockSize^2, nBy*nBx);
    pixIdx = pixByTile(:, insideTiles);
end

if numel(insideTiles) < 2
    % Nothing to shuffle (mask smaller than two blocks) -- no null.
    return
end

% In-mask pixels not covered by a whole inside tile (the mask's ragged
% edge) are shuffled individually among themselves: left in place they
% would keep their true ch1/ch2 pairing in every iteration and bias the
% null r upwards.
edgeMask = mask;
edgeMask(pixIdx) = false;
edgeIdx = find(edgeMask);

for iIter = 1:nIterations
    ch2Shuf = ch2;
    order = randperm(numel(insideTiles));
    ch2Shuf(pixIdx) = ch2(pixIdx(:, order));
    ch2Shuf(edgeIdx) = ch2(edgeIdx(randperm(numel(edgeIdx))));
    rNull(iIter) = colocPearson(ch1, ch2Shuf, mask);
end

pValue = mean(rNull >= robs);

end % colocCostesRandomization
