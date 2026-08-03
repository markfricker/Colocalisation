function [scoreImage, gatedMask] = colocScoreImage(ch1, ch2, mask, T1, T2, scoreType)
%COLOCSCOREIMAGE  Per-pixel colocalisation score, gated by threshold.
%
%   [scoreImage, gatedMask] = colocScoreImage(ch1, ch2, mask, T1, T2, scoreType)
%
% Produces a per-pixel colocalisation "map" for spatial interpretation of
% where a scalar coefficient (e.g. Manders M1/M2) is actually coming from
% -- a single number can't show whether colocalisation is concentrated in
% one region or spread throughout. Pixels failing the mask or the
% threshold gate (ch1<=T1 or ch2<=T2) are NaN; render these as grey/
% background, not part of the colour-mapped score.
%
% INPUTS
%   ch1, ch2  - equal-size 2D numeric intensity images.
%   mask      - (optional) logical, same size, restricting the pixel
%               population (e.g. cell boundary, with saturated pixels
%               already excluded by the caller). Default: all finite
%               pixels in both channels.
%   T1, T2    - intensity thresholds for ch1/ch2 (e.g. from
%               colocCostesThreshold, or a manual/background fallback).
%               A pixel is gated in only if ch1>T1 AND ch2>T2.
%   scoreType - 'geomean' (default) or 'pdm':
%     'geomean' - sqrt(ch1n .* ch2n), where ch1n/ch2n are each channel
%                 linearly rescaled 0..1 between its own threshold and
%                 its own max over the gated population. Always >=0;
%                 reads like a brightness image of "how strongly both
%                 channels agree here". Matches the Imaris/Zen
%                 "colocalisation channel" display convention.
%     'pdm'     - product of the differences from the mean,
%                 (ch1-mean(ch1(mask))).*(ch2-mean(ch2(mask))), evaluated
%                 at gated pixels for display (the means themselves are
%                 computed over the whole mask, matching the standard
%                 whole-image PDM definition). Signed: positive = this
%                 pixel pulls Pearson's r upward, negative = pulls it
%                 down. Use a diverging colormap centred at 0.
%
% OUTPUTS
%   scoreImage - same size as ch1, NaN outside the gated population.
%   gatedMask  - logical, same size, the pixels scoreImage is defined at.
%
% REFERENCES
%   Bolte S, Cordelieres FP (2006) A guided tour into subcellular
%   colocalization analysis in light microscopy. J Microsc 224:213-232.
%   [PDM / product-of-differences-from-the-mean colocalisation map,
%   popularised by the JACoP ImageJ plugin]

if nargin < 3 || isempty(mask)
    mask = true(size(ch1));
end
if nargin < 6 || isempty(scoreType)
    scoreType = 'geomean';
end

mask = mask & isfinite(ch1) & isfinite(ch2);
gatedMask = mask & ch1 > T1 & ch2 > T2;

scoreImage = nan(size(ch1));

if ~any(gatedMask, 'all')
    return
end

switch lower(scoreType)
    case 'geomean'
        ch1g = double(ch1(gatedMask));
        ch2g = double(ch2(gatedMask));
        range1 = max(ch1g) - T1;
        range2 = max(ch2g) - T2;
        if range1 <= 0 || range2 <= 0
            return
        end
        ch1n = (ch1g - T1) ./ range1;
        ch2n = (ch2g - T2) ./ range2;
        scoreImage(gatedMask) = sqrt(ch1n .* ch2n);

    case 'pdm'
        m1 = mean(double(ch1(mask)));
        m2 = mean(double(ch2(mask)));
        ch1g = double(ch1(gatedMask));
        ch2g = double(ch2(gatedMask));
        scoreImage(gatedMask) = (ch1g - m1) .* (ch2g - m2);

    otherwise
        error('colocScoreImage:unknownScoreType', ...
            'Unknown scoreType ''%s'' -- expected ''geomean'' or ''pdm''.', scoreType);
end

end % colocScoreImage
