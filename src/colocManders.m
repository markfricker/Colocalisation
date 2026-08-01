function [M1, M2] = colocManders(ch1, ch2, mask, T1, T2)
%COLOCMANDERS  Manders' overlap coefficients M1/M2 for two channels.
%
%   [M1, M2] = colocManders(ch1, ch2, mask, T1, T2)
%
% M1 is the fraction of channel 1's above-threshold intensity that
% occurs at pixels where channel 2 is also above its own threshold, and
% M2 is the symmetric quantity for channel 2. Unlike Pearson's r, M1/M2
% distinguish which channel drives an asymmetric overlap and tolerate
% very different absolute intensities/gains between channels -- but the
% result depends entirely on the choice of T1/T2, so use
% colocCostesThreshold rather than a hand-picked cutoff wherever
% possible.
%
% INPUTS
%   ch1, ch2 - equal-size 2D numeric intensity images.
%   mask     - (optional) logical, same size, restricting the pixel
%              population (e.g. cell boundary). Default: all finite
%              pixels in both channels.
%   T1, T2   - intensity thresholds for ch1/ch2 respectively. A pixel
%              counts as "above" when its value is strictly greater than
%              the threshold. Typically the Costes automatic thresholds
%              from colocCostesThreshold; default 0 (classic "any
%              non-zero pixel" Manders coefficient) if omitted.
%
% OUTPUTS
%   M1 - fraction of ch1's above-threshold intensity colocalising with
%        above-threshold ch2, 0..1. NaN if ch1 has no above-threshold
%        signal at all.
%   M2 - fraction of ch2's above-threshold intensity colocalising with
%        above-threshold ch1, 0..1. NaN if ch2 has no above-threshold
%        signal at all.
%
% REFERENCES
%   Manders EMM, Verbeek FJ, Aten JA (1993) Measurement of co-localization
%   of objects in dual-colour confocal images. J Microsc 169:375-382.
%   Costes SV, Daelemans D, Cho EH, Dobbin Z, Pavlakis G, Lockett S (2004)
%   Automatic and quantitative measurement of protein-protein
%   colocalization in live cells. Biophys J 86:3993-4003.

if nargin < 3 || isempty(mask)
    mask = true(size(ch1));
end
if nargin < 4 || isempty(T1)
    T1 = 0;
end
if nargin < 5 || isempty(T2)
    T2 = 0;
end

mask = mask & isfinite(ch1) & isfinite(ch2);

x = double(ch1(mask));
y = double(ch2(mask));

above1 = x > T1;
above2 = y > T2;
both   = above1 & above2;

sumX = sum(x(above1));
if sumX == 0
    M1 = NaN;
else
    M1 = sum(x(both)) / sumX;
end

sumY = sum(y(above2));
if sumY == 0
    M2 = NaN;
else
    M2 = sum(y(both)) / sumY;
end

end % colocManders
