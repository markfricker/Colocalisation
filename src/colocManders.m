function [M1, M2] = colocManders(ch1, ch2, mask, T1, T2)
%COLOCMANDERS  Manders' split coefficients M1/M2 for two channels.
%
%   [M1, M2] = colocManders(ch1, ch2, mask, T1, T2)
%
% M1 is the fraction of channel 1's TOTAL intensity (all masked pixels)
% that lies at pixels where channel 2 is present (ch2 > 0) and at or
% above its threshold T2; M2 is the same for channel 2 against T1:
%
%   M1 = sum(ch1 where ch2 > 0 & ch2 >= T2) / sum(ch1)
%   M2 = sum(ch2 where ch1 > 0 & ch1 >= T1) / sum(ch2)
%
% This is Manders et al. (1993)'s definition with the "other channel
% present" test raised from > 0 to the threshold, exactly as Fiji's Coloc2
% computes its thresholded tM1/tM2 (MandersColocalization, ThresholdMode
% Above), so values are directly comparable with Coloc2. With T1 = T2 = 0
% it is the classic unthresholded M1/M2. Unlike Pearson's r, M1/M2
% distinguish which channel drives an asymmetric overlap and tolerate
% very different absolute intensities/gains between channels -- but the
% result depends entirely on the choice of T1/T2, so use
% colocCostesThreshold rather than a hand-picked cutoff wherever
% possible. Intensities are assumed background-subtracted (Costes 2004).
%
% INPUTS
%   ch1, ch2 - equal-size 2D numeric intensity images.
%   mask     - (optional) logical, same size, restricting the pixel
%              population (e.g. cell boundary). Default: all finite
%              pixels in both channels.
%   T1, T2   - intensity thresholds for ch1/ch2 respectively (a pixel
%              passes when its value is >= the threshold, as in Coloc2).
%              Typically the Costes automatic thresholds from
%              colocCostesThreshold; default 0 (classic "any non-zero
%              pixel" Manders coefficient) if omitted.
%
% OUTPUTS
%   M1 - fraction of ch1's total intensity colocalising with
%        above-threshold ch2, 0..1. NaN if ch1's total intensity is 0.
%   M2 - fraction of ch2's total intensity colocalising with
%        above-threshold ch1, 0..1. NaN if ch2's total intensity is 0.
%
% REFERENCES
%   Manders EMM, Verbeek FJ, Aten JA (1993) Measurement of co-localization
%   of objects in dual-colour confocal images. J Microsc 169:375-382.
%   Costes SV, Daelemans D, Cho EH, Dobbin Z, Pavlakis G, Lockett S (2004)
%   Automatic and quantitative measurement of protein-protein
%   colocalization in live cells. Biophys J 86:3993-4003.
%   Fiji Coloc2: sc.fiji.coloc.algorithms.MandersColocalization
%   (github.com/fiji/Colocalisation_Analysis).

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

with2 = y > 0 & y >= T2;   % ch2 present and above its threshold
with1 = x > 0 & x >= T1;   % ch1 present and above its threshold

sumX = sum(x);
if sumX == 0
    M1 = NaN;
else
    M1 = sum(x(with2)) / sumX;
end

sumY = sum(y);
if sumY == 0
    M2 = NaN;
else
    M2 = sum(y(with1)) / sumY;
end

end % colocManders
