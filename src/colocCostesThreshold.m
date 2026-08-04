function [Tr, Tg, rAtThreshold, converged, fullyConverged] = colocCostesThreshold(ch1, ch2, mask)
%COLOCCOSTESTHRESHOLD  Costes' automatic colocalisation threshold.
%
%   [Tr, Tg, rAtThreshold, converged, fullyConverged] = colocCostesThreshold(ch1, ch2, mask)
%
% Finds the intensity thresholds Tr (channel 1) and Tg (channel 2) that
% objectively separate background/coincidental overlap from genuine
% colocalising signal, removing the subjectivity of a hand-dragged
% Manders' threshold.
%
% Costes' original procedure: fit the least-squares regression line
% ch2 = a*ch1 + b over the whole (masked) population, then walk the
% threshold DOWN from the brightest pixel along that line, recomputing
% Pearson's r for the population still BELOW the current threshold pair
% at each step, and stopping at the highest threshold where that
% below-threshold population's r first drops to <=0 (i.e. indistinguishable
% from, or below, an uncorrelated/background population).
%
% This implementation computes the identical result exactly, but in
% O(n log n) rather than iterating literal threshold steps: pixels are
% sorted once by channel-1 intensity, and r of every possible
% below-threshold prefix is obtained from cumulative sums in a single
% pass (threshold Tr = xs(k) <-> below-threshold population = the k
% dimmest pixels).
%
% Requiring r to actually cross <=0 assumes there is a genuinely
% uncorrelated background population somewhere in the sampled pixels.
% That assumption can fail -- e.g. when the mask already excludes
% true background (a tight cell-boundary ROI), or when two adjacent
% frames of the same live-cell data simply look very similar throughout
% their whole intensity range -- in which case r asymptotes to some
% small positive value and never reaches zero. Rather than discard the
% fit, this falls back to the threshold at MINIMUM r (still read off the
% same regression line/pixel ranking Costes' method already computed),
% flagging that result as only partially converged via fullyConverged.
%
% INPUTS
%   ch1, ch2 - equal-size 2D numeric intensity images.
%   mask     - (optional) logical, same size, restricting the pixel
%              population. Default: all finite pixels in both channels.
%
% OUTPUTS
%   Tr, Tg       - channel-1/channel-2 thresholds. NaN if the regression
%                  slope isn't positive (Costes' method assumes the two
%                  channels trend together) or if fewer than 3 pixels are
%                  available.
%   rAtThreshold - Pearson's r of the below-threshold population at the
%                  chosen threshold (<=0 when fullyConverged; the
%                  achieved minimum, still >0, otherwise).
%   converged    - true if a usable threshold was found at all (full or
%                  partial/min-r fallback).
%   fullyConverged - true only if r actually reached <=0 (the strict
%                  Costes criterion); false if the min-r fallback was
%                  used instead, or if converged is false.
%
% REFERENCES
%   Costes SV, Daelemans D, Cho EH, Dobbin Z, Pavlakis G, Lockett S (2004)
%   Automatic and quantitative measurement of protein-protein
%   colocalization in live cells. Biophys J 86:3993-4003.

if nargin < 3 || isempty(mask)
    mask = true(size(ch1));
end
mask = mask & isfinite(ch1) & isfinite(ch2);

x = double(ch1(mask));
y = double(ch2(mask));
n = numel(x);

Tr = NaN;
Tg = NaN;
rAtThreshold = NaN;
converged = false;
fullyConverged = false;

if n < 3
    return
end

fit = [x, ones(n,1)] \ y;   % ordinary least squares: y = a*x + b
a = fit(1);
b = fit(2);

if ~(a > 0)
    % Not exceptional -- callers already get this via converged=false and
    % the GUI already surfaces it; no console warning (was pure noise
    % once colocalisation started running per-cell, firing once per
    % weakly-correlated cell on every Run).
    return
end

[xs, order] = sort(x, 'ascend');
ys = y(order);

cumX  = cumsum(xs);
cumY  = cumsum(ys);
cumXY = cumsum(xs .* ys);
cumX2 = cumsum(xs .^ 2);
cumY2 = cumsum(ys .^ 2);

k = (1:n)';
num = k .* cumXY - cumX .* cumY;
den = sqrt((k .* cumX2 - cumX.^2) .* (k .* cumY2 - cumY.^2));

rPrefix = num ./ den;   % rPrefix(k) = Pearson's r of the k dimmest (by ch1) pixels
rPrefix(k < 2) = NaN;

kThresh = find(rPrefix <= 0, 1, 'last');
if ~isempty(kThresh)
    Tr = xs(kThresh);
    Tg = a * Tr + b;
    rAtThreshold = rPrefix(kThresh);
    converged = true;
    fullyConverged = true;
    return
end

[rMin, kMin] = min(rPrefix, [], 'omitnan');
if isempty(rMin) || isnan(rMin)
    return
end

Tr = xs(kMin);
Tg = a * Tr + b;
rAtThreshold = rMin;
converged = true;
fullyConverged = false;

end % colocCostesThreshold
