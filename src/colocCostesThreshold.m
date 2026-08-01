function [Tr, Tg, rAtThreshold, converged] = colocCostesThreshold(ch1, ch2, mask)
%COLOCCOSTESTHRESHOLD  Costes' automatic colocalisation threshold.
%
%   [Tr, Tg, rAtThreshold, converged] = colocCostesThreshold(ch1, ch2, mask)
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
% INPUTS
%   ch1, ch2 - equal-size 2D numeric intensity images.
%   mask     - (optional) logical, same size, restricting the pixel
%              population. Default: all finite pixels in both channels.
%
% OUTPUTS
%   Tr, Tg       - channel-1/channel-2 thresholds. NaN if the regression
%                  slope isn't positive (Costes' method assumes the two
%                  channels trend together) or if r never reaches <=0
%                  while lowering the threshold.
%   rAtThreshold - Pearson's r of the below-threshold population at the
%                  chosen threshold (<=0 by construction when converged).
%   converged    - true if a valid threshold was found.
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

if n < 3
    return
end

fit = [x, ones(n,1)] \ y;   % ordinary least squares: y = a*x + b
a = fit(1);
b = fit(2);

if ~(a > 0)
    warning('colocCostesThreshold:nonPositiveSlope', ...
        'Regression slope is not positive; Costes'' threshold is undefined for this channel pair.');
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
if isempty(kThresh)
    warning('colocCostesThreshold:noConvergence', ...
        'Pearson''s r never dropped to <=0 while lowering the threshold; no valid Costes threshold found.');
    return
end

Tr = xs(kThresh);
Tg = a * Tr + b;
rAtThreshold = rPrefix(kThresh);
converged = true;

end % colocCostesThreshold
