function [T1, T2, rAtThreshold, converged, fullyConverged, slope, intercept] = colocCostesThreshold(ch1, ch2, mask)
%COLOCCOSTESTHRESHOLD  Costes' automatic colocalisation threshold.
%
%   [T1, T2, rAtThreshold, converged, fullyConverged, slope, intercept] = ...
%       colocCostesThreshold(ch1, ch2, mask)
%
% Finds the intensity thresholds T1 (channel 1) and T2 (channel 2) that
% objectively separate background/coincidental overlap from genuine
% colocalising signal, removing the subjectivity of a hand-dragged
% Manders' threshold.
%
% Costes' procedure (Costes et al. 2004, "Quantification of colocalization
% and identification of colocalized pixels" + Fig. 1), as implemented in
% Fiji's Coloc2 (AutoThresholdRegression, PearsonsCorrelation):
%   1. Fit the line ch2 = a*ch1 + b over the whole (masked) population by
%      ORTHOGONAL regression (total least squares) -- symmetric in the two
%      channels, unlike ordinary least squares of ch2 on ch1, whose slope
%      is pulled towards 0 when r is low.
%   2. Walk the threshold pair (T, a*T + b) DOWN along that line from the
%      brightest pixel. The working threshold runs on channel 1 when the
%      line is shallower than the diagonal (a < 1) and on channel 2
%      otherwise (Coloc2's rule), so swapping the channels gives the same
%      threshold pair.
%   3. At each step, compute Pearson's r of the pixels BELOW the threshold
%      pair, i.e. every pixel that is not above both thresholds
%      (ch1 < T1 OR ch2 < T2 -- Coloc2's ThresholdMode.Below).
%   4. Stop at the highest threshold where that r first drops to <= 0.
%      Pixels with ch1 >= T1 AND ch2 >= T2 are the colocalised population.
%
% This implementation computes the walk exactly rather than in fixed
% intensity steps: for positive slope a pixel stays above the threshold
% pair for every working threshold t <= s, where s = min(ch1, (ch2-b)/a)
% (walking on ch1) or s = min(a*ch1+b, ch2) (walking on ch2). Sorting the
% pixels once by s turns every possible below-threshold population into a
% prefix (the k pixels with the smallest s), whose r comes from cumulative
% sums in a single O(n log n) pass. The threshold returned is the s of the
% first pixel above it, so the "above" test is ch >= T, matching Coloc2.
%
% Requiring r to actually cross <=0 assumes there is a genuinely
% uncorrelated background population somewhere in the sampled pixels.
% That assumption can fail -- e.g. when the mask already excludes
% true background (a tight cell-boundary ROI), or when two adjacent
% frames of the same live-cell data simply look very similar throughout
% their whole intensity range -- in which case r asymptotes to some
% small positive value and never reaches zero. Rather than discard the
% fit, this falls back to the threshold at MINIMUM r (still read off the
% same regression line/pixel ranking), flagging that result as only
% partially converged via fullyConverged. (Coloc2 instead stops at the
% first step where r rises; with exact per-pixel steps that would fire on
% the first noise-level wobble, so the global minimum is used.)
%
% INPUTS
%   ch1, ch2 - equal-size 2D numeric intensity images.
%   mask     - (optional) logical, same size, restricting the pixel
%              population. Default: all finite pixels in both channels.
%
% OUTPUTS
%   T1, T2       - channel-1/channel-2 thresholds (named for the channel
%                  index, not colour -- AnalyzER's two compared channels
%                  are user-selectable, not necessarily red/green). A
%                  pixel is above threshold when ch1 >= T1 AND ch2 >= T2.
%                  NaN if the regression slope isn't positive (Costes'
%                  method assumes the two channels trend together) or if
%                  fewer than 3 pixels are available.
%   rAtThreshold - Pearson's r of the below-threshold population at the
%                  chosen threshold (<=0 when fullyConverged; the
%                  achieved minimum, still >0, otherwise).
%   converged    - true if a usable threshold was found at all (full or
%                  partial/min-r fallback).
%   fullyConverged - true only if r actually reached <=0 (the strict
%                  Costes criterion); false if the min-r fallback was
%                  used instead, or if converged is false.
%   slope, intercept - the orthogonal regression line ch2 = slope*ch1 +
%                  intercept (NaN if fewer than 3 pixels or no variance).
%
% REFERENCES
%   Costes SV, Daelemans D, Cho EH, Dobbin Z, Pavlakis G, Lockett S (2004)
%   Automatic and quantitative measurement of protein-protein
%   colocalization in live cells. Biophys J 86:3993-4003.
%   Fiji Coloc2: sc.fiji.coloc.algorithms.AutoThresholdRegression,
%   PearsonsCorrelation (github.com/fiji/Colocalisation_Analysis).

if nargin < 3 || isempty(mask)
    mask = true(size(ch1));
end
mask = mask & isfinite(ch1) & isfinite(ch2);

x = double(ch1(mask));
y = double(ch2(mask));
n = numel(x);

T1 = NaN;
T2 = NaN;
rAtThreshold = NaN;
converged = false;
fullyConverged = false;
slope = NaN;
intercept = NaN;

if n < 3
    return
end

% Orthogonal regression (Coloc2's closed form): slope of the major axis of
% the ch1/ch2 scatter.
mx = mean(x);
my = mean(y);
varX = sum((x - mx).^2) / (n - 1);
varY = sum((y - my).^2) / (n - 1);
covXY = sum((x - mx) .* (y - my)) / (n - 1);

if ~(covXY > 0)
    % Not exceptional -- callers already get this via converged=false and
    % the GUI already surfaces it; no console warning (was pure noise
    % once colocalisation started running per-cell, firing once per
    % weakly-correlated cell on every Run). covXY > 0 <=> slope > 0.
    if covXY < 0
        slope = (varY - varX - sqrt((varY - varX)^2 + 4*covXY^2)) / (2*covXY);
        intercept = my - slope*mx;
    end
    return
end
a = (varY - varX + sqrt((varY - varX)^2 + 4*covXY^2)) / (2*covXY);
b = my - a*mx;
slope = a;
intercept = b;

% Working threshold on the channel the line leans towards (Coloc2).
walkOnCh1 = a < 1;
if walkOnCh1
    s = min(x, (y - b) / a);
else
    s = min(a*x + b, y);
end

[ss, order] = sort(s, 'ascend');
% Centre before the cumulative sums to keep them well conditioned.
xs = x(order) - mx;
ys = y(order) - my;

cumX  = cumsum(xs);
cumY  = cumsum(ys);
cumXY = cumsum(xs .* ys);
cumX2 = cumsum(xs .^ 2);
cumY2 = cumsum(ys .^ 2);

k = (1:n)';
num = k .* cumXY - cumX .* cumY;
den = sqrt((k .* cumX2 - cumX.^2) .* (k .* cumY2 - cumY.^2));

rPrefix = num ./ den;   % rPrefix(k) = Pearson's r of the k pixels with the smallest s
rPrefix(k < 3) = NaN;
rPrefix(den <= 0) = NaN;
% A prefix is a real below-threshold population only if the next pixel
% (the first one above) has a strictly larger s; ties can't be split.
% The full population (k = n) has no threshold above it.
valid = [ss(1:end-1) < ss(2:end); false];
rPrefix(~valid) = NaN;

kThresh = find(rPrefix <= 0, 1, 'last');
fullyConverged = ~isempty(kThresh);
if ~fullyConverged
    [rMin, kThresh] = min(rPrefix, [], 'omitnan');
    if isempty(rMin) || isnan(rMin)
        return
    end
end

t = ss(kThresh + 1);
if walkOnCh1
    T1 = t;
    T2 = a*t + b;
else
    T2 = t;
    T1 = (t - b) / a;
end
rAtThreshold = rPrefix(kThresh);
converged = true;

end % colocCostesThreshold
