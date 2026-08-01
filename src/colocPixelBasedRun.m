function results = colocPixelBasedRun(ch1, ch2, mask, p)
%COLOCPIXELBASEDRUN  Full pixel-based colocalisation pipeline for one
% channel-1/channel-2 image pair.
%
%   results = colocPixelBasedRun(ch1, ch2, mask, p)
%
% Runs the standard "Coloc2-equivalent" pixel-based pipeline: Pearson's
% r, Costes' automatic threshold, Manders' M1/M2 at that threshold, and
% a Costes randomisation significance test on the whole-image r.
% Registration and background subtraction are assumed already handled
% upstream of this function (e.g. by AnalyzER's process module).
%
% INPUTS
%   ch1, ch2 - equal-size 2D numeric intensity images.
%   mask     - (optional) logical, same size, restricting the pixel
%              population (e.g. cell boundary). Default: all finite
%              pixels in both channels.
%   p        - (optional) parameter struct, see colocParamsDefault. Any
%              missing field falls back to its factory default.
%
% OUTPUTS
%   results - struct:
%     .pearsonR          - whole-population Pearson's r.
%     .pearsonN          - pixel count used for pearsonR.
%     .costesTr/.costesTg - Costes threshold pair (NaN if not converged).
%     .costesConverged   - logical.
%     .costesR           - r of the below-threshold population (<=0 when
%                           converged).
%     .manders1/.manders2 - M1/M2 at the Costes threshold, or at
%              p.manualThreshold1/2 if Costes' threshold didn't converge.
%     .randPValue        - Costes randomisation p-value on the
%              whole-population r (see colocCostesRandomization).
%     .randNullR         - [p.nIterations x 1] null distribution.
%     .pixelValues       - [n x 2] masked [ch1 ch2] intensity pairs, kept
%              for the planned ch1/ch2 scatter-plot GUI so callers don't
%              have to recompute the mask.
%     .maskLinearIdx     - [n x 1] linear indices into ch1/ch2 for each
%              row of pixelValues, for image<->plot brushing.
%
% REFERENCES
%   See colocPearson, colocManders, colocCostesThreshold,
%   colocCostesRandomization for the individual method references.

if nargin < 3 || isempty(mask)
    mask = true(size(ch1));
end
if nargin < 4 || isempty(p)
    p = colocParamsDefault();
end

mask = mask & isfinite(ch1) & isfinite(ch2);

results = struct();

[results.pearsonR, results.pearsonN] = colocPearson(ch1, ch2, mask);

[results.costesTr, results.costesTg, results.costesR, results.costesConverged] = ...
    colocCostesThreshold(ch1, ch2, mask);

if results.costesConverged
    T1 = results.costesTr;
    T2 = results.costesTg;
else
    T1 = p.manualThreshold1;
    T2 = p.manualThreshold2;
end
[results.manders1, results.manders2] = colocManders(ch1, ch2, mask, T1, T2);

[results.randPValue, ~, results.randNullR] = colocCostesRandomization( ...
    ch1, ch2, mask, p.nIterations, p.blockSize, p.rngSeed);

linIdx = find(mask);
results.pixelValues   = [double(ch1(linIdx)), double(ch2(linIdx))];
results.maskLinearIdx = linIdx;

end % colocPixelBasedRun
