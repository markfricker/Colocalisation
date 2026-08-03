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
%     .mandersThreshold1/.mandersThreshold2 - the threshold pair actually
%              used for manders1/manders2 above, whichever source it came
%              from -- callers should not assume that's costesTr/costesTg,
%              since those are NaN whenever costesConverged is false.
%     .randPValue        - Costes randomisation p-value on the
%              whole-population r (see colocCostesRandomization).
%     .randNullR         - [p.nIterations x 1] null distribution.
%     .pixelValues       - [n x 2] masked [ch1 ch2] intensity pairs, kept
%              for the planned ch1/ch2 scatter-plot GUI so callers don't
%              have to recompute the mask.
%     .maskLinearIdx     - [n x 1] linear indices into ch1/ch2 for each
%              row of pixelValues, for image<->plot brushing.
%     .scoreImageGeomean/.scoreImagePdm - both per-pixel colocalisation
%              score images (see colocScoreImage), NaN outside the
%              mandersThreshold1/2 gate -- same threshold Manders itself
%              used, whichever source it came from, so the maps and the
%              M1/M2 numbers are always consistent with each other. Both
%              are always computed (cheap relative to the randomisation
%              test) so a caller can offer either without re-running.
%
% REFERENCES
%   See colocPearson, colocManders, colocCostesThreshold,
%   colocCostesRandomization, colocScoreImage for the individual method
%   references.

if nargin < 3 || isempty(mask)
    mask = true(size(ch1));
end
if nargin < 4 || isempty(p)
    p = colocParamsDefault();
end

mask = mask & isfinite(ch1) & isfinite(ch2);

% Exclude saturated pixels -- clipped intensities distort the true
% intensity relationship and specifically corrupt Costes' regression fit.
% Default Inf (no exclusion); see colocParamsDefault for why this is
% never auto-detected.
if isfield(p, 'saturationValue1') && isfinite(p.saturationValue1)
    mask = mask & ch1 < p.saturationValue1;
end
if isfield(p, 'saturationValue2') && isfinite(p.saturationValue2)
    mask = mask & ch2 < p.saturationValue2;
end

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
results.mandersThreshold1 = T1;
results.mandersThreshold2 = T2;

[results.randPValue, ~, results.randNullR] = colocCostesRandomization( ...
    ch1, ch2, mask, p.nIterations, p.blockSize, p.rngSeed);

linIdx = find(mask);
results.pixelValues   = [double(ch1(linIdx)), double(ch2(linIdx))];
results.maskLinearIdx = linIdx;

results.scoreImageGeomean = colocScoreImage(ch1, ch2, mask, T1, T2, 'geomean');
results.scoreImagePdm     = colocScoreImage(ch1, ch2, mask, T1, T2, 'pdm');

end % colocPixelBasedRun
