function p = colocParamsDefault()
%COLOCPARAMSDEFAULT  Factory default parameters for the pixel-based
% colocalisation library (colocPixelBasedRun).
%
%   p = colocParamsDefault()
%
%   colocCostesRandomization
%   ----------------------------------------------------------------------
%   nIterations   Costes randomisation iterations.
%   blockSize     Randomisation block edge length, px -- set to roughly
%                 the imaging PSF FWHM in pixels.
%   rngSeed       Randomisation seed, [] for non-reproducible (each call
%                 reseeds from system entropy).
%
%   colocManders fallback
%   ----------------------------------------------------------------------
%   manualThreshold1/2   Fallback M1/M2 thresholds used only when Costes'
%                        threshold fails to converge (e.g. non-positive
%                        regression slope) -- default 0 (classic "any
%                        non-zero pixel" Manders convention).
%
%   Saturation exclusion
%   ----------------------------------------------------------------------
%   saturationValue1/2   Pixels at or above this value in ch1/ch2 are
%                        excluded from every statistic and the score
%                        image -- clipped intensities distort the true
%                        intensity relationship, which specifically
%                        corrupts Costes' regression fit. Default Inf (no
%                        exclusion): only the user knows their instrument's
%                        real ceiling and typical background level, so
%                        this function itself never auto-detects a value
%                        -- a caller wanting an auto-detect policy (e.g.
%                        "exclude pixels at a channel's own observed max")
%                        computes a concrete value and passes it in here
%                        like any other threshold.
%
%   colocPixelBasedRun always computes both colocScoreImage variants
%   ('geomean' and 'pdm') -- no params field needed to select one.

p.nIterations      = 100;
p.blockSize        = 3;    % px, ~PSF FWHM
p.rngSeed          = [];
p.manualThreshold1 = 0;
p.manualThreshold2 = 0;

p.saturationValue1 = Inf;
p.saturationValue2 = Inf;

end % colocParamsDefault
