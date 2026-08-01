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

p.nIterations      = 100;
p.blockSize        = 3;    % px, ~PSF FWHM
p.rngSeed          = [];
p.manualThreshold1 = 0;
p.manualThreshold2 = 0;

end % colocParamsDefault
