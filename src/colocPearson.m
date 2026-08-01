function [r, n] = colocPearson(ch1, ch2, mask)
%COLOCPEARSON  Pearson's correlation coefficient between two channels.
%
%   [r, n] = colocPearson(ch1, ch2, mask)
%
% Linear correlation of paired pixel intensities across two channels,
% restricted to MASK. This is the workhorse pixel-based colocalisation
% metric reported by every major tool (Coloc2, JACoP, CellProfiler,
% Imaris) -- see REFERENCES. Registration and background subtraction are
% assumed already done upstream of this function.
%
% Deliberately returns no p-value: neighbouring pixels are spatially
% autocorrelated (via the imaging point-spread function), so a
% naive IID-based significance test on pixel data is badly overoptimistic.
% Use colocCostesRandomization for a valid significance test instead.
%
% INPUTS
%   ch1, ch2 - equal-size 2D numeric intensity images (single/double).
%   mask     - (optional) logical, same size as ch1/ch2, restricting the
%              pixel population (e.g. cell boundary). Default: all
%              finite pixels in both channels.
%
% OUTPUTS
%   r - Pearson's correlation coefficient, -1..1. NaN if fewer than 2
%       pixels are in the masked population, or if either channel has
%       zero variance over that population.
%   n - number of pixels used.
%
% REFERENCES
%   Manders EMM, Stap J, Brakenhoff GJ, van Driel R, Aten JA (1992)
%   Dynamics of three-dimensional replication patterns during the
%   S-phase, analysed by double labelling of DNA and confocal microscopy.
%   J Cell Sci 103:857-862. [first microscopy application of Pearson's r]
%   Adler J, Parmryd I (2010) Quantifying colocalization by correlation:
%   the Pearson correlation coefficient is superior to the Manders'
%   overlap coefficient. Cytometry A 77:733-742.
%   Bolte S, Cordelieres FP (2006) A guided tour into subcellular
%   colocalization analysis in light microscopy. J Microsc 224:213-232.

if nargin < 3 || isempty(mask)
    mask = true(size(ch1));
end
mask = mask & isfinite(ch1) & isfinite(ch2);

x = double(ch1(mask));
y = double(ch2(mask));
n = numel(x);

if n < 2
    r = NaN;
    return
end

xc = x - mean(x);
yc = y - mean(y);
denom = sqrt(sum(xc.^2) * sum(yc.^2));

if denom == 0
    r = NaN;
else
    r = sum(xc .* yc) / denom;
end

end % colocPearson
