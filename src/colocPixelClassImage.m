function classImage = colocPixelClassImage(ch1, ch2, mask, T1, T2)
%COLOCPIXELCLASSIMAGE  Label every analysed pixel by which channels are
% above their threshold: neither / channel 1 only / channel 2 only / both.
%
%   classImage = colocPixelClassImage(ch1, ch2, mask, T1, T2)
%
% The score images (colocScoreImage) are defined only on the colocalised
% population, so pixels where just ONE channel is present -- e.g. ER that
% does not overlap any mitochondrion -- are blank there although they
% count fully in Pearson's r and in the Manders denominators. This map
% shows every analysed pixel, using the same test as colocManders and
% colocScoreImage (value > 0 and >= its threshold), so its classes match
% the M1/M2 numbers:
%
%   class 0  not analysed (outside the mask, non-finite, or saturated)
%   class 1  neither channel above its threshold
%   class 2  channel 1 only
%   class 3  channel 2 only
%   class 4  both
%
%   M1 = sum(ch1 in classes 3,4) / sum(ch1 in classes 1-4)
%   M2 = sum(ch2 in classes 2,4) / sum(ch2 in classes 1-4)
%
% (M1 counts channel 1 wherever channel 2 passes its threshold; channel 1
% itself may be below T1 there, so class 3 is included.)
%
% INPUTS
%   ch1, ch2 - equal-size 2D numeric intensity images.
%   mask     - (optional) logical, same size, the analysed population
%              (cell region, saturated pixels already removed). Default:
%              all finite pixels in both channels.
%   T1, T2   - intensity thresholds (the pair Manders used).
%
% OUTPUT
%   classImage - uint8, same size as ch1, values 0..4 as above. Colours
%                and legend labels: colocPixelClassColours.

if nargin < 3 || isempty(mask)
    mask = true(size(ch1));
end

mask = mask & isfinite(ch1) & isfinite(ch2);
in1 = ch1 > 0 & ch1 >= T1;
in2 = ch2 > 0 & ch2 >= T2;

classImage = zeros(size(ch1), 'uint8');
classImage(mask) = 1;
classImage(mask &  in1 & ~in2) = 2;
classImage(mask & ~in1 &  in2) = 3;
classImage(mask &  in1 &  in2) = 4;

end % colocPixelClassImage
