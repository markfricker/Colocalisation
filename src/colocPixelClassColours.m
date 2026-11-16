function [cmap, labels] = colocPixelClassColours(name1, name2)
%COLOCPIXELCLASSCOLOURS  Palette + legend labels for colocPixelClassImage.
%
%   [cmap, labels] = colocPixelClassColours(name1, name2)
%
% Single source of the pixel class colours, shared by the drawing code and
% the app's colour-bar legend so the two cannot disagree. Row k of cmap is
% class k (class 0, not analysed, takes the display background colour):
%   1 neither          (grey)
%   2 channel 1 only   (magenta)
%   3 channel 2 only   (green)
%   4 both             (yellow)
%
% name1 / name2 (default 'ch 1' / 'ch 2') only change the legend wording,
% e.g. colocPixelClassColours('ER', 'mito').

if nargin < 1 || isempty(name1), name1 = 'ch 1'; end
if nargin < 2 || isempty(name2), name2 = 'ch 2'; end

cmap = [0.55 0.55 0.55
        0.90 0.20 0.80
        0.15 0.80 0.30
        1.00 0.85 0.10];

labels = {'neither', [name1 ' only'], [name2 ' only'], 'both'};
end
