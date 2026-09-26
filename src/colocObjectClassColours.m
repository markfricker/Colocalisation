function [cmap, labels] = colocObjectClassColours(nameA, nameB)
%COLOCOBJECTCLASSCOLOURS  Palette + legend labels for colocObjectClassImage.
%
%   [cmap, labels] = colocObjectClassColours(nameA, nameB)
%
% Single source of truth for the object-colocalisation class image colours,
% shared by colocObjectClassImage (drawing) and the app's colour-bar legend,
% so the two can never disagree. Row k of cmap is class k:
%   1 A outline            (cyan)     e.g. ordinary ER cisterna
%   2 A outline, stream    (red)      A object flagged isStream
%   3 B isolated           (grey)     no partner within the contact distance
%   4 B in contact         (yellow)   within contactDistance, no overlap
%   5 B overlaps A         (green)    overlaps an ordinary A object
%   6 B overlaps stream    (magenta)  overlaps a stream A object (takes
%                                     precedence: the case to be sceptical of)
%   7 overlap pixels       (white)    A and B pixels coincide
%
% nameA / nameB (default 'A' / 'B') only change the legend wording, e.g.
% colocObjectClassColours('cisterna', 'organelle').

if nargin < 1 || isempty(nameA), nameA = 'A'; end
if nargin < 2 || isempty(nameB), nameB = 'B'; end

cmap = [0   1   1
        1   0   0
        0.5 0.5 0.5
        1   1   0
        0   1   0
        1   0   1
        1   1   1];

labels = {[nameA ' outline'], ...
          [nameA ' outline (stream)'], ...
          [nameB ' isolated'], ...
          [nameB ' in contact'], ...
          [nameB ' overlaps ' nameA], ...
          [nameB ' overlaps stream'], ...
          'overlap pixels'};
end
