function lines = colocObjectLinks(objectOut, population)
%COLOCOBJECTLINKS  Nearest-partner link lines (DiAna-style) for the app's
% ROI overlay: one segment per object from its own closest pixel to the
% closest pixel of its nearest partner.
%
%   lines = colocObjectLinks(objectOut, population)
%
% INPUTS
%   objectOut  - {1 x nZ x nT} per-object tables from colocObjectOverlap
%                (needs population, linkX1/Y1/X2/Y2, nearestPartnerDistance).
%   population - 'A' or 'B': whose links to draw (e.g. 'B' = organelle ->
%                nearest ER cisterna).
%
% OUTPUT
%   lines - {1 x nZ x nT} cell array of [* x 3] polylines [x y distance],
%           NaN-separated segments with the distance (microns) repeated at
%           both ends, starting with a NaN row -- the same convention as
%           the ER-distance radial lines, so the app's existing overlay
%           code (incl. colour-by-distance) can draw them. Objects that
%           overlap their partner (distance 0) or have none are skipped.

[~, nZ, nT] = size(objectOut);
lines = cell(1, nZ, nT);
for iT = 1:nT
    for iZ = 1:nZ
        O = objectOut{1, iZ, iT};
        seg = zeros(0, 3);
        if istable(O) && ~isempty(O) && all(ismember({'population','linkX1','linkY1','linkX2','linkY2'}, O.Properties.VariableNames))
            O = O(strcmp(O.population, population), :);
            d = O.nearestPartnerDistance;
            ok = isfinite(O.linkX1) & isfinite(O.linkX2) & d > 0;
            n = nnz(ok);
            seg = nan(3 * n, 3);
            seg(1:3:end, :) = [O.linkX1(ok) O.linkY1(ok) d(ok)];
            seg(2:3:end, :) = [O.linkX2(ok) O.linkY2(ok) d(ok)];
        end
        lines{1, iZ, iT} = [nan nan nan; seg];
    end
end
end
