function rgb = colocObjectClassImage(objectOut, statsA, statsB, chA, chB, imSize)
%COLOCOBJECTCLASSIMAGE  RGB class image of an object-colocalisation run, for
% checking results by eye (e.g. whether an organelle's overlap with an ER
% "stream" is real or a segmentation artefact).
%
%   rgb = colocObjectClassImage(objectOut, statsA, statsB, chA, chB, imSize)
%
% INPUTS
%   objectOut - {1 x nZ x nT} per-object tables from colocObjectOverlap
%               (needs population, cellID, organelleID, colocalised,
%               inContact; overlapAreaStream used when present).
%   statsA    - {.. x nZ x nT} population-A tables (organellePixelIdxList,
%               organelleID, cellID; isStream optional -> stream outline).
%   statsB    - same for population B (may be the same array as statsA).
%   chA, chB  - channel slots of statsA / statsB.
%   imSize    - [nY nX].
%
% OUTPUT
%   rgb - uint8 [nY x nX x 1 x nZ x nT x 3] (the app's RGB image layout).
%         B objects are filled by status, A objects are outlined, and
%         pixels where A and B coincide are white -- colours and meaning
%         from colocObjectClassColours. Only objects that took part in the
%         run (i.e. appear in objectOut) are drawn.

cmap = uint8(round(255 * colocObjectClassColours()));
[~, nZ, nT] = size(objectOut);
nY = imSize(1); nX = imSize(2);
rgb = zeros(nY, nX, 1, nZ, nT, 3, 'uint8');

for iT = 1:nT
    for iZ = 1:nZ
        O = objectOut{1, iZ, iT};
        if ~istable(O) || isempty(O)
            continue
        end
        SA = slot(statsA, chA, iZ, iT);
        SB = slot(statsB, chB, iZ, iT);
        cls = zeros(nY, nX, 'uint8');

        % --- B objects, filled by status -----------------------------------
        OB = O(strcmp(O.population, 'B'), :);
        hasStreamCol = ismember('overlapAreaStream', OB.Properties.VariableNames);
        [pixB, inB] = matchedPixels(SB, OB);
        Bmask = false(nY, nX);
        for k = 1:numel(pixB)
            r = OB(inB(k), :);
            if hasStreamCol && r.overlapAreaStream > 0
                c = 6;
            elseif r.colocalised > 0
                c = 5;
            elseif r.inContact > 0
                c = 4;
            else
                c = 3;
            end
            cls(pixB{k}) = c;
            Bmask(pixB{k}) = true;
        end

        % --- A objects, outlined (stream in red) ---------------------------
        OA = O(strcmp(O.population, 'A'), :);
        [pixA, ~, rowsA] = matchedPixels(SA, OA);
        LA = zeros(nY, nX);
        isS = false(numel(pixA), 1);
        if ismember('isStream', SA.Properties.VariableNames)
            isS = SA.isStream(rowsA) == 1;
        end
        Amask = false(nY, nX);
        for k = 1:numel(pixA)
            LA(pixA{k}) = k;
            Amask(pixA{k}) = true;
        end
        edgeA = boundary(LA);
        streamLab = find(isS);
        cls(edgeA & ~ismember(LA, streamLab)) = 1;
        cls(edgeA &  ismember(LA, streamLab)) = 2;

        % --- overlap pixels ------------------------------------------------
        cls(Amask & Bmask) = 7;

        for ch = 1:3
            lut = [0; cmap(:, ch)];
            rgb(:, :, 1, iZ, iT, ch) = reshape(lut(double(cls) + 1), nY, nX);
        end
    end
end
end % colocObjectClassImage


% =========================================================================
function S = slot(stats, ch, iZ, iT)
S = table.empty;
if iscell(stats) && ch <= size(stats, 1)
    S = stats{ch, min(iZ, size(stats, 2)), min(iT, size(stats, 3))};
end
if ~istable(S)
    S = table.empty;
end
end

function [pix, inRun, rows] = matchedPixels(S, O)
% pixel lists of the rows of S that appear in O (same cellID + organelleID);
% inRun(k) = row of O, rows(k) = row of S
pix = {}; inRun = []; rows = [];
if isempty(S) || isempty(O) || ~all(ismember({'cellID','organelleID','organellePixelIdxList'}, S.Properties.VariableNames))
    return
end
[tf, loc] = ismember([double(S.cellID) double(S.organelleID)], [double(O.cellID) double(O.organelleID)], 'rows');
rows  = find(tf & ~cellfun(@isempty, S.organellePixelIdxList));
inRun = loc(rows);
pix   = cellfun(@(x) double(x(:)), S.organellePixelIdxList(rows), 'UniformOutput', false);
end

function E = boundary(L)
% labelled pixels with a 4-neighbour of a different label (incl. background)
P = padarray(L, [1 1], 0);
C = P(2:end-1, 2:end-1);
E = C > 0 & (P(1:end-2,2:end-1) ~= C | P(3:end,2:end-1) ~= C | ...
             P(2:end-1,1:end-2) ~= C | P(2:end-1,3:end) ~= C);
end
