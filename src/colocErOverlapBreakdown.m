function [objectOut, summaryOut] = colocErOverlapBreakdown(objectOut, summaryOut, statsB, chB, cisternaeObj, chA, cellIDIn, calibration)
%COLOCEROVERLAPBREAKDOWN  Split each organelle's overlap with the whole ER
% (colocErObjects as population A) into the part on tubules, on cisternae
% and on stream cisternae.
%
%   [objectOut, summaryOut] = colocErOverlapBreakdown(objectOut, summaryOut, ...
%       statsB, chB, cisternaeObj, chA, cellIDIn, calibration)
%
% With the ER as one combined object, colocObjectOverlap no longer knows
% which part of the network an organelle sits on. This adds it back as
% columns. It only looks up pixels (no shuffles, no distance transforms),
% so it is cheap: when the streams threshold changes, call it again on the
% stored outputs to refresh the stream columns instead of re-running
% colocObjectOverlap. Existing breakdown columns are overwritten.
%
% INPUTS
%   objectOut, summaryOut - {1 x nZ x nT} outputs of colocObjectOverlap with
%                    A = colocErObjects, B = organelles.
%   statsB         - {nCb x nZ x nT} organelle stats tables (population B;
%                    organellePixelIdxList, organelleID, cellID).
%   chB            - channel slot of statsB.
%   cisternaeObj   - {nC x nZ x nT} cisternae tables from
%                    colocCisternaeObjects(..., 'all', thr) (needs
%                    organellePixelIdxList, cellID, isStream;
%                    cisternaeSpeedMax when optical flow was run), or {}
%                    when no cisternae analysis was run -> columns NaN.
%   chA            - ER channel (slot of cisternaeObj, plane of cellIDIn).
%   cellIDIn       - [nY x nX x cC x cZ x cT] cell label image; overlap is
%                    counted inside the organelle's own cell, as in
%                    colocObjectOverlap.
%   calibration    - microns per pixel.
%
% OUTPUTS
%   objectOut  - organelle (B) rows gain (ER (A) rows get NaN):
%                threshold-free
%                  overlapAreaTubules, overlapAreaCisternae (um^2; they sum
%                  to overlapArea; with the ER skeleton as the tubules,
%                  overlapAreaTubules is skeleton pixels x pixel area,
%                  i.e. centreline length x one pixel), nPartnersCisternae,
%                  maxOverlapPartnerCisternaeSpeedMax (fastest overlapping
%                  cisterna: redo the stream split at any threshold)
%                depend on the streams threshold
%                  overlapAreaStream (part of overlapAreaCisternae),
%                  nPartnersStream
%   summaryOut - per cell: nCisternae, fracBOverlappingTubules,
%                fracBOverlappingCisternae (threshold-free); nStreams,
%                fracBOverlappingStream (threshold-dependent).

haveCis = ~isempty(cisternaeObj);
[~, nZ, nT] = size(objectOut);
[nY, nX, cC, cZ, cT] = size(cellIDIn);
cal2 = calibration^2;

for iT = 1:nT
    for iZ = 1:nZ
        O = objectOut{1, iZ, iT};
        if ~istable(O) || isempty(O)
            continue
        end
        n = height(O);
        [ovT, ovC, ovS, nPC, nPS, vMax] = deal(nan(n, 1));
        SC = slot(cisternaeObj, chA, iZ, iT);
        LC = zeros(nY, nX);                 % row of SC at each cisterna pixel
        isS = false(height(SC), 1);
        speed = nan(height(SC), 1);
        if ~isempty(SC)
            for k = 1:height(SC)
                LC(SC.organellePixelIdxList{k}) = k;
            end
            isS = SC.isStream == 1;
            if ismember('cisternaeSpeedMax', SC.Properties.VariableNames)
                speed = double(SC.cisternaeSpeedMax);
            end
        end
        cellPlane = cellIDIn(:,:, min(cC,chA), min(cZ,iZ), min(cT,iT));

        if haveCis
            isB = find(strcmp(O.population, 'B'));
            SB = slot(statsB, chB, iZ, iT);
            [pix, inRun] = matchedPixels(SB, O(isB, :));
            for k = 1:numel(pix)
                row = isB(inRun(k));
                p = pix{k}(cellPlane(pix{k}) == O.cellID(row));
                lab = LC(p);
                lab = lab(lab > 0);
                nC = numel(lab);
                partners = unique(lab);
                ovC(row) = nC * cal2;
                ovS(row) = nnz(isS(lab)) * cal2;
                ovT(row) = max(0, round(O.overlapArea(row) / cal2) - nC) * cal2;
                nPC(row) = numel(partners);
                nPS(row) = nnz(isS(partners));
                if ~isempty(partners)
                    vMax(row) = max(speed(partners));
                end
            end
        end
        O.overlapAreaTubules                 = ovT;
        O.overlapAreaCisternae               = ovC;
        O.nPartnersCisternae                 = nPC;
        O.maxOverlapPartnerCisternaeSpeedMax = vMax;
        O.overlapAreaStream                  = ovS;
        O.nPartnersStream                    = nPS;
        objectOut{1, iZ, iT} = O;

        S = summaryOut{1, iZ, iT};
        if ~istable(S) || isempty(S)
            continue
        end
        m = height(S);
        [nCis, nStr, fT, fC, fS] = deal(nan(m, 1));
        if haveCis
            for i = 1:m
                g = S.cellID(i);
                if ~isempty(SC)
                    inG = SC.cellID == g;
                    nCis(i) = nnz(inG);
                    nStr(i) = nnz(inG & isS);
                else
                    [nCis(i), nStr(i)] = deal(0);
                end
                b = strcmp(O.population, 'B') & O.cellID == g;
                if any(b)
                    fT(i) = mean(ovT(b) > 0);
                    fC(i) = mean(ovC(b) > 0);
                    fS(i) = mean(ovS(b) > 0);
                end
            end
        end
        S.nCisternae                = nCis;
        S.fracBOverlappingTubules   = fT;
        S.fracBOverlappingCisternae = fC;
        S.nStreams                  = nStr;
        S.fracBOverlappingStream    = fS;
        summaryOut{1, iZ, iT} = S;
    end
end
end % colocErOverlapBreakdown


% =========================================================================
function S = slot(stats, ch, iZ, iT)
S = table.empty;
if iscell(stats) && ~isempty(stats) && ch <= size(stats, 1)
    S = stats{ch, min(iZ, size(stats, 2)), min(iT, size(stats, 3))};
end
if ~istable(S)
    S = table.empty;
end
end

function [pix, inRun] = matchedPixels(S, O)
% pixel lists of the rows of S that appear in O (same cellID + organelleID);
% inRun(k) = row of O
pix = {}; inRun = [];
if isempty(S) || isempty(O) || ~all(ismember({'cellID','organelleID','organellePixelIdxList'}, S.Properties.VariableNames))
    return
end
[tf, loc] = ismember([double(S.cellID) double(S.organelleID)], [double(O.cellID) double(O.organelleID)], 'rows');
rows  = find(tf & ~cellfun(@isempty, S.organellePixelIdxList));
inRun = loc(rows);
pix   = cellfun(@(x) double(x(:)), S.organellePixelIdxList(rows), 'UniformOutput', false);
end
