function objStats = colocErObjects(erImage, cellIDIn, chER, cisternaeStats)
%COLOCEROBJECTS  The whole ER (tubules and cisternae together) as object
% tables for colocObjectOverlap: one object per connected piece of the ER
% within each cell.
%
%   objStats = colocErObjects(erImage, cellIDIn, chER)
%   objStats = colocErObjects(erImage, cellIDIn, chER, cisternaeStats)
%
% The ER is one continuous network, so there is no reason to split it into
% tubules and cisternae for object colocalisation: an organelle overlapping
% a tubule is as much "on the ER" as one overlapping a cisterna. Each cell's
% part of the ER image is labelled with 8-connectivity; in a well-traced
% cell that is usually one large object plus a few detached fragments. The
% tubule / cisterna / stream split of the overlap is made afterwards by
% colocErOverlapBreakdown.
%
% INPUTS
%   erImage        - [nY x nX x 1 x nZ x nT] logical ER image: the ER
%                    skeleton (tubule centrelines) OR the cisternae, as for
%                    organelleD2ErCompute -- not the ER mask, which is a
%                    segmentation step towards those.
%   cellIDIn       - [nY x nX x cC x cZ x cT] cell label image; the plane
%                    for channel chER is used. Pixels outside every cell
%                    (label 0) are ignored.
%   chER           - ER channel; the tables go in this channel slot.
%   cisternaeStats - (optional) {nC x nZ x nT} cisternae tables with
%                    .cisternaePixelIdxList; their pixels are added to the
%                    ER image so every analysed cisterna is part of an ER
%                    object.
%
% OUTPUT
%   objStats - {chER x nZ x nT} cell array; slot chER holds tables with the
%              columns colocObjectOverlap expects: organellePixelIdxList
%              (linear indices into the [nY nX] plane), organelleID (1..n,
%              unique within the plane) and cellID. Other slots are empty.

if nargin < 4
    cisternaeStats = {};
end
[nY, nX, ~, nZ, nT] = size(erImage);
[~, ~, cC, cZ, cT]  = size(cellIDIn);
objStats = cell(chER, nZ, nT);

for iT = 1:nT
    for iZ = 1:nZ
        m = logical(erImage(:,:,1,iZ,iT));
        if ~isempty(cisternaeStats) && chER <= size(cisternaeStats, 1)
            S = cisternaeStats{chER, min(iZ, size(cisternaeStats,2)), min(iT, size(cisternaeStats,3))};
            if istable(S) && ~isempty(S) && ismember('cisternaePixelIdxList', S.Properties.VariableNames)
                m(cat(1, S.cisternaePixelIdxList{:}, zeros(0,1))) = true;
            end
        end
        C = double(cellIDIn(:,:, min(cC,chER), min(cZ,iZ), min(cT,iT)));
        C(~m) = 0;
        cellPix = label2idx(C);

        pix = cell(0, 1); cells = zeros(0, 1);
        for g = find(~cellfun(@isempty, cellPix))
            % label the cell's part of the mask in its bounding box
            idx = cellPix{g}(:);
            [r, c] = ind2sub([nY nX], idx);
            r0 = min(r); c0 = min(c);
            h = max(r) - r0 + 1; w = max(c) - c0 + 1;
            sub = false(h, w);
            sub(sub2ind([h w], r - r0 + 1, c - c0 + 1)) = true;
            cc = bwconncomp(sub, 8);
            for k = 1:cc.NumObjects
                [rr, cq] = ind2sub([h w], cc.PixelIdxList{k});
                pix{end+1, 1} = sub2ind([nY nX], rr + r0 - 1, cq + c0 - 1); %#ok<AGROW>
                cells(end+1, 1) = g; %#ok<AGROW>
            end
        end
        ids = (1:numel(pix))';
        T = table();
        T.organellePixelIdxList = pix;
        T.organelleID = ids;
        T.cellID = cells;
        objStats{chER, iZ, iT} = T;
    end
end
end
