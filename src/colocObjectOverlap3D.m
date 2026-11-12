function [objectOut, pairOut, summaryOut, neighbourOut] = colocObjectOverlap3D(labelsIn, statsIn, cellIDIn, chA, chB, spacing, contactDistance, nShuffles, code, opts)
%COLOCOBJECTOVERLAP3D  Object-based colocalisation between two segmented 3D
% organelle channels: the volumetric counterpart of colocObjectOverlap.
%
%   [objectOut, pairOut, summaryOut, neighbourOut] = colocObjectOverlap3D( ...
%       labelsIn, statsIn, cellIDIn, chA, chB, spacing, contactDistance, ...
%       nShuffles, code, opts)
%
% Same measures and the same table layout as colocObjectOverlap (see there
% for the meaning of every column), computed on 3D objects (e.g. Organelle
% 3D labels) with the anisotropic voxel spacing:
%   - areas become VOLUMES: the columns are named volume / overlapVolume
%     (um^3) instead of area / overlapArea, so a 3D sheet is never read as
%     um^2. fracOfA/B, jaccard, overlapFraction, objM1/2 are volume ratios.
%   - distances (nearestPartnerDistance, edgeDistance, centreDistance,
%     contactDistance) are true Euclidean distances in um between voxel
%     centres: an exact separable distance transform (in-plane bwdist per
%     section, then the minimum over sections with the Z offset added), so
%     a Z step counts dz and an XY step dx. XY voxels are assumed square
%     (dx == dy, as after the app's resampling).
%   - contactFraction is over the object's surface voxels (a 6-neighbour
%     of a different label or background).
%   - the shuffle null model re-places B objects at random in X, Y and Z
%     inside the cell VOLUME (shuffleObjectsInWindow with a 3D window:
%     shape kept, no B-B overlaps, may overlap A; dihedral flips in XY only).
%   - one result per frame for the whole stack: outputs are {1 x 1 x nT}
%     and the section column is 1. linkZ1 / linkZ2 (section index) are
%     added to the link columns.
%
% INPUTS
%   labelsIn        - [nY nX lC nZ nT] label volume (0 = background).
%                     Channel c uses label channel min(lC, c); with lC == 1
%                     both channels would be the same objects (error).
%   statsIn         - {nC x 1 x nT} object tables from
%                     analyzerFeatureAnalysis3D: objectID (= label value)
%                     and cellID per object. Only the objects listed take
%                     part (so the run matches the morphology rows).
%   cellIDIn        - [nY nX cC cZ cT] cell labels, 2D (cZ == 1, extruded
%                     through every section) or 3D (cZ == nZ); [] = no
%                     cells: the whole stack is one window (cellID 0).
%   chA, chB        - channel indices of the two populations (must differ).
%   spacing         - [dx dy dz] um (or [dxy dz]).
%   contactDistance - um, surface-to-surface distance counted as contact.
%   nShuffles       - shuffles per cell (0 = descriptive only).
%   code            - filename string stamped into every row.
%   opts            - (optional) struct: .kNearest (default 3), .rngSeed,
%                     .progressFcn @(fracDone, message) -> cancelled,
%                     .dihedral (default false), .labelColumns.
%
% OUTPUTS (each {1 x 1 x nT} cell array of tables; um, um^3)
%   objectOut, pairOut, summaryOut, neighbourOut - as colocObjectOverlap,
%   with volume / overlapVolume in place of area / overlapArea, and
%   linkZ1 / linkZ2 on objectOut.
%
% See also: colocObjectOverlap, colocObjectClassImage3D,
%           shuffleObjectsInWindow, analyzerFeatureAnalysis3D

% copyright (c) Mark Fricker 2026-
% Department of Biology, University of Oxford
%
% Nov 2026 - original version (3D counterpart of colocObjectOverlap)

if nargin < 10 || isempty(opts)
    opts = struct();
end
if chA == chB
    error('colocObjectOverlap3D:sameChannel', 'chA and chB must be different channels.');
end
[nY, nX, lC, nZ, nT] = size(labelsIn);
if lC == 1
    error('colocObjectOverlap3D:sharedLabels', ...
        'The label volume has one channel: channels %d and %d would be the same objects.', chA, chB);
end
if numel(spacing) == 2
    spacing = [spacing(1) spacing(1) spacing(2)];
end
spacing = double(spacing(:)).';
shOpts = struct('dihedral', isfield(opts,'dihedral') && ~isempty(opts.dihedral) && opts.dihedral);
progressFcn = [];
if isfield(opts, 'progressFcn'), progressFcn = opts.progressFcn; end
kNearest = 3;
if isfield(opts, 'kNearest') && ~isempty(opts.kNearest), kNearest = opts.kNearest; end
labelColumns = struct();
if isfield(opts, 'labelColumns') && isstruct(opts.labelColumns)
    labelColumns = opts.labelColumns;
end
if isfield(opts, 'rngSeed') && ~isempty(opts.rngSeed)
    prevState  = rng(opts.rngSeed);
    restoreRng = onCleanup(@() rng(prevState));
end

objectOut    = cell(1, 1, nT);
pairOut      = cell(1, 1, nT);
summaryOut   = cell(1, 1, nT);
neighbourOut = cell(1, 1, nT);
if max(chA, chB) > size(statsIn, 1)
    return
end
sz = [nY nX nZ];
hasCells = ~isempty(cellIDIn);

% cell jobs for progress reporting
jobs = zeros(0, 2);
for iT = 1:min(nT, size(statsIn, 3))
    g = cellsIn(statsIn{chA,1,iT}, statsIn{chB,1,iT}, hasCells);
    jobs = [jobs; repmat(iT, numel(g), 1), g(:)]; %#ok<AGROW>
end
nJobs = size(jobs, 1);

iTDone = 0;
for iJob = 1:nJobs
    iT = jobs(iJob, 1); g = jobs(iJob, 2);
    if iT ~= iTDone       % label lists of this frame, once
        idxA = label2idx(reshape(labelsIn(:,:,min(lC,chA),:,iT), sz));
        idxB = label2idx(reshape(labelsIn(:,:,min(lC,chB),:,iT), sz));
        if hasCells
            [~, ~, cC, cZ, cT] = size(cellIDIn);
            cid = full(reshape(cellIDIn(:,:,min(cC,chA),:,min(cT,iT)), nY, nX, cZ));
        end
        iTDone = iT;
    end
    [pixA, idA] = cellObjects(statsIn{chA,1,iT}, idxA, g);
    [pixB, idB] = cellObjects(statsIn{chB,1,iT}, idxB, g);
    if hasCells
        window = cid == g;
        if size(window, 3) == 1
            window = repmat(window, 1, 1, nZ);
        end
    else
        window = true(sz);
    end
    meta = struct('code', code, 'chA', chA, 'chB', chB, 'iZ', 1, 'iT', iT, 'g', g);
    [O, P, S, N] = cellOverlap3D(pixA, idA, pixB, idB, window, sz, spacing, ...
        contactDistance, nShuffles, shOpts, meta, kNearest);
    objectOut{1,1,iT}    = [objectOut{1,1,iT};    addLabels(O, labelColumns)];
    pairOut{1,1,iT}      = [pairOut{1,1,iT};      addLabels(P, labelColumns)];
    summaryOut{1,1,iT}   = [summaryOut{1,1,iT};   addLabels(S, labelColumns)];
    neighbourOut{1,1,iT} = [neighbourOut{1,1,iT}; addLabels(N, labelColumns)];
    if ~isempty(progressFcn)
        if progressFcn(iJob / nJobs, sprintf('3D object overlap: cell %d of %d', iJob, nJobs))
            return
        end
    end
end
end % colocObjectOverlap3D


% =========================================================================
function g = cellsIn(SA, SB, hasCells)
% cells holding objects of either population; without cell labels the
% whole stack is one window, cellID 0
g = [];
for S = {SA, SB}
    s = S{1};
    if istable(s) && ~isempty(s) && ismember('cellID', s.Properties.VariableNames)
        g = [g; double(s.cellID(:))]; %#ok<AGROW>
    end
end
if hasCells
    g = unique(g(g > 0));
else
    g = unique(g(g == 0));
end
end

function [pix, ids] = cellObjects(S, idx, g)
% voxel lists and object IDs of the table rows in cell g
pix = cell(0, 1); ids = zeros(0, 1);
if ~istable(S) || isempty(S) || ~all(ismember({'cellID','objectID'}, S.Properties.VariableNames))
    return
end
ids = double(S.objectID(S.cellID == g));
ids = ids(ids >= 1 & ids <= numel(idx));
pix = reshape(idx(ids), [], 1);
keep = ~cellfun(@isempty, pix);
pix = cellfun(@(x) double(x(:)), pix(keep), 'UniformOutput', false);
ids = ids(keep);
end

function T = addLabels(T, labels)
if ~istable(T) || width(T) == 0
    return
end
n = height(T);
for f = fieldnames(labels)'
    val = labels.(f{1});
    if ischar(val) || isstring(val)
        T.(f{1}) = repmat({char(val)}, n, 1);
    else
        T.(f{1}) = repmat(double(val), n, 1);
    end
end
end


% =========================================================================
function [O, P, S, N] = cellOverlap3D(pixA, idA, pixB, idB, window, sz, spacing, dContact, nShuffles, shOpts, meta, kNearest)
%CELLOVERLAP3D  Single-cell worker -- see colocObjectOverlap3D.
nA = numel(pixA); nB = numel(pixB);
vVox = prod(spacing);

% --- crop to the cell volume + all its objects ----------------------------
[rM, cM, zM] = ind2sub(sz, find(window));
[rO, cO, zO] = ind2sub(sz, cat(1, pixA{:}, pixB{:}, zeros(0,1)));
rr = [rM; rO]; cc = [cM; cO]; zz = [zM; zO];
if isempty(rr)
    rr = 1; cc = 1; zz = 1;
end
lo = max([1 1 1], [min(rr) min(cc) min(zz)] - 1);
hi = min(sz,      [max(rr) max(cc) max(zz)] + 1);
szc = hi - lo + 1;
toCrop = @(idx) cropIdx(idx, sz, lo, szc);
cpA = cellfun(toCrop, pixA, 'UniformOutput', false);
cpB = cellfun(toCrop, pixB, 'UniformOutput', false);
win = window(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3));

LA = labelFrom(cpA, szc);
LB = labelFrom(cpB, szc);
volA = reshape(cellfun(@numel, cpA), [], 1) * vVox;
volB = reshape(cellfun(@numel, cpB), [], 1) * vVox;
idA = reshape(idA, [], 1);
idB = reshape(idB, [], 1);

% --- pairs -----------------------------------------------------------------
ovMask = LA > 0 & LB > 0;
if any(ovMask(:))
    [upr, ~, j] = unique([LA(ovMask) LB(ovMask)], 'rows');
    ovVox = accumarray(j, 1);
else
    upr = zeros(0, 2); ovVox = zeros(0, 1);
end
ovVol = ovVox * vVox;
nPairs = size(upr, 1);
fracOfA = ovVol ./ volA(upr(:,1));
fracOfB = ovVol ./ volB(upr(:,2));
jac     = ovVol ./ (volA(upr(:,1)) + volB(upr(:,2)) - ovVol);
P = table(repmat({meta.code}, nPairs, 1), repmat(meta.chA, nPairs, 1), ...
    repmat(meta.chB, nPairs, 1), repmat(meta.iZ, nPairs, 1), repmat(meta.iT, nPairs, 1), ...
    repmat(meta.g, nPairs, 1), idA(upr(:,1)), idB(upr(:,2)), ovVol, fracOfA, fracOfB, jac, ...
    'VariableNames', {'filename','channel','partnerChannel','section','frame','cellID', ...
    'organelleIDA','organelleIDB','overlapVolume','fracOfA','fracOfB','jaccard'});

% --- per object, both directions -------------------------------------------
edgeA = surfaceVoxels(LA);
edgeB = surfaceVoxels(LB);
[OA, nearA] = objectRows(cpA, idA, volA, LB > 0, edgeA, upr(:,1), ovVol, spacing, dContact, meta, meta.chA, meta.chB, lo, szc, 'A');
[OB, nearB] = objectRows(cpB, idB, volB, LA > 0, edgeB, upr(:,2), ovVol, spacing, dContact, meta, meta.chB, meta.chA, lo, szc, 'B');
O = [OA; OB];

% --- k-nearest partners: surface-to-surface distance matrix ----------------
N = table.empty;
if kNearest > 0 && nA > 0 && nB > 0
    E = zeros(nA, nB);
    if nA <= nB
        qB = cat(1, cpB{:}); ownB = ownerOf(cpB);
        for a = 1:nA
            E(a, :) = accumarray(ownB, distQuery(LA == a, qB, szc, spacing), [nB 1], @min).';
        end
    else
        qA = cat(1, cpA{:}); ownA = ownerOf(cpA);
        for b = 1:nB
            E(:, b) = accumarray(ownA, distQuery(LB == b, qA, szc, spacing), [nA 1], @min);
        end
    end
    cenA = centroidsUm(cpA, szc, spacing);
    cenB = centroidsUm(cpB, szc, spacing);
    C = sqrt((cenA(:,1) - cenB(:,1)').^2 + (cenA(:,2) - cenB(:,2)').^2 + (cenA(:,3) - cenB(:,3)').^2);
    N = [knnRows(idA, idB, E,  C,  kNearest, meta, meta.chA, meta.chB); ...
         knnRows(idB, idA, E', C', kNearest, meta, meta.chB, meta.chA)];
end

% --- per cell summary + shuffle test ---------------------------------------
totA = sum(volA);  totB = sum(volB);
ovTot = nnz(ovMask) * vVox;
objM1 = ovTot / totA;
objM2 = ovTot / totB;
fracAContact = mean(nearA <= dContact);
fracBContact = mean(nearB <= dContact);

[pMoreOv, pLessOv, m1Null, pMoreCt, pLessCt, ctNull, fbFrac] = deal(NaN);
if nShuffles > 0 && nA > 0 && nB > 0 && any(win(:))
    Amask = LA > 0;
    qA = cat(1, cpA{:}); ownA = ownerOf(cpA);
    simM1 = zeros(nShuffles, 1);
    simCt = zeros(nShuffles, 1);
    nFb = 0;
    for s = 1:nShuffles
        [placed, fb] = shuffleObjectsInWindow(cpB, win, shOpts);
        nFb = nFb + nnz(fb);
        Bm = false(szc);
        Bm(cat(1, placed{:})) = true;
        simM1(s) = nnz(Amask & Bm) * vVox / totA;
        nearSim = accumarray(ownA, distQuery(Bm, qA, szc, spacing), [nA 1], @min);
        simCt(s) = mean(nearSim <= dContact);
    end
    [pMoreOv, pLessOv] = mcP(objM1, simM1);
    [pMoreCt, pLessCt] = mcP(fracAContact, simCt);
    m1Null = mean(simM1);
    ctNull = mean(simCt);
    fbFrac = nFb / (nShuffles * nB);
end
if nA == 0 || nB == 0
    [objM1, objM2, fracAContact, fracBContact] = deal(NaN);
    if nA == 0, objM2 = 0; end
    if nB == 0, objM1 = 0; end
end

S = table({meta.code}, meta.chA, meta.chB, meta.iZ, meta.iT, meta.g, nA, nB, nPairs, ...
    'VariableNames', {'filename','channel','partnerChannel','section','frame','cellID','nA','nB','nPairs'});
S.fracAOverlapping = safeMean(ismember((1:nA)', upr(:,1)), nA);
S.fracBOverlapping = safeMean(ismember((1:nB)', upr(:,2)), nB);
S.objM1 = objM1;
S.objM2 = objM2;
S.meanPairJaccard = safeMean(jac, nPairs);
S.fracAContact = fracAContact;
S.fracBContact = fracBContact;
S.pMoreOverlap = pMoreOv;
S.pLessOverlap = pLessOv;
S.objM1NullMean = m1Null;
S.pMoreContact = pMoreCt;
S.pLessContact = pLessCt;
S.fracAContactNullMean = ctNull;
S.fallbackFraction = fbFrac;
S.nShuffles = nShuffles;
end % cellOverlap3D


% =========================================================================
function [O, nearest] = objectRows(cp, ids, vols, partnerMask, edgeMask, pairOwner, ovVol, spacing, dContact, meta, ch, partnerCh, lo, szc, population)
n = numel(cp);
nearest = nan(n, 1);
contact = nan(n, 1);
link = nan(n, 6);      % [x y z] own closest voxel, [x y z] partner's closest voxel (full stack)
if any(partnerMask(:)) && n > 0
    q = cat(1, cp{:});
    [d, ni] = distQuery(partnerMask, q, szc, spacing);
    onEdge = edgeMask(q);
    first = [0; cumsum(cellfun(@numel, cp(:)))];
    for k = 1:n
        rows = first(k)+1:first(k+1);
        [nearest(k), j] = min(d(rows));
        e = rows(onEdge(rows));
        if isempty(e), e = rows; end
        contact(k) = mean(d(e) <= dContact);
        [ro, co, zo] = ind2sub(szc, q(rows(j)));
        [rp, cq, zp] = ind2sub(szc, ni(rows(j)));
        link(k, :) = [co ro zo cq rp zp] + [lo([2 1 3]) lo([2 1 3])] - 1;
    end
end
ovObj = accumarray(pairOwner(:), ovVol(:), [n 1]);
nPart = accumarray(pairOwner(:), 1, [n 1]);
O = table(repmat({meta.code}, n, 1), repmat(ch, n, 1), repmat(partnerCh, n, 1), ...
    repmat(meta.iZ, n, 1), repmat(meta.iT, n, 1), repmat(meta.g, n, 1), ids, vols, ...
    ovObj, ovObj ./ vols, nPart, nearest, contact, double(ovObj > 0), double(nearest <= dContact), ...
    'VariableNames', {'filename','channel','partnerChannel','section','frame','cellID', ...
    'organelleID','volume','overlapVolume','overlapFraction','nPartners', ...
    'nearestPartnerDistance','contactFraction','colocalised','inContact'});
O.population = repmat({population}, n, 1);
O.linkX1 = link(:,1);  O.linkY1 = link(:,2);  O.linkZ1 = link(:,3);
O.linkX2 = link(:,4);  O.linkY2 = link(:,5);  O.linkZ2 = link(:,6);
end % objectRows


% =========================================================================
function [d, nearIdx] = distQuery(mask, q, szc, spacing)
%DISTQUERY  Exact anisotropic Euclidean distance (um) from each query voxel
% q (linear indices into an szc volume) to the nearest true voxel of MASK,
% and that voxel's linear index. Separable: the in-plane distance transform
% of every section holding mask voxels, then the minimum over those
% sections of (in-plane distance)^2 + (Z offset * dz)^2. Inf / NaN when
% MASK is empty.
q = q(:);
d = inf(numel(q), 1);
nearIdx = nan(numel(q), 1);
if isempty(q), return, end
nP = szc(1) * szc(2);
qp = mod(q - 1, nP) + 1;               % in-plane index
qz = floor((q - 1) / nP) + 1;          % section
qr = mod(qp - 1, szc(1)) + 1;
qc = floor((qp - 1) / szc(1)) + 1;
dx = spacing(1); dy = spacing(2); dz = spacing(3);
ks = find(reshape(any(any(mask, 1), 2), 1, []));
for k = ks
    % nearest in-plane mask pixel (bwdist's feature transform); the
    % distance itself is recomputed in double from the coordinates
    % (bwdist's distance is single precision)
    [~, Ik] = bwdist(mask(:,:,k));
    ip = double(Ik(qp));
    ir = mod(ip - 1, szc(1)) + 1;
    ic = floor((ip - 1) / szc(1)) + 1;
    d2 = ((qr - ir) * dy).^2 + ((qc - ic) * dx).^2 + ((qz - k) * dz).^2;
    better = d2 < d;
    d(better) = d2(better);
    nearIdx(better) = ip(better) + (k - 1) * nP;
end
d = sqrt(d);
end

function own = ownerOf(cp)
% object index of every voxel of cat(1, cp{:})
own = reshape(repelem((1:numel(cp))', cellfun(@numel, cp(:))), [], 1);  % repelem(1,n) is a row
end

function cen = centroidsUm(cp, szc, spacing)
% [x y z] centroids in um (crop frame; only differences are used)
cen = zeros(numel(cp), 3);
for k = 1:numel(cp)
    [r, c, z] = ind2sub(szc, cp{k});
    cen(k, :) = [mean(c) * spacing(1), mean(r) * spacing(2), mean(z) * spacing(3)];
end
end

function E = surfaceVoxels(L)
% labelled voxels with a 6-neighbour of a different label (incl. background)
P = padarray(L, [1 1 1], 0);
C = P(2:end-1, 2:end-1, 2:end-1);
E = C > 0 & (P(1:end-2,2:end-1,2:end-1) ~= C | P(3:end,2:end-1,2:end-1) ~= C | ...
             P(2:end-1,1:end-2,2:end-1) ~= C | P(2:end-1,3:end,2:end-1) ~= C | ...
             P(2:end-1,2:end-1,1:end-2) ~= C | P(2:end-1,2:end-1,3:end) ~= C);
end

function c = cropIdx(idx, sz, lo, szc)
[r, cc, z] = ind2sub(sz, idx);
c = sub2ind(szc, r - lo(1) + 1, cc - lo(2) + 1, z - lo(3) + 1);
end

function L = labelFrom(cp, sz)
L = zeros(sz);
for k = 1:numel(cp)
    L(cp{k}) = k;
end
end

function [N, pidx] = knnRows(ids, partnerIds, E, C, k, meta, ch, partnerCh)
% one row per (object, rank) for the k nearest partners by surface distance
% (ties broken by centre distance) -- as colocObjectOverlap's knnRows
n  = numel(ids);
kk = min(k, numel(partnerIds));
rows = cell(n, 1);
pidxCell = cell(n, 1);
for i = 1:n
    [~, order] = sortrows([E(i,:)' C(i,:)']);
    j = order(1:kk);
    pidxCell{i} = j(:);
    rows{i} = table(repmat(ids(i), kk, 1), (1:kk)', partnerIds(j), E(i,j)', C(i,j)', ...
        'VariableNames', {'organelleID','rank','partnerID','edgeDistance','centreDistance'});
end
pidx = cat(1, pidxCell{:}, zeros(0, 1));
N = cat(1, rows{:});
m = height(N);
N = [table(repmat({meta.code}, m, 1), repmat(ch, m, 1), repmat(partnerCh, m, 1), ...
    repmat(meta.iZ, m, 1), repmat(meta.iT, m, 1), repmat(meta.g, m, 1), ...
    'VariableNames', {'filename','channel','partnerChannel','section','frame','cellID'}), N];
end

function [pMore, pLess] = mcP(obs, sims)
n = numel(sims);
pMore = (1 + sum(sims >= obs)) / (n + 1);
pLess = (1 + sum(sims <= obs)) / (n + 1);
end

function m = safeMean(x, n)
if n == 0, m = NaN; else, m = mean(double(x)); end
end
