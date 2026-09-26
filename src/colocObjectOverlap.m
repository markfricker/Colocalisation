function [objectOut, pairOut, summaryOut, neighbourOut] = colocObjectOverlap(morphologyStats, cellIDIn, chA, chB, calibration, contactDistance, nShuffles, code, opts)
%COLOCOBJECTOVERLAP  Object-based colocalisation between two segmented
% organelle channels: per-pair overlap, per-object overlap/contact/nearest
% distance, k-nearest partners, per-cell summary, and an object-shuffle
% significance test.
%
%   [objectOut, pairOut, summaryOut, neighbourOut] = colocObjectOverlap( ...
%       morphologyStats, cellIDIn, chA, chB, calibration, contactDistance, ...
%       nShuffles, code, opts)
%
% Object-based counterpart to the pixel-based core (colocPixelBasedRun):
% rather than asking whether two channels' INTENSITIES co-vary, it asks
% which segmented object in channel A overlaps / touches which object in
% channel B, and by how much -- the numbers DiAna (Gilles et al. 2017,
% Methods 115:55-64) reports, e.g. "40% of peroxisomes touch a
% mitochondrion".
%
% How this differs from organelleD2ErCompute's organelleErOverlapArea:
% that measures one organelle population against the ER, a single
% continuous network represented by its 1-px skeleton centreline (+
% cisternae), so it is one-sided (fraction of the organelle only) and has
% no partner identity. Here BOTH populations are discrete, full-area
% segmented objects, so every overlapping pair has an identity and a
% fraction on each side (fracOfA, fracOfB, Jaccard), and every object has
% a partner count.
%
% SIGNIFICANCE: within each cell, channel B's objects are re-placed at
% random (shape kept, no B-B overlaps, may overlap A; shuffleObjectsInWindow
% from OrganelleDistances_sandbox) nShuffles times, and two per-cell
% statistics are recomputed:
%   - objM1        : overlapping area / total A area (object-level Manders)
%   - fracAContact : fraction of A objects within contactDistance of any B
% Monte Carlo p = (1 + #{sim >= obs}) / (nShuffles + 1) for "more than
% chance" (and <= for "less than chance"). Overlap area is symmetric, so
% one shuffle direction suffices for the overlap test.
%
% INPUTS
%   morphologyStats - {nC x nZ x nT} organelle stats tables (as produced by
%                     analyzerFeatureAnalysis). Each needs
%                     .organellePixelIdxList (linear indices into an
%                     [nY nX] plane), .organelleID and .cellID.
%   cellIDIn        - [nY x nX x cC x cZ x cT] cell label image. The plane
%                     for channel chA is used for both populations.
%   chA, chB        - channel indices of the two populations (must differ).
%   calibration     - microns per pixel.
%   contactDistance - microns; edge-to-edge distance counted as "in
%                     contact" (0 = touching/overlapping only).
%   nShuffles       - shuffles per cell for the significance test (0 =
%                     descriptive only, p-values NaN).
%   code            - filename string stamped into every row.
%   opts            - (optional) struct:
%     .dihedral    - random 90-degree rotation/flip of shuffled objects
%                    (default false).
%     .rngSeed     - reproducible shuffles (global RNG restored after).
%     .progressFcn - @(fracDone, message) -> cancelled, once per cell.
%     .kNearest    - number of nearest partners listed per object in
%                    neighbourOut (default 3; 0 = skip).
%     .labelColumns - struct of constants appended as columns to every output
%                    table, recording the run settings (e.g.
%                    struct('cisternaeClass','all','streamsThreshold',0.5)),
%                    so saved sheets from different runs stay distinguishable.
%     .partnerAttributes - cellstr of per-object columns of population A
%                    (e.g. {'isStream','cisternaeSpeedMax'} for ER cisternae)
%                    carried into the outputs so a class split can be made
%                    downstream (e.g. in R, at any threshold) instead of
%                    choosing a class before the run. For each name X:
%                      pairOut      X of the pair's A object
%                      objectOut    X (A rows, own value), nearestPartnerX
%                                   and maxOverlapPartnerX (B rows: value of
%                                   the nearest A partner by edge distance /
%                                   the max over overlapping A partners)
%                      neighbourOut partnerX (B-direction rows)
%                    If X is 'isStream', B rows also get overlapAreaStream,
%                    overlapAreaOrdinary, nPartnersStream, and summaryOut gets
%                    nStreams, fracBOverlappingStream, fracBOverlappingOrdinary.
%                    Default: {'isStream'} when any A table has isStream.
%                    NaN where not applicable.
%     .statsB      - {nCb x nZ x nT} stats array for population B, when it
%                    does not live in morphologyStats (e.g. A = ER
%                    cisternae via colocCisternaeObjects, B = organelle
%                    morphology). chB then indexes statsB, and chA may equal
%                    chB. B is the population that gets shuffled.
%
% OUTPUTS (each {1 x nZ x nT} cell array of tables; distances in microns,
% areas in microns^2)
%   objectOut  - one row per object of EITHER population: filename,
%                channel, partnerChannel, section, frame, cellID,
%                organelleID, area, overlapArea (with any partner),
%                overlapFraction, nPartners, nearestPartnerDistance
%                (edge-to-edge, 0 if overlapping, NaN if no partner in the
%                cell), contactFraction (perimeter pixels within
%                contactDistance of the partner population),
%                colocalised (overlapArea > 0), inContact
%                (nearestPartnerDistance <= contactDistance), population
%                ('A'/'B' -- unambiguous even when chA == chB), and
%                linkX1/linkY1 -> linkX2/linkY2: full-image pixel [x y] of
%                the object's own closest pixel and of the closest partner
%                pixel (the nearest-partner link; equal points when
%                overlapping; NaN without partners). See colocObjectLinks.
%   pairOut    - one row per overlapping (A,B) pair: filename, channel
%                (=chA), partnerChannel (=chB), section, frame, cellID, organelleIDA,
%                organelleIDB, overlapArea, fracOfA, fracOfB, jaccard.
%   summaryOut - one row per cell: filename, channel (=chA),
%                partnerChannel (=chB), section, frame, cellID, nA, nB,
%                nPairs, fracAOverlapping, fracBOverlapping, objM1, objM2,
%                meanPairJaccard, fracAContact, fracBContact,
%                pMoreOverlap, pLessOverlap, objM1NullMean, pMoreContact,
%                pLessContact, fracAContactNullMean, fallbackFraction,
%                nShuffles.
%   neighbourOut - DiAna-style k-nearest partners: for every object of
%                EITHER population, one row per rank 1..min(k, #partners
%                in the cell): filename, channel, partnerChannel, section,
%                frame, cellID, organelleID, rank, partnerID, edgeDistance
%                (closest pixel-to-pixel distance, 0 if overlapping; rank
%                order), centreDistance (centroid-to-centroid). Rank-1
%                edgeDistance equals objectOut's nearestPartnerDistance.
%                Direction-agnostic (distance transform), unlike the
%                normal-direction ray-cast of organelleD2OrganelleCompute.

if nargin < 9 || isempty(opts)
    opts = struct();
end
% Population B normally comes from the same stats array (a different
% channel); opts.statsB supplies it from a separate array instead (e.g. A =
% ER cisternae, B = organelle morphology), in which case chA == chB is fine.
statsB = morphologyStats;
if isfield(opts, 'statsB') && ~isempty(opts.statsB)
    statsB = opts.statsB;
elseif chA == chB
    error('colocObjectOverlap:sameChannel', 'chA and chB must be different channels.');
end
shOpts = struct('dihedral', isfield(opts,'dihedral') && ~isempty(opts.dihedral) && opts.dihedral);
progressFcn = [];
if isfield(opts, 'progressFcn'), progressFcn = opts.progressFcn; end
kNearest = 3;
if isfield(opts, 'kNearest') && ~isempty(opts.kNearest), kNearest = opts.kNearest; end
if isfield(opts, 'rngSeed') && ~isempty(opts.rngSeed)
    prevState  = rng(opts.rngSeed);
    restoreRng = onCleanup(@() rng(prevState)); %#ok<NASGU>
end

[nC_in, nZ, nT]    = size(morphologyStats);
[nY, nX, cC, cZ, cT] = size(cellIDIn);
objectOut    = cell(1, nZ, nT);
pairOut      = cell(1, nZ, nT);
summaryOut   = cell(1, nZ, nT);
neighbourOut = cell(1, nZ, nT);
if chA > nC_in || chB > size(statsB, 1)
    return
end

% Per-object attributes of population A carried into the outputs (see
% opts.partnerAttributes). Default: isStream if any A table has it, so the
% ER cisternae flag is always carried.
attrNames = {};
if isfield(opts, 'partnerAttributes') && ~isempty(opts.partnerAttributes)
    attrNames = cellstr(opts.partnerAttributes);
elseif any(cellfun(@(t) istable(t) && ismember('isStream', t.Properties.VariableNames), ...
        reshape(morphologyStats(chA,:,:), [], 1)))
    attrNames = {'isStream'};
end
labelColumns = struct();
if isfield(opts, 'labelColumns') && isstruct(opts.labelColumns)
    labelColumns = opts.labelColumns;
end

% cell count for progress reporting
jobs = zeros(0, 3);
for iT = 1:nT
    for iZ = 1:nZ
        g = cellsInPlane(morphologyStats{chA,iZ,iT}, statsB{chB,min(iZ,size(statsB,2)),min(iT,size(statsB,3))});
        jobs = [jobs; repmat([iZ iT], numel(g), 1), g(:)]; %#ok<AGROW>
    end
end
nJobs = size(jobs, 1);

for iJob = 1:nJobs
    iZ = jobs(iJob,1); iT = jobs(iJob,2); g = jobs(iJob,3);
    SA = morphologyStats{chA,iZ,iT};
    SB = statsB{chB,min(iZ,size(statsB,2)),min(iT,size(statsB,3))};
    SA = cellRows(SA, g);
    SB = cellRows(SB, g);
    cellMask = cellIDIn(:,:, min(cC,chA), min(cZ,iZ), min(cT,iT)) == g;

    meta = struct('code', code, 'chA', chA, 'chB', chB, 'iZ', iZ, 'iT', iT, 'g', g);
    [O, P, S, N] = cellOverlap(SA, SB, cellMask, [nY nX], calibration, ...
        contactDistance, nShuffles, shOpts, meta, kNearest, attrNames);
    O = addLabels(O, labelColumns);
    P = addLabels(P, labelColumns);
    S = addLabels(S, labelColumns);
    N = addLabels(N, labelColumns);
    objectOut{1,iZ,iT}    = [objectOut{1,iZ,iT};    O];
    pairOut{1,iZ,iT}      = [pairOut{1,iZ,iT};      P];
    summaryOut{1,iZ,iT}   = [summaryOut{1,iZ,iT};   S];
    neighbourOut{1,iZ,iT} = [neighbourOut{1,iZ,iT}; N];

    if ~isempty(progressFcn)
        if progressFcn(iJob / nJobs, sprintf('Object overlap: cell %d of %d', iJob, nJobs))
            return
        end
    end
end

end % colocObjectOverlap


% =========================================================================
function v = objColumn(S, name)
% column NAME of S for the rows kept by pixLists/objIDs (non-empty pixel
% lists), as double; NaN when S lacks the column
if ~istable(S) || isempty(S)
    v = zeros(0, 1);
    return
end
keep = ~cellfun(@isempty, S.organellePixelIdxList);
if ismember(name, S.Properties.VariableNames)
    v = double(S.(name)(keep));
else
    v = nan(nnz(keep), 1);
end
v = v(:);
end % objColumn

% =========================================================================
function T = addLabels(T, labels)
% append constant columns (e.g. cisternaeClass, streamsThreshold) to every
% table that has variables; a bare 0x0 table.empty is left alone so it
% still vertically concatenates with anything
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
end % addLabels

% =========================================================================
function S = cellRows(S, g)
% rows of stats table S in cell g; a missing/empty slot becomes an empty table
if istable(S) && ~isempty(S) && ismember('cellID', S.Properties.VariableNames)
    S = S(S.cellID == g, :);
else
    S = table.empty;
end
end % cellRows

% =========================================================================
function g = cellsInPlane(SA, SB)
g = [];
for S = {SA, SB}
    s = S{1};
    if istable(s) && ~isempty(s) && ismember('cellID', s.Properties.VariableNames)
        g = [g; s.cellID(:)]; %#ok<AGROW>
    end
end
g = unique(g(g > 0));
end % cellsInPlane


% =========================================================================
function [O, P, S, N] = cellOverlap(SA, SB, cellMask, imSize, cal, dContact, nShuffles, shOpts, meta, kNearest, attrNames)
%CELLOVERLAP  Single-cell worker -- see colocObjectOverlap for docs.

pixA = pixLists(SA);  idA = objIDs(SA);  nA = numel(pixA);
pixB = pixLists(SB);  idB = objIDs(SB);  nB = numel(pixB);

% --- crop to the cell + all its objects ----------------------------------
[rM, cM] = find(cellMask);
[rO, cO] = ind2sub(imSize, cat(1, pixA{:}, pixB{:}, zeros(0,1)));
rr = [rM; rO]; cc = [cM; cO];
r0 = max(1, min(rr) - 1);  r1 = min(imSize(1), max(rr) + 1);
c0 = max(1, min(cc) - 1);  c1 = min(imSize(2), max(cc) + 1);
nYc = r1 - r0 + 1;  nXc = c1 - c0 + 1;
toCrop = @(idx) cropIdx(idx, imSize, r0, c0, nYc);
cpA = cellfun(toCrop, pixA, 'UniformOutput', false);
cpB = cellfun(toCrop, pixB, 'UniformOutput', false);
window = cellMask(r0:r1, c0:c1);

LA = labelFrom(cpA, [nYc nXc]);
LB = labelFrom(cpB, [nYc nXc]);
edgeA = boundaryPixels(LA);
edgeB = boundaryPixels(LB);
areaA = reshape(cellfun(@numel, cpA), [], 1) * cal^2;
areaB = reshape(cellfun(@numel, cpB), [], 1) * cal^2;
idA = reshape(idA, [], 1);
idB = reshape(idB, [], 1);

% --- pairs ---------------------------------------------------------------
ovMask = LA > 0 & LB > 0;
if any(ovMask(:))
    pr = [LA(ovMask) LB(ovMask)];
    [upr, ~, j] = unique(pr, 'rows');
    ovPix = accumarray(j, 1);
else
    upr = zeros(0, 2); ovPix = zeros(0, 1);
end
ovArea = ovPix * cal^2;
nPairs = size(upr, 1);
fracOfA = ovArea ./ areaA(upr(:,1));
fracOfB = ovArea ./ areaB(upr(:,2));
jac     = ovArea ./ (areaA(upr(:,1)) + areaB(upr(:,2)) - ovArea);
P = table(repmat({meta.code}, nPairs, 1), repmat(meta.chA, nPairs, 1), ...
    repmat(meta.chB, nPairs, 1), repmat(meta.iZ, nPairs, 1), repmat(meta.iT, nPairs, 1), ...
    repmat(meta.g, nPairs, 1), idA(upr(:,1)), idB(upr(:,2)), ovArea, fracOfA, fracOfB, jac, ...
    'VariableNames', {'filename','channel','partnerChannel','section','frame','cellID', ...
    'organelleIDA','organelleIDB','overlapArea','fracOfA','fracOfB','jaccard'});

% --- per object, both directions ------------------------------------------
off = struct('r0', r0, 'c0', c0);
[OA, nearA] = objectRows(cpA, idA, areaA, LB > 0, edgeA, upr(:,1), ovArea, cal, dContact, meta, meta.chA, meta.chB, off, 'A');
[OB, nearB] = objectRows(cpB, idB, areaB, LA > 0, edgeB, upr(:,2), ovArea, cal, dContact, meta, meta.chB, meta.chA, off, 'B');
O = [OA; OB];

% --- k-nearest partners (both directions) ----------------------------------
% Edge distance E(a,b) = closest pixel-to-pixel distance, symmetric, so one
% distance transform per object of the SMALLER population fills the matrix.
N = table.empty;
E = zeros(nA, nB);
if (kNearest > 0 || ~isempty(attrNames)) && nA > 0 && nB > 0
    if nA <= nB
        for a = 1:nA
            D = bwdist(LA == a);
            E(a, :) = cellfun(@(x) min(D(x)), cpB);
        end
    else
        for b = 1:nB
            D = bwdist(LB == b);
            E(:, b) = cellfun(@(x) min(D(x)), cpA);
        end
    end
    E = double(E) * cal;   % bwdist returns single
end
pidxBA = zeros(0, 1);      % for B-direction neighbour rows: index of the A partner
nRowsAB = 0;
if kNearest > 0 && nA > 0 && nB > 0
    cenA = cropCentroidsRC(cpA, nYc);
    cenB = cropCentroidsRC(cpB, nYc);
    C = sqrt((cenA(:,1) - cenB(:,1)').^2 + (cenA(:,2) - cenB(:,2)').^2) * cal;
    NAB = knnRows(idA, idB, E,  C,  kNearest, meta, meta.chA, meta.chB);
    [NBA, pidxBA] = knnRows(idB, idA, E', C', kNearest, meta, meta.chB, meta.chA);
    nRowsAB = height(NAB);
    N = [NAB; NBA];
end

% --- population-A attributes carried to partners (e.g. ER cisterna isStream /
% cisternaeSpeedMax), so the stream split can be made downstream at any
% threshold. A rows get their own value; B rows get the value of their
% nearest / fastest-overlapping A partner; pairs and B-direction neighbours
% get the A partner's value.
cap = @(s) [upper(s(1)) s(2:end)];
streamInfo = [];
if nB > 0 && nA > 0
    [~, iNearA] = min(E, [], 1);        % nearest A partner of each B (edge distance)
else
    iNearA = zeros(1, nB);
end
for iAt = 1:numel(attrNames)
    name = attrNames{iAt};
    a = objColumn(SA, name);            % nA x 1, kept-row order (NaN if absent)
    P.(name) = a(upr(:,1));
    nearVal = nan(nB, 1);
    if nA > 0 && nB > 0
        nearVal = a(iNearA(:));
    end
    maxOv = nan(nB, 1);
    for b = 1:nB
        partners = upr(upr(:,2) == b, 1);
        if ~isempty(partners)
            maxOv(b) = max(a(partners));
        end
    end
    O.(name)                                = [a; nan(nB, 1)];
    O.(['nearestPartner' cap(name)])        = [nan(nA, 1); nearVal];
    O.(['maxOverlapPartner' cap(name)])     = [nan(nA, 1); maxOv];
    if istable(N) && width(N) > 0
        N.(['partner' cap(name)]) = [nan(nRowsAB, 1); a(pidxBA)];
    end
    if strcmp(name, 'isStream')
        isS = a == 1;
        ovS = accumarray(upr(:,2), ovArea .* isS(upr(:,1)),  [nB 1]);
        ovO = accumarray(upr(:,2), ovArea .* ~isS(upr(:,1)), [nB 1]);
        nS  = accumarray(upr(:,2), double(isS(upr(:,1))),   [nB 1]);
        O.overlapAreaStream   = [nan(nA, 1); ovS];
        O.overlapAreaOrdinary = [nan(nA, 1); ovO];
        O.nPartnersStream     = [nan(nA, 1); nS];
        streamInfo = struct('nStreams', nnz(isS), ...
            'fracBOverlappingStream',   safeMean(ovS > 0, nB), ...
            'fracBOverlappingOrdinary', safeMean(ovO > 0, nB));
    end
end

% --- per cell summary + shuffle test --------------------------------------
totA = sum(areaA);  totB = sum(areaB);
ovTot = nnz(ovMask) * cal^2;
objM1 = ovTot / totA;
objM2 = ovTot / totB;
fracAContact = mean(nearA <= dContact);
fracBContact = mean(nearB <= dContact);

[pMoreOv, pLessOv, m1Null, pMoreCt, pLessCt, ctNull, fbFrac] = deal(NaN);
if nShuffles > 0 && nA > 0 && nB > 0 && any(window(:))
    Amask = LA > 0;
    simM1 = zeros(nShuffles, 1);
    simCt = zeros(nShuffles, 1);
    nFb = 0;
    for s = 1:nShuffles
        [placed, fb] = shuffleObjectsInWindow(cpB, window, shOpts);
        nFb = nFb + nnz(fb);
        Bm = false(nYc, nXc);
        Bm(cat(1, placed{:})) = true;
        simM1(s) = nnz(Amask & Bm) * cal^2 / totA;
        DB = double(bwdist(Bm)) * cal;
        nearSim = cellfun(@(x) min(DB(x)), cpA);
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
if ~isempty(streamInfo)
    S.nStreams                 = streamInfo.nStreams;
    S.fracBOverlappingStream   = streamInfo.fracBOverlappingStream;
    S.fracBOverlappingOrdinary = streamInfo.fracBOverlappingOrdinary;
end

end % cellOverlap


% =========================================================================
function [O, nearest] = objectRows(cp, ids, areas, partnerMask, edgeMask, pairOwner, ovArea, cal, dContact, meta, ch, partnerCh, off, population)
n = numel(cp);
nearest = nan(n, 1);
contact = nan(n, 1);
link = nan(n, 4);      % [x y] of own closest pixel, [x y] of partner's closest pixel (full-image px)
if any(partnerMask(:))
    [Ds, IDX] = bwdist(partnerMask);
    D = double(Ds) * cal;
    szc = size(partnerMask);
    for k = 1:n
        [nearest(k), j] = min(D(cp{k}));
        pk = cp{k}(edgeMask(cp{k}));
        contact(k) = mean(D(pk) <= dContact);
        [ro, co] = ind2sub(szc, cp{k}(j));
        [rp, cq] = ind2sub(szc, double(IDX(cp{k}(j))));
        link(k, :) = [co + off.c0 - 1, ro + off.r0 - 1, cq + off.c0 - 1, rp + off.r0 - 1];
    end
end
ovObj = accumarray(pairOwner(:), ovArea(:), [n 1]);  % sum over this object's pairs
nPart = accumarray(pairOwner(:), 1, [n 1]);
O = table(repmat({meta.code}, n, 1), repmat(ch, n, 1), repmat(partnerCh, n, 1), ...
    repmat(meta.iZ, n, 1), repmat(meta.iT, n, 1), repmat(meta.g, n, 1), ids, areas, ...
    ovObj, ovObj ./ areas, nPart, nearest, contact, double(ovObj > 0), double(nearest <= dContact), ...
    'VariableNames', {'filename','channel','partnerChannel','section','frame','cellID', ...
    'organelleID','area','overlapArea','overlapFraction','nPartners', ...
    'nearestPartnerDistance','contactFraction','colocalised','inContact'});
O.population = repmat({population}, n, 1);
O.linkX1 = link(:,1);  O.linkY1 = link(:,2);
O.linkX2 = link(:,3);  O.linkY2 = link(:,4);
end % objectRows


% =========================================================================
function [N, pidx] = knnRows(ids, partnerIds, E, C, k, meta, ch, partnerCh)
%KNNROWS  One row per (object, rank) for the k nearest partners by edge
% distance (ties broken by centre distance). pidx = partner index (into
% partnerIds) of every row, for attaching partner attributes.
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
end % knnRows

function rc = cropCentroidsRC(cp, nYc)
rc = zeros(numel(cp), 2);
for k = 1:numel(cp)
    rc(k, :) = [mean(mod(cp{k} - 1, nYc) + 1), mean(floor((cp{k} - 1) / nYc) + 1)];
end
end

% =========================================================================
function p = pixLists(S)
if istable(S) && ~isempty(S)
    p = cellfun(@(x) double(x(:)), S.organellePixelIdxList, 'UniformOutput', false);
    keep = ~cellfun(@isempty, p);
    p = p(keep);
else
    p = cell(0, 1);
end
end

function id = objIDs(S)
if istable(S) && ~isempty(S)
    keep = ~cellfun(@isempty, S.organellePixelIdxList);
    if ismember('organelleID', S.Properties.VariableNames)
        id = double(S.organelleID(keep));
    else
        id = find(keep);
    end
else
    id = zeros(0, 1);
end
end

function c = cropIdx(idx, imSize, r0, c0, nYc)
[r, cc] = ind2sub(imSize, idx);
c = (r - r0 + 1) + (cc - c0) * nYc;
end

function L = labelFrom(cp, sz)
L = zeros(sz);
for k = 1:numel(cp)
    L(cp{k}) = k;
end
end

function E = boundaryPixels(L)
% object pixels with a 4-neighbour carrying a different label (incl. 0)
P = padarray(L, [1 1], 0);
C = P(2:end-1, 2:end-1);
E = C > 0 & (P(1:end-2,2:end-1) ~= C | P(3:end,2:end-1) ~= C | ...
             P(2:end-1,1:end-2) ~= C | P(2:end-1,3:end) ~= C);
end

function [pMore, pLess] = mcP(obs, sims)
n = numel(sims);
pMore = (1 + sum(sims >= obs)) / (n + 1);
pLess = (1 + sum(sims <= obs)) / (n + 1);
end

function m = safeMean(x, n)
if n == 0, m = NaN; else, m = mean(double(x)); end
end
