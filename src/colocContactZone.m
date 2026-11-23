function [cellOut, objectOut, zoneOut] = colocContactZone(SA, chA, SB, chB, test, cellIDIn, calibration, params, code, opts)
%COLOCCONTACTZONE  Three-way contact-zone test: is a test population (or a
% test channel's intensity) enriched where two other populations meet?
%
%   [cellOut, objectOut, zoneOut] = colocContactZone(SA, chA, SB, chB, ...
%       test, cellIDIn, calibration, params, code, opts)
%
% The contact zone of populations A and B (e.g. ER and mitochondria) is
%
%   Z = pixels within distanceA of A  AND  within distanceB of B,
%
% inside the cell. Distances are pixel-centre to nearest object pixel
% (bwdist), in microns, so distance 0 makes Z the overlap of A and B -- the
% same "0 = touching/overlapping only" rule as colocObjectOverlap's
% contactDistance. The cell is split into four disjoint classes: zone,
% near A only, near B only, neither.
%
% Two modes:
%   'objects'   - test = a segmented population (e.g. contact-site puncta).
%                 Statistic: fraction of test objects touching Z.
%   'intensity' - test = one intensity image (no segmentation needed).
%                 Statistic: mean test intensity in Z / mean in the rest
%                 of the cell (enrichment).
%
% Null model (params.shuffle):
%   'partner' (default) - the B objects are re-placed at random inside the
%                 cell (shape kept, no B-B overlap; shuffleObjectsInWindow
%                 from OrganelleDistances_sandbox) and Z is rebuilt. A and
%                 the test stay put, so a test that merely follows A (e.g.
%                 an ER membrane protein) is NOT called enriched: the
%                 shuffled zones are next to A too. Needs B to be separate
%                 objects -- put a network population (ER) in A.
%   'test'      - objects mode only: the test objects are re-placed at
%                 random inside the cell, Z fixed. Controls for test
%                 density, not for association with A.
% Monte Carlo p = (1 + #{sim >= obs}) / (nShuffles + 1) for "more than
% chance" (pMore) and <= for pLess.
%
% INPUTS
%   SA, SB      - {nC x nZ x nT} stats arrays (organellePixelIdxList,
%                 organelleID, cellID), as for colocObjectOverlap; chA / chB
%                 pick the slot. Z/T slots of SB and the test are clamped
%                 (min) to their own size.
%   test        - 'objects': {nC x nZ x nT} stats array, slot
%                 params.testChannel.
%                 'intensity': numeric image [nY nX nZ nT] (or
%                 [nY nX 1 nZ nT]) of the test channel.
%   cellIDIn    - [nY x nX x cC x cZ x cT] cell label image; the plane of
%                 channel chA defines the cells.
%   calibration - microns per pixel.
%   params      - struct:
%     .distanceA, .distanceB - microns (>= 0)
%     .mode        - 'objects' | 'intensity'
%     .testChannel - slot of the test stats array (objects mode; also
%                    stamped as testChannel in intensity mode, default NaN)
%     .nShuffles   - shuffles per cell (0 = descriptive only, p NaN)
%     .shuffle     - 'partner' (default) | 'test'
%   code        - filename string stamped into every row.
%   opts        - (optional) struct:
%     .dihedral     - random 90-degree rotation/flip of shuffled objects
%     .rngSeed      - reproducible shuffles (global RNG restored after)
%     .progressFcn  - @(fracDone, message) -> cancelled, once per cell
%     .labelColumns - struct of constants appended as columns to every
%                     output table (e.g. populationA/B/test names)
%
% OUTPUTS
%   cellOut   - {1 x nZ x nT} tables, one row per cell (same columns in
%               both modes, NaN where a column does not apply):
%               filename, channelA, channelB, testChannel, section, frame,
%               cellID, mode, shuffle, distanceA, distanceB, nA, nB, nTest,
%               cellArea, zoneArea, zoneFraction (zone / cell area),
%               nearAOnlyArea, nearBOnlyArea,
%               -- objects mode --
%               fracTestInZone (test objects touching Z), fracTestAreaInZone
%               (test pixels in Z / all test pixels), areaEnrichment
%               (fracTestAreaInZone / zoneFraction: 1 = as if uniform),
%               -- intensity mode --
%               meanZone, meanNearAOnly, meanNearBOnly, meanNeither,
%               meanOutsideZone, fracSignalInZone (summed test intensity
%               in Z / in the cell),
%               -- test (both modes) --
%               statistic (fracTestInZone or zone/outside intensity ratio),
%               statisticNullMean, enrichment (statistic / null mean in
%               objects mode; meanZone / meanOutsideZone in intensity
%               mode), pMore, pLess, fallbackFraction (shuffled objects
%               left in place), nShuffles.
%               Areas in microns^2.
%   objectOut - {1 x nZ x nT} tables, one row per test object (objects
%               mode; empty tables in intensity mode): filename, channelA,
%               channelB, testChannel, section, frame, cellID, organelleID,
%               area, fracInZone, inContactZone (any pixel in Z),
%               distanceToA, distanceToB (closest object pixel to the
%               nearest A / B pixel, microns; 0 if overlapping; NaN if the
%               cell has no A / B).
%   zoneOut   - [nY x nX x nZ x nT] uint8 class image: 0 outside the
%               cells, 1 neither, 2 near A only, 3 near B only, 4 zone.

if nargin < 10 || isempty(opts)
    opts = struct();
end
P = normaliseParams(params);
if strcmp(P.mode, 'intensity') && strcmp(P.shuffle, 'test')
    error('colocContactZone:shuffleTestIntensity', ...
        'Shuffle ''test'' needs test objects; intensity mode shuffles the partner (B).');
end
shOpts = struct('dihedral', isfield(opts,'dihedral') && ~isempty(opts.dihedral) && opts.dihedral);
progressFcn = [];
if isfield(opts, 'progressFcn'), progressFcn = opts.progressFcn; end
labelColumns = struct();
if isfield(opts, 'labelColumns') && isstruct(opts.labelColumns)
    labelColumns = opts.labelColumns;
end
if isfield(opts, 'rngSeed') && ~isempty(opts.rngSeed)
    prevState  = rng(opts.rngSeed);
    restoreRng = onCleanup(@() rng(prevState));
end

[nY, nX, cC, cZ, cT] = size(cellIDIn);
[~, nZ, nT] = size(SA);
isObj = strcmp(P.mode, 'objects');
if ~isObj
    I = test;
    if ndims(I) == 5 || (size(I, 3) == 1 && ndims(I) > 3)
        I = reshape(I, nY, nX, size(I, 4), size(I, 5));
    end
    if size(I, 1) ~= nY || size(I, 2) ~= nX
        error('colocContactZone:imageSize', 'The test image must be %d x %d.', nY, nX);
    end
end

cellOut   = cell(1, nZ, nT);
objectOut = cell(1, nZ, nT);
zoneOut   = zeros(nY, nX, nZ, nT, 'uint8');

% (plane, cell) jobs: every labelled cell of the plane
jobs = zeros(0, 3);
for iT = 1:nT
    for iZ = 1:nZ
        L = cellIDIn(:,:, min(cC,chA), min(cZ,iZ), min(cT,iT));
        g = unique(L(L > 0));
        jobs = [jobs; repmat([iZ iT], numel(g), 1), double(g(:))]; %#ok<AGROW>
    end
end
nJobs = size(jobs, 1);

for iJob = 1:nJobs
    iZ = jobs(iJob,1); iT = jobs(iJob,2); g = jobs(iJob,3);
    tA = cellRows(slot(SA, chA, iZ, iT), g);
    tB = cellRows(slot(SB, chB, iZ, iT), g);
    if isObj
        tT = cellRows(slot(test, P.testChannel, iZ, iT), g);
        Ic = [];
    else
        tT = table.empty;
        Ic = double(I(:,:, min(size(I,3),iZ), min(size(I,4),iT)));
    end
    cellMask = cellIDIn(:,:, min(cC,chA), min(cZ,iZ), min(cT,iT)) == g;
    meta = struct('code', code, 'chA', chA, 'chB', chB, 'chT', P.testChannel, ...
        'iZ', iZ, 'iT', iT, 'g', g);
    [C, O, cls, box] = cellZone(tA, tB, tT, Ic, cellMask, [nY nX], calibration, P, shOpts, meta);
    plane = zoneOut(:,:,iZ,iT);
    sub = plane(box(1):box(2), box(3):box(4));
    sub(cls > 0) = cls(cls > 0);
    plane(box(1):box(2), box(3):box(4)) = sub;
    zoneOut(:,:,iZ,iT) = plane;
    cellOut{1,iZ,iT}   = [cellOut{1,iZ,iT};   addLabels(C, labelColumns)];
    objectOut{1,iZ,iT} = [objectOut{1,iZ,iT}; addLabels(O, labelColumns)];
    if ~isempty(progressFcn)
        if progressFcn(iJob / nJobs, sprintf('Contact zone: cell %d of %d', iJob, nJobs))
            return
        end
    end
end

end % colocContactZone


% =========================================================================
function [C, O, cls, box] = cellZone(tA, tB, tT, I, cellMask, imSize, cal, P, shOpts, meta)
%CELLZONE  Single-cell worker -- see colocContactZone for docs.

pixA = pixLists(tA);  pixB = pixLists(tB);  pixT = pixLists(tT);
idT  = objIDs(tT);
nA = numel(pixA);  nB = numel(pixB);  nTest = numel(pixT);
isObj = strcmp(P.mode, 'objects');

% --- crop to the cell + all its objects ----------------------------------
[rM, cM] = find(cellMask);
[rO, cO] = ind2sub(imSize, cat(1, pixA{:}, pixB{:}, pixT{:}, zeros(0,1)));
rr = [rM; rO]; cc = [cM; cO];
r0 = max(1, min(rr) - 1);  r1 = min(imSize(1), max(rr) + 1);
c0 = max(1, min(cc) - 1);  c1 = min(imSize(2), max(cc) + 1);
box = [r0 r1 c0 c1];
szc = [r1 - r0 + 1, c1 - c0 + 1];
toCrop = @(idx) (mod(idx - 1, imSize(1)) + 1 - r0 + 1) + (floor((idx - 1) / imSize(1)) + 1 - c0) * szc(1);
cpA = cellfun(toCrop, pixA, 'UniformOutput', false);
cpB = cellfun(toCrop, pixB, 'UniformOutput', false);
cpT = cellfun(toCrop, pixT, 'UniformOutput', false);
win = cellMask(r0:r1, c0:c1);
if ~isempty(I), I = I(r0:r1, c0:c1); end

Am = maskFrom(cpA, szc);
Bm = maskFrom(cpB, szc);
DA = distMap(Am, cal);
DB = distMap(Bm, cal);
nearA = DA <= P.distanceA;
nearB = DB <= P.distanceB;
Z = nearA & nearB & win;

cls = zeros(szc, 'uint8');
cls(win) = 1;
cls(win & nearA & ~nearB) = 2;
cls(win & nearB & ~nearA) = 3;
cls(Z) = 4;

cellArea = nnz(win) * cal^2;
zoneArea = nnz(Z) * cal^2;
zoneFrac = nnz(Z) / max(nnz(win), 1);
if ~any(win(:)), zoneFrac = NaN; end

% --- observed statistic + per-object rows -------------------------------
c = struct('fracTestInZone', NaN, 'fracTestAreaInZone', NaN, 'areaEnrichment', NaN, ...
    'meanZone', NaN, 'meanNearAOnly', NaN, 'meanNearBOnly', NaN, 'meanNeither', NaN, ...
    'meanOutsideZone', NaN, 'fracSignalInZone', NaN);
O = table.empty;
if isObj
    areaT = reshape(cellfun(@numel, cpT), [], 1) * cal^2;
    fracIn = reshape(cellfun(@(x) mean(Z(x)), cpT), [], 1);
    dA = reshape(cellfun(@(x) min(DA(x)), cpT), [], 1);
    dB = reshape(cellfun(@(x) min(DB(x)), cpT), [], 1);
    dA(isinf(dA)) = NaN;   % no A / B in the cell
    dB(isinf(dB)) = NaN;
    inZ = double(fracIn > 0);
    if nTest > 0
        c.fracTestInZone = mean(inZ);
        allT = cat(1, cpT{:});
        c.fracTestAreaInZone = mean(Z(allT));
        c.areaEnrichment = c.fracTestAreaInZone / zoneFrac;
    end
    stat = c.fracTestInZone;
    n = nTest;
    O = table(repmat({meta.code}, n, 1), repmat(meta.chA, n, 1), repmat(meta.chB, n, 1), ...
        repmat(meta.chT, n, 1), repmat(meta.iZ, n, 1), repmat(meta.iT, n, 1), ...
        repmat(meta.g, n, 1), idT, areaT, fracIn, inZ, dA, dB, ...
        'VariableNames', {'filename','channelA','channelB','testChannel','section','frame', ...
        'cellID','organelleID','area','fracInZone','inContactZone','distanceToA','distanceToB'});
else
    ok = win & ~isnan(I);
    c.meanZone        = safeMean(I(Z & ok));
    c.meanNearAOnly   = safeMean(I(cls == 2 & ok));
    c.meanNearBOnly   = safeMean(I(cls == 3 & ok));
    c.meanNeither     = safeMean(I(cls == 1 & ok));
    c.meanOutsideZone = safeMean(I(ok & ~Z));
    tot = sum(I(ok));
    if tot ~= 0, c.fracSignalInZone = sum(I(Z & ok)) / tot; end
    stat = c.meanZone / c.meanOutsideZone;
end

% --- shuffle test ---------------------------------------------------------
[nullMean, pMore, pLess, fbFrac] = deal(NaN);
movers = cpB;
if strcmp(P.shuffle, 'test'), movers = cpT; end
canShuffle = P.nShuffles > 0 && ~isnan(stat) && ~isempty(movers) && any(win(:)) ...
    && nA > 0 && nB > 0 && (~isObj || nTest > 0);
if canShuffle
    sims = nan(P.nShuffles, 1);
    nFb = 0;
    for s = 1:P.nShuffles
        [placed, fb] = shuffleObjectsInWindow(movers, win, shOpts);
        nFb = nFb + nnz(fb);
        if strcmp(P.shuffle, 'test')
            Zs = Z;
            cpTs = placed;
        else
            Bs = maskFrom(placed, szc);
            Zs = nearA & (distMap(Bs, cal) <= P.distanceB) & win;
            cpTs = cpT;
        end
        if isObj
            sims(s) = mean(cellfun(@(x) any(Zs(x)), cpTs));
        else
            ok = win & ~isnan(I);
            sims(s) = safeMean(I(Zs & ok)) / safeMean(I(ok & ~Zs));
        end
    end
    sims = sims(~isnan(sims));
    if ~isempty(sims)
        nullMean = mean(sims);
        pMore = (1 + sum(sims >= stat)) / (numel(sims) + 1);
        pLess = (1 + sum(sims <= stat)) / (numel(sims) + 1);
    end
    fbFrac = nFb / (P.nShuffles * numel(movers));
end
if isObj
    enrichment = stat / nullMean;
else
    enrichment = stat;
end

nTestCol = NaN;
if isObj, nTestCol = nTest; end
C = table({meta.code}, meta.chA, meta.chB, meta.chT, meta.iZ, meta.iT, meta.g, ...
    {P.mode}, {P.shuffle}, P.distanceA, P.distanceB, nA, nB, nTestCol, ...
    cellArea, zoneArea, zoneFrac, nnz(cls == 2) * cal^2, nnz(cls == 3) * cal^2, ...
    'VariableNames', {'filename','channelA','channelB','testChannel','section','frame', ...
    'cellID','mode','shuffle','distanceA','distanceB','nA','nB','nTest', ...
    'cellArea','zoneArea','zoneFraction','nearAOnlyArea','nearBOnlyArea'});
for f = fieldnames(c)'
    C.(f{1}) = c.(f{1});
end
C.statistic         = stat;
C.statisticNullMean = nullMean;
C.enrichment        = enrichment;
C.pMore             = pMore;
C.pLess             = pLess;
C.fallbackFraction  = fbFrac;
C.nShuffles         = P.nShuffles;

end % cellZone


% =========================================================================
function P = normaliseParams(p)
P = struct('distanceA', 0, 'distanceB', 0, 'mode', 'objects', 'testChannel', NaN, ...
    'nShuffles', 0, 'shuffle', 'partner');
for f = fieldnames(p)'
    if isfield(P, f{1}), P.(f{1}) = p.(f{1}); end
end
P.mode    = lower(char(P.mode));
P.shuffle = lower(char(P.shuffle));
if ~ismember(P.mode, {'objects', 'intensity'})
    error('colocContactZone:mode', 'params.mode must be ''objects'' or ''intensity''.');
end
if ~ismember(P.shuffle, {'partner', 'test'})
    error('colocContactZone:shuffle', 'params.shuffle must be ''partner'' or ''test''.');
end
if strcmp(P.mode, 'objects') && (isempty(P.testChannel) || isnan(P.testChannel))
    error('colocContactZone:testChannel', 'Objects mode needs params.testChannel.');
end
if P.distanceA < 0 || P.distanceB < 0
    error('colocContactZone:distance', 'Distances must be >= 0.');
end
P.distanceA = double(P.distanceA);
P.distanceB = double(P.distanceB);
P.nShuffles = double(P.nShuffles);
P.testChannel = double(P.testChannel);
end

function D = distMap(M, cal)
% distance (microns) from every pixel to the nearest true pixel of M; Inf
% everywhere when M is empty (so "near" is never true)
if any(M(:))
    D = double(bwdist(M)) * cal;
else
    D = inf(size(M));
end
end

function M = maskFrom(cp, sz)
M = false(sz);
M(cat(1, cp{:}, zeros(0,1))) = true;
end

function S = slot(X, ch, iZ, iT)
S = table.empty;
if ~iscell(X) || isempty(X) || ch > size(X, 1)
    return
end
S = X{ch, min(iZ, size(X,2)), min(iT, size(X,3))};
end

function S = cellRows(S, g)
if istable(S) && ~isempty(S) && ismember('cellID', S.Properties.VariableNames)
    S = S(S.cellID == g, :);
else
    S = table.empty;
end
end

function p = pixLists(S)
if istable(S) && ~isempty(S)
    p = cellfun(@(x) double(x(:)), S.organellePixelIdxList, 'UniformOutput', false);
    p = p(~cellfun(@isempty, p));
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
id = id(:);
end

function m = safeMean(x)
if isempty(x), m = NaN; else, m = mean(x); end
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
