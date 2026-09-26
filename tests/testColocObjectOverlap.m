classdef testColocObjectOverlap < matlab.unittest.TestCase
%TESTCOLOCOBJECTOVERLAP  Unit tests for colocObjectOverlap (object-based
% colocalisation between two segmented channels).
%
% USAGE
%   results = runtests('tests/testColocObjectOverlap');
%
% COVERAGE
%   Known pair geometry (overlap area, fracOfA/B, Jaccard, per-object
%     rows, nearest distance and contact of a lone object)          — 1 test
%   One A object with two B partners                                 — 1 test
%   k-nearest partners: ranks, edge/centre distances, symmetry,
%     rank-1 = nearestPartnerDistance, kNearest=0 skips              — 1 test
%   Shuffle test: associated populations -> pMoreOverlap/Contact small;
%     independent populations -> false-positive rate controlled      — 2 tests
%   Empty partner population, same-channel error, per-cell
%     separation (no cross-cell pairs), reproducibility              — 4 tests
%
% REQUIREMENTS: OrganelleDistances_sandbox/src (shuffleObjectsInWindow),
% added automatically if found next to this repo.

    properties (Constant)
        Sz  = [100 100]
        Cal = 0.1
    end

    methods (TestClassSetup)
        function addPaths(tc) %#ok<MANU>
            root = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(root, 'src'));
            sib = fullfile(root, '..', 'OrganelleDistances_sandbox', 'src');
            if isfolder(sib), addpath(sib); end
        end
    end

    methods (Static, Access = private)
        function idx = rect(rows, cols)
            [R, C] = ndgrid(rows, cols);
            idx = sub2ind(testColocObjectOverlap.Sz, R(:), C(:));
        end

        function T = stats(pix, cellIDs)
            n = numel(pix);
            if nargin < 2, cellIDs = ones(n,1); end
            T = table(pix(:), (1:n)', cellIDs(:), ...
                'VariableNames', {'organellePixelIdxList','organelleID','cellID'});
        end

        function [O, P, S] = runPair(pixA, pixB, nSh, cellID, cellA, cellB, opts)
            sz = testColocObjectOverlap.Sz;
            if nargin < 4 || isempty(cellID), cellID = ones(sz); end
            if nargin < 5, cellA = ones(numel(pixA),1); cellB = ones(numel(pixB),1); end
            if nargin < 7, opts = struct(); end
            ms = cell(2,1,1);
            ms{1} = testColocObjectOverlap.stats(pixA, cellA);
            ms{2} = testColocObjectOverlap.stats(pixB, cellB);
            [Oc, Pc, Sc] = colocObjectOverlap(ms, cellID, 1, 2, ...
                testColocObjectOverlap.Cal, 0.2, nSh, 'f', opts);
            O = Oc{1}; P = Pc{1}; S = Sc{1};
        end

        function [pix, occ] = discs(n, radius, occ, centres, spread)
            % n non-overlapping discs; uniform if centres empty, else near centres
            sz = testColocObjectOverlap.Sz;
            [dR, dC] = ndgrid(-radius:radius);
            in = dR.^2 + dC.^2 <= radius^2; dR = dR(in); dC = dC(in);
            pix = cell(n,1); k = 0;
            while k < n
                if isempty(centres)
                    c = randi([radius+2, 99-radius], 1, 2);
                else
                    c = round(centres(k+1,:) + spread*randn(1,2));
                end
                r = c(1) + dR; cc = c(2) + dC;
                if any(r < 1 | r > sz(1) | cc < 1 | cc > sz(2)), continue, end
                idx = sub2ind(sz, r, cc);
                if any(occ(idx)), continue, end
                occ(idx) = true; k = k + 1; pix{k} = idx;
            end
        end
    end

    methods (Test)
        function testKnownPairGeometry(tc)
            A = {testColocObjectOverlap.rect(21:30, 21:30); ...   % overlaps B1 by 10x5
                 testColocObjectOverlap.rect(21:30, 61:70)};      % lone, 5 px from B2
            B = {testColocObjectOverlap.rect(21:30, 26:35); ...
                 testColocObjectOverlap.rect(36:45, 61:70)};
            [O, P, S] = testColocObjectOverlap.runPair(A, B, 0);
            cal2 = testColocObjectOverlap.Cal^2;

            tc.verifyEqual(height(P), 1);
            tc.verifyEqual([P.organelleIDA P.organelleIDB], [1 1]);
            tc.verifyEqual(P.overlapArea, 50*cal2, 'AbsTol', 1e-12);
            tc.verifyEqual([P.fracOfA P.fracOfB], [0.5 0.5], 'AbsTol', 1e-12);
            tc.verifyEqual(P.jaccard, 50/150, 'AbsTol', 1e-12);

            OA = O(O.channel == 1, :);
            tc.verifyEqual(OA.nPartners, [1; 0]);
            tc.verifyEqual(OA.colocalised, [1; 0]);
            tc.verifyEqual(OA.nearestPartnerDistance(1), 0);
            tc.verifyEqual(OA.nearestPartnerDistance(2), 0.6, 'AbsTol', 1e-9);  % rows 30->36
            tc.verifyEqual(OA.inContact, [1; 0]);                              % d = 0.2 um

            tc.verifyEqual(S.nPairs, 1);
            tc.verifyEqual(S.fracAOverlapping, 0.5);
            tc.verifyEqual(S.objM1, 50/200, 'AbsTol', 1e-12);
            tc.verifyTrue(isnan(S.pMoreOverlap));            % nShuffles = 0
        end

        function testKNearestPartners(tc)
            % A1 overlaps B1; B2 is 5 px below A1; B3 far away.
            A = {testColocObjectOverlap.rect(21:30, 21:30)};
            B = {testColocObjectOverlap.rect(21:30, 26:35); ...     % overlaps A1
                 testColocObjectOverlap.rect(36:40, 21:30); ...     % rows 30->36: 6 px
                 testColocObjectOverlap.rect(81:85, 81:85)};        % far
            ms = {testColocObjectOverlap.stats(A); testColocObjectOverlap.stats(B)};
            [Oc, ~, ~, Nc] = colocObjectOverlap(ms, ones(testColocObjectOverlap.Sz), 1, 2, ...
                testColocObjectOverlap.Cal, 0.2, 0, 'f', struct('kNearest', 2));
            N = Nc{1}; O = Oc{1};
            NA = N(N.channel == 1, :);
            tc.verifyEqual(NA.rank, [1; 2]);                         % k = 2 of 3 partners
            tc.verifyEqual(NA.partnerID, [1; 2]);
            tc.verifyEqual(NA.edgeDistance, [0; 0.6], 'AbsTol', 1e-9);
            tc.verifyEqual(NA.centreDistance(1), 0.5, 'AbsTol', 1e-9);  % centroids 5 px apart
            % rank-1 edge distance matches objectOut's nearestPartnerDistance
            OA = O(O.channel == 1, :);
            tc.verifyEqual(NA.edgeDistance(1), OA.nearestPartnerDistance(1), 'AbsTol', 1e-9);
            % reverse direction: each B lists its single A partner (only 1 A)
            NB = N(N.channel == 2, :);
            tc.verifyEqual(height(NB), 3);
            tc.verifyEqual(sortrows(NB(:, {'organelleID','edgeDistance'})).edgeDistance(1:2), [0; 0.6], 'AbsTol', 1e-9);
            % symmetric: A1->B2 distance equals B2->A1 distance
            tc.verifyEqual(NB.edgeDistance(NB.organelleID == 2), NA.edgeDistance(2), 'AbsTol', 1e-9);
            % kNearest = 0 skips the table
            [~, ~, ~, N0] = colocObjectOverlap(ms, ones(testColocObjectOverlap.Sz), 1, 2, 0.1, 0.2, 0, 'f', struct('kNearest', 0));
            tc.verifyTrue(isempty(N0{1}));
        end

        function testTwoPartners(tc)
            A = {testColocObjectOverlap.rect(41:60, 41:60)};
            B = {testColocObjectOverlap.rect(41:45, 38:45); testColocObjectOverlap.rect(56:60, 56:63)};
            [O, P] = testColocObjectOverlap.runPair(A, B, 0);
            OA = O(O.channel == 1, :);
            tc.verifyEqual(OA.nPartners, 2);
            tc.verifyEqual(OA.overlapArea, sum(P.overlapArea), 'AbsTol', 1e-12);
            tc.verifyEqual(height(P), 2);
        end

        function testAssociatedPopulationsSignificant(tc)
            rng(1);
            occA = false(testColocObjectOverlap.Sz);
            [A, occA] = testColocObjectOverlap.discs(15, 3, occA, [], 0); %#ok<ASGLU>
            cA = cellfun(@(x) mean(testColocObjectOverlap.sub(x),1), A, 'UniformOutput', false);
            [B] = testColocObjectOverlap.discs(15, 2, false(testColocObjectOverlap.Sz), cat(1, cA{:}), 1.5);
            [~, ~, S] = testColocObjectOverlap.runPair(A, B, 99, [], ones(15,1), ones(15,1), struct('rngSeed', 2));
            tc.verifyLessThanOrEqual(S.pMoreOverlap, 0.02);
            tc.verifyLessThanOrEqual(S.pMoreContact, 0.02);
            tc.verifyGreaterThan(S.objM1, S.objM1NullMean);
        end

        function testIndependentPopulationsNotFlagged(tc)
            rng(3);
            nRep = 40; rej = zeros(nRep, 2);
            for rep = 1:nRep
                [A, occ] = testColocObjectOverlap.discs(12, 3, false(testColocObjectOverlap.Sz), [], 0); %#ok<ASGLU>
                B = testColocObjectOverlap.discs(12, 2, false(testColocObjectOverlap.Sz), [], 0);
                [~, ~, S] = testColocObjectOverlap.runPair(A, B, 39);
                rej(rep,:) = [S.pMoreOverlap <= 0.05, S.pMoreContact <= 0.05];
            end
            tc.verifyLessThanOrEqual(mean(rej(:,1)), 0.15);
            tc.verifyLessThanOrEqual(mean(rej(:,2)), 0.15);
        end

        function testEmptyPartnerPopulation(tc)
            A = {testColocObjectOverlap.rect(21:30, 21:30)};
            sz = testColocObjectOverlap.Sz;
            ms = {testColocObjectOverlap.stats(A); testColocObjectOverlap.stats(cell(0,1))};
            [Oc, Pc, Sc] = colocObjectOverlap(ms, ones(sz), 1, 2, 0.1, 0.2, 19, 'f');
            tc.verifyEqual(height(Pc{1}), 0);
            tc.verifyEqual(Sc{1}.nB, 0);
            tc.verifyTrue(isnan(Sc{1}.pMoreOverlap));
            tc.verifyTrue(isnan(Oc{1}.nearestPartnerDistance(1)));
        end

        function testSameChannelErrors(tc)
            ms = {testColocObjectOverlap.stats({testColocObjectOverlap.rect(1:3,1:3)})};
            tc.verifyError(@() colocObjectOverlap(ms, ones(testColocObjectOverlap.Sz), 1, 1, 0.1, 0.2, 0, 'f'), ...
                'colocObjectOverlap:sameChannel');
        end

        function testPerCellSeparation(tc)
            cellID = ones(testColocObjectOverlap.Sz); cellID(:, 51:end) = 2;
            % A in cell 1 touches the boundary, B in cell 2 overlapping pixels
            % would be impossible (disjoint), so check cross-cell proximity is ignored:
            A = {testColocObjectOverlap.rect(21:30, 41:50)};              % cell 1
            B = {testColocObjectOverlap.rect(21:30, 51:60)};              % cell 2, adjacent
            [O, P, S] = testColocObjectOverlap.runPair(A, B, 0, cellID, 1, 2);
            tc.verifyEqual(height(P), 0);
            tc.verifyEqual(sort(S.cellID), [1; 2]);
            OA = O(O.channel == 1, :);
            tc.verifyTrue(isnan(OA.nearestPartnerDistance), 'partner in another cell must not count');
        end

        function testReproducibleWithSeed(tc)
            rng(4);
            A = testColocObjectOverlap.discs(10, 3, false(testColocObjectOverlap.Sz), [], 0);
            B = testColocObjectOverlap.discs(10, 2, false(testColocObjectOverlap.Sz), [], 0);
            [~, ~, S1] = testColocObjectOverlap.runPair(A, B, 29, [], ones(10,1), ones(10,1), struct('rngSeed', 5));
            [~, ~, S2] = testColocObjectOverlap.runPair(A, B, 29, [], ones(10,1), ones(10,1), struct('rngSeed', 5));
            tc.verifyEqual(S1, S2);
        end
    end

    methods (Static)
        function rc = sub(idx)
            [r, c] = ind2sub(testColocObjectOverlap.Sz, idx);
            rc = [r c];
        end
    end
end
