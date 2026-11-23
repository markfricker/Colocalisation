classdef testColocContactZone < matlab.unittest.TestCase
%TESTCOLOCCONTACTZONE  Unit tests for colocContactZone (three-way
% contact-zone test: test population / intensity vs the zone where A and B
% meet).
%
% USAGE
%   results = runtests('tests/testColocContactZone');
%
% COVERAGE
%   Zone geometry = bwdist rule; distance 0 = overlap of A and B     — 2 tests
%   Per-object fracInZone / inContactZone / distanceToA/B            — 1 test
%   Objects mode: puncta at the junction -> pMore small;
%     independent puncta -> false-positive rate controlled           — 2 tests
%   Intensity mode: junction signal -> pMore small; signal that only
%     follows A -> enriched vs the cell but NOT significant with the
%     partner shuffle                                                — 2 tests
%   Shuffle 'test' (objects); no B in the cell; per-cell separation;
%     argument errors; label columns + reproducibility               — 5 tests
%
% REQUIREMENTS: OrganelleDistances_sandbox/src (shuffleObjectsInWindow),
% added automatically if found next to this repo.

    properties (Constant)
        Sz  = [120 120]
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
            idx = sub2ind(testColocContactZone.Sz, R(:), C(:));
        end

        function S = stats(pix, cellIDs)
            n = numel(pix);
            if nargin < 2, cellIDs = ones(n,1); end
            T = table(pix(:), (1:n)', cellIDs(:), ...
                'VariableNames', {'organellePixelIdxList','organelleID','cellID'});
            S = {T};
        end

        function [pixA, pixB, Z] = network(dA, dB)
            % A = 4 horizontal "tubules"; B = 8 squares beside them (gap 2 px)
            sz = testColocContactZone.Sz;
            pixA = {};
            for r = [20 45 70 95]
                pixA{end+1} = testColocContactZone.rect(r:r+1, 10:110); %#ok<AGROW>
            end
            pixB = {};
            for r = [20 45 70 95]
                for c0 = [25 75]
                    pixB{end+1} = testColocContactZone.rect(r+4:r+9, c0:c0+5); %#ok<AGROW>
                end
            end
            Am = false(sz); Am(cat(1, pixA{:})) = true;
            Bm = false(sz); Bm(cat(1, pixB{:})) = true;
            cal = testColocContactZone.Cal;
            Z = double(bwdist(Am))*cal <= dA & double(bwdist(Bm))*cal <= dB;
        end

        function [C, O, Zc] = runZone(pixA, pixB, test, P, cellID, opts)
            sz = testColocContactZone.Sz;
            if nargin < 5 || isempty(cellID), cellID = ones(sz); end
            if nargin < 6, opts = struct(); end
            SA = testColocContactZone.stats(pixA);
            SB = testColocContactZone.stats(pixB);
            if iscell(test), test = testColocContactZone.stats(test); end
            if ~isfield(P, 'testChannel'), P.testChannel = 1; end
            [Cc, Oc, Zc] = colocContactZone(SA, 1, SB, 1, test, cellID, ...
                testColocContactZone.Cal, P, 'f', opts);
            C = Cc{1}; O = Oc{1};
        end

        function pix = puncta(centres)
            % 3x3 puncta at the given [row col] centres
            pix = cell(size(centres,1), 1);
            for k = 1:size(centres,1)
                pix{k} = testColocContactZone.rect(centres(k,1)+(-1:1), centres(k,2)+(-1:1));
            end
        end
    end

    methods (Test)
        function zoneMatchesDistanceRule(tc)
            [pixA, pixB, Z] = testColocContactZone.network(0.3, 0.3);
            P = struct('distanceA', 0.3, 'distanceB', 0.3, 'mode', 'objects', 'nShuffles', 0);
            [C, ~, Zc] = testColocContactZone.runZone(pixA, pixB, {}, P);
            tc.verifyEqual(Zc == 4, Z);
            tc.verifyEqual(C.zoneArea, nnz(Z) * 0.01, 'AbsTol', 1e-9);
            tc.verifyEqual(C.cellArea, prod(testColocContactZone.Sz) * 0.01, 'AbsTol', 1e-9);
            tc.verifyEqual(C.zoneFraction, nnz(Z) / prod(testColocContactZone.Sz), 'AbsTol', 1e-12);
            % classes are disjoint and cover the cell
            tc.verifyTrue(all(Zc(:) >= 1 & Zc(:) <= 4));
            tc.verifyEqual(C.nA, 4);
            tc.verifyEqual(C.nB, 8);
        end

        function zeroDistanceIsOverlap(tc)
            pixA = {testColocContactZone.rect(50:51, 10:110)};
            pixB = {testColocContactZone.rect(45:55, 30:35), testColocContactZone.rect(80:85, 30:35)};
            P = struct('distanceA', 0, 'distanceB', 0, 'mode', 'objects', 'nShuffles', 0);
            [C, ~, Zc] = testColocContactZone.runZone(pixA, pixB, {}, P);
            tc.verifyEqual(nnz(Zc == 4), 2 * 6);   % rows 50:51 x cols 30:35
            tc.verifyEqual(C.zoneArea, 12 * 0.01, 'AbsTol', 1e-9);
        end

        function perObjectRows(tc)
            pixA = {testColocContactZone.rect(50:51, 10:110)};
            pixB = {testColocContactZone.rect(54:59, 40:45)};
            % punctum 1 in the gap between A and B; punctum 2 far away
            T = testColocContactZone.puncta([53 42; 100 100]);
            P = struct('distanceA', 0.3, 'distanceB', 0.3, 'mode', 'objects', 'nShuffles', 0);
            [C, O] = testColocContactZone.runZone(pixA, pixB, T, P);
            tc.verifyEqual(height(O), 2);
            tc.verifyEqual(O.inContactZone, [1; 0]);
            tc.verifyGreaterThan(O.fracInZone(1), 0);
            tc.verifyEqual(O.fracInZone(2), 0);
            % punctum 1 rows 52:54 -> touches A (row 51) at distance 1 px, overlaps B (row 54)
            tc.verifyEqual(O.distanceToA(1), 0.1, 'AbsTol', 1e-9);
            tc.verifyEqual(O.distanceToB(1), 0, 'AbsTol', 1e-9);
            % punctum 2 rows 99:101 -> nearest A row 51: 48 px
            tc.verifyEqual(O.distanceToA(2), 4.8, 'AbsTol', 1e-9);
            tc.verifyEqual(C.nTest, 2);
            tc.verifyEqual(C.fracTestInZone, 0.5);
            tc.verifyEqual(C.statistic, 0.5);
            tc.verifyTrue(isnan(C.pMore));
            tc.verifyTrue(isnan(C.meanZone));   % intensity columns NaN in objects mode
        end

        function objectsJunctionEnriched(tc)
            [pixA, pixB] = testColocContactZone.network(0.3, 0.3);
            % one punctum in each A-B gap (B squares at rows r+4.., A at r:r+1)
            cen = [];
            for r = [20 45 70 95]
                for c0 = [25 75]
                    cen(end+1, :) = [r+3, c0+2]; %#ok<AGROW>
                end
            end
            T = testColocContactZone.puncta(cen);
            P = struct('distanceA', 0.3, 'distanceB', 0.3, 'mode', 'objects', 'nShuffles', 99);
            C = testColocContactZone.runZone(pixA, pixB, T, P, [], struct('rngSeed', 1));
            tc.verifyEqual(C.fracTestInZone, 1);
            tc.verifyLessThan(C.pMore, 0.05);
            tc.verifyGreaterThan(C.enrichment, 2);
            tc.verifyGreaterThan(C.areaEnrichment, 2);
        end

        function objectsIndependentFalsePositiveRate(tc)
            [pixA, pixB] = testColocContactZone.network(0.3, 0.3);
            P = struct('distanceA', 0.3, 'distanceB', 0.3, 'mode', 'objects', 'nShuffles', 49);
            rng(7);
            nRep = 20;
            p = zeros(nRep, 1);
            for k = 1:nRep
                % puncta placed uniformly, then B shuffled for the null:
                % under H0 (B position independent of the puncta) p ~ U(0,1)
                pixBr = shuffleObjectsInWindow(pixB(:), true(testColocContactZone.Sz));
                cen = randi([3 118], 12, 2);
                T = testColocContactZone.puncta(cen);
                C = testColocContactZone.runZone(pixA, pixBr, T, P);
                p(k) = C.pMore;
            end
            tc.verifyLessThanOrEqual(mean(p(~isnan(p)) < 0.05), 0.2);
        end

        function intensityJunctionSignificant(tc)
            [pixA, pixB, Z] = testColocContactZone.network(0.3, 0.3);
            I = ones(testColocContactZone.Sz) + 4 * Z;      % bright only at the junctions
            P = struct('distanceA', 0.3, 'distanceB', 0.3, 'mode', 'intensity', 'nShuffles', 99);
            [C, O] = testColocContactZone.runZone(pixA, pixB, I, P, [], struct('rngSeed', 2));
            tc.verifyEqual(C.meanZone, 5, 'AbsTol', 1e-9);
            tc.verifyEqual(C.meanOutsideZone, 1, 'AbsTol', 1e-9);
            tc.verifyEqual(C.enrichment, 5, 'AbsTol', 1e-9);
            tc.verifyLessThan(C.pMore, 0.05);
            tc.verifyTrue(isnan(C.nTest));
            tc.verifyEqual(height(O), 0);
        end

        function intensityFollowingANotSignificant(tc)
            % signal on and near A everywhere (an ER membrane protein): the
            % zone is brighter than the rest of the cell, but shuffled B
            % zones are next to A too, so the partner shuffle is not fooled
            [pixA, pixB] = testColocContactZone.network(0.3, 0.3);
            Am = false(testColocContactZone.Sz); Am(cat(1, pixA{:})) = true;
            I = ones(testColocContactZone.Sz) + 4 * (bwdist(Am) <= 3);
            P = struct('distanceA', 0.3, 'distanceB', 0.3, 'mode', 'intensity', 'nShuffles', 99);
            C = testColocContactZone.runZone(pixA, pixB, I, P, [], struct('rngSeed', 3));
            tc.verifyGreaterThan(C.enrichment, 2);       % naive ratio looks enriched
            tc.verifyGreaterThan(C.pMore, 0.05);         % but not vs the partner shuffle
            tc.verifyGreaterThan(C.meanNearAOnly, C.meanNeither);
        end

        function shuffleTestObjects(tc)
            [pixA, pixB] = testColocContactZone.network(0.3, 0.3);
            cen = [];
            for r = [20 45 70 95]
                cen(end+1, :) = [r+3, 27]; %#ok<AGROW>
            end
            T = testColocContactZone.puncta(cen);
            P = struct('distanceA', 0.3, 'distanceB', 0.3, 'mode', 'objects', ...
                'nShuffles', 99, 'shuffle', 'test');
            C = testColocContactZone.runZone(pixA, pixB, T, P, [], struct('rngSeed', 4));
            tc.verifyEqual(C.shuffle, {'test'});
            tc.verifyLessThan(C.pMore, 0.05);
            tc.verifyEqual(C.fallbackFraction, 0);
        end

        function noPartnerInCell(tc)
            pixA = {testColocContactZone.rect(50:51, 10:110)};
            T = testColocContactZone.puncta([53 42]);
            P = struct('distanceA', 0.3, 'distanceB', 0.3, 'mode', 'objects', 'nShuffles', 19);
            [C, O, Zc] = testColocContactZone.runZone(pixA, {}, T, P);
            tc.verifyEqual(C.nB, 0);
            tc.verifyEqual(C.zoneArea, 0);
            tc.verifyTrue(isnan(C.pMore));
            tc.verifyTrue(isnan(O.distanceToB));
            tc.verifyFalse(any(Zc(:) == 4 | Zc(:) == 3));
        end

        function perCellSeparation(tc)
            sz = testColocContactZone.Sz;
            cellID = ones(sz); cellID(:, 61:end) = 2;
            pixA = {testColocContactZone.rect(50:51, 10:50), testColocContactZone.rect(50:51, 70:110)};
            pixB = {testColocContactZone.rect(54:59, 30:35)};      % only in cell 1
            SA = {table(pixA(:), [1;2], [1;2], 'VariableNames', {'organellePixelIdxList','organelleID','cellID'})};
            SB = {table(pixB(:), 1, 1, 'VariableNames', {'organellePixelIdxList','organelleID','cellID'})};
            P = struct('distanceA', 0.3, 'distanceB', 0.3, 'mode', 'intensity', 'nShuffles', 0);
            [Cc, ~, Zc] = colocContactZone(SA, 1, SB, 1, ones(sz), cellID, 0.1, P, 'f');
            C = Cc{1};
            tc.verifyEqual(C.cellID, [1; 2]);
            tc.verifyEqual(C.nA, [1; 1]);
            tc.verifyEqual(C.nB, [1; 0]);
            tc.verifyGreaterThan(C.zoneArea(1), 0);
            tc.verifyEqual(C.zoneArea(2), 0);
            tc.verifyFalse(any(any(Zc(:, 61:end) == 4)));
        end

        function argumentErrors(tc)
            S = {table({1}, 1, 1, 'VariableNames', {'organellePixelIdxList','organelleID','cellID'})};
            I = ones(testColocContactZone.Sz);
            cid = ones(testColocContactZone.Sz);
            tc.verifyError(@() colocContactZone(S, 1, S, 1, I, cid, 0.1, ...
                struct('mode', 'intensity', 'shuffle', 'test'), 'f'), ...
                'colocContactZone:shuffleTestIntensity');
            tc.verifyError(@() colocContactZone(S, 1, S, 1, S, cid, 0.1, ...
                struct('mode', 'objects'), 'f'), 'colocContactZone:testChannel');
            tc.verifyError(@() colocContactZone(S, 1, S, 1, I, cid, 0.1, ...
                struct('mode', 'pixels'), 'f'), 'colocContactZone:mode');
            tc.verifyError(@() colocContactZone(S, 1, S, 1, I, cid, 0.1, ...
                struct('mode', 'intensity', 'distanceA', -1), 'f'), 'colocContactZone:distance');
        end

        function labelsAndReproducibility(tc)
            [pixA, pixB, Z] = testColocContactZone.network(0.3, 0.3);
            I = ones(testColocContactZone.Sz) + Z;
            P = struct('distanceA', 0.3, 'distanceB', 0.3, 'mode', 'intensity', 'nShuffles', 19);
            opts = struct('rngSeed', 11, 'labelColumns', struct('populationA', 'ER whole', 'run', 3));
            C1 = testColocContactZone.runZone(pixA, pixB, I, P, [], opts);
            C2 = testColocContactZone.runZone(pixA, pixB, I, P, [], opts);
            tc.verifyEqual(C1, C2);
            tc.verifyEqual(C1.populationA, {'ER whole'});
            tc.verifyEqual(C1.run, 3);
        end
    end
end
