classdef testColocObjectOverlap3D < matlab.unittest.TestCase
%TESTCOLOCOBJECTOVERLAP3D  Unit tests for colocObjectOverlap3D (3D
% object-based colocalisation) and colocObjectClassImage3D.
%
% USAGE
%   results = runtests('tests/testColocObjectOverlap3D');
%
% COVERAGE
%   Known box pair: overlap volume, fractions, Jaccard, volume columns  - 1 test
%   Distances: anisotropic Z gap, brute-force check of every distance   - 2 tests
%   Shuffle: associated -> small pMoreOverlap; independent -> not       - 1 test
%   3D shuffle keeps shapes inside the window, no B-B overlap           - 1 test
%   Cells: 2D cell labels extruded through Z, no cross-cell pairs       - 1 test
%   Shared single-channel labels error; class image size and colours    - 2 tests
%
% REQUIREMENTS: OrganelleDistances_sandbox/src (shuffleObjectsInWindow),
% added automatically if found next to this repo.

    properties (Constant)
        Sz  = [40 40 12]
        Sp  = [0.1 0.1 0.3]       % dx dy dz (um), anisotropic like a confocal stack
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
        function L = box(L, ch, id, r, c, z)
            L(r, c, ch, z) = id;
        end

        function T = stats(L, ch, cellIDs)
            ids = unique(L(:,:,ch,:));
            ids = double(ids(ids > 0));
            if nargin < 3, cellIDs = zeros(size(ids)); end
            T = table(cellIDs(:), ids(:), 'VariableNames', {'cellID','objectID'});
        end

        function P = voxUm(V, ids, sp)
            % voxel centres [x y z] in um of the labels ids in V [nY nX 1 nZ]
            [r, c, ~, z] = ind2sub(size(V), find(ismember(V, ids)));
            P = [c * sp(1), r * sp(2), z * sp(3)];
        end

        function [O, P, S, N] = runPair(L, nSh, cellID, opts, cellA, cellB)
            if nargin < 3, cellID = []; end
            if nargin < 4, opts = struct(); end
            ms = cell(2, 1, 1);
            if nargin < 5
                ms{1} = testColocObjectOverlap3D.stats(L, 1);
                ms{2} = testColocObjectOverlap3D.stats(L, 2);
            else
                ms{1} = testColocObjectOverlap3D.stats(L, 1, cellA);
                ms{2} = testColocObjectOverlap3D.stats(L, 2, cellB);
            end
            [Oc, Pc, Sc, Nc] = colocObjectOverlap3D(L, ms, cellID, 1, 2, ...
                testColocObjectOverlap3D.Sp, 0.2, nSh, 'f', opts);
            O = Oc{1}; P = Pc{1}; S = Sc{1}; N = Nc{1};
        end
    end

    methods (Test)
        function knownBoxPair(tc)
            L = zeros([40 40 2 12]);
            L = testColocObjectOverlap3D.box(L, 1, 1, 5:14, 5:14, 3:6);    % 10x10x4 = 400 vox
            L = testColocObjectOverlap3D.box(L, 2, 1, 10:19, 5:14, 5:8);   % 10x10x4, overlap 5x10x2 = 100
            [O, P] = testColocObjectOverlap3D.runPair(L, 0);
            vV = prod(testColocObjectOverlap3D.Sp);
            tc.verifyEqual(height(P), 1);
            tc.verifyEqual(P.overlapVolume, 100 * vV, 'AbsTol', 1e-12);
            tc.verifyEqual(P.fracOfA, 0.25, 'AbsTol', 1e-12);
            tc.verifyEqual(P.jaccard, 100 / 700, 'AbsTol', 1e-12);
            tc.verifyTrue(all(ismember({'volume','overlapVolume','linkZ1','linkZ2'}, O.Properties.VariableNames)));
            tc.verifyFalse(ismember('area', O.Properties.VariableNames));
            tc.verifyEqual(O.volume, [400; 400] * vV, 'AbsTol', 1e-12);
            tc.verifyEqual(O.nearestPartnerDistance, [0; 0]);
            tc.verifyEqual(O.colocalised, [1; 1]);
        end

        function zGapCountsDz(tc)
            % same XY footprint, sections 2-3 and 7-8: gap of 4 sections,
            % centre-to-centre of the closest voxels = 4 * dz
            L = zeros([40 40 2 12]);
            L = testColocObjectOverlap3D.box(L, 1, 1, 10:15, 10:15, 2:3);
            L = testColocObjectOverlap3D.box(L, 2, 1, 10:15, 10:15, 7:8);
            [O, ~, ~, N] = testColocObjectOverlap3D.runPair(L, 0);
            tc.verifyEqual(O.nearestPartnerDistance, [4; 4] * 0.3, 'AbsTol', 1e-9);
            tc.verifyEqual(O.inContact, [0; 0]);
            tc.verifyEqual(O.linkZ1(1), 3);
            tc.verifyEqual(O.linkZ2(1), 7);
            tc.verifyEqual(N.edgeDistance(N.rank == 1), [1.2; 1.2], 'AbsTol', 1e-9);
            tc.verifyEqual(N.centreDistance(1), 5 * 0.3, 'AbsTol', 1e-9);
        end

        function distancesMatchBruteForce(tc)
            rng(3)
            L = zeros([40 40 2 12]);
            for k = 1:4
                r = randi(30); c = randi(30); z = randi(9);
                L = testColocObjectOverlap3D.box(L, 1, k, r:r+3, c:c+2, z:z+1);
            end
            for k = 1:5
                r = randi(30); c = randi(30); z = randi(9);
                L(r:r+2, c:c+3, 2, z:z+2) = k;
            end
            [O, ~, ~, N] = testColocObjectOverlap3D.runPair(L, 0, [], struct('kNearest', 5));
            sp = testColocObjectOverlap3D.Sp;
            xyz = @(ch, id) testColocObjectOverlap3D.voxUm(L(:,:,ch,:), id, sp);
            idsA = unique(L(:,:,1,:)); idsA = idsA(idsA > 0);
            idsB = unique(L(:,:,2,:)); idsB = idsB(idsB > 0);
            allB = testColocObjectOverlap3D.voxUm(L(:,:,2,:), idsB, sp);
            for a = idsA'
                pa = xyz(1, a);
                d = min(sqrt(sum((permute(pa, [1 3 2]) - permute(allB, [3 1 2])).^2, 3)), [], 'all');
                % overlap -> 0
                if any(ismember(find(L(:,:,1,:) == a), find(L(:,:,2,:) > 0))), d = 0; end
                row = O(strcmp(O.population, 'A') & O.organelleID == a, :);
                tc.verifyEqual(row.nearestPartnerDistance, d, 'AbsTol', 1e-9);
                for b = idsB'
                    pb = xyz(2, b);
                    e = min(sqrt(sum((permute(pa, [1 3 2]) - permute(pb, [3 1 2])).^2, 3)), [], 'all');
                    if any(ismember(find(L(:,:,1,:) == a), find(L(:,:,2,:) == b))), e = 0; end
                    nr = N(N.channel == 1 & N.organelleID == a & N.partnerID == b, :);
                    tc.verifyEqual(nr.edgeDistance, e, 'AbsTol', 1e-9);
                end
            end
        end

        function shuffleSeparatesAssociation(tc)
            % B sits inside A (associated) vs B in random places
            L = zeros([40 40 2 12]);
            L2 = L;
            rng(1)
            for k = 1:6
                r = 3 + 6*(k-1); c = 5; z = 2 + mod(k, 3) * 3;
                L = testColocObjectOverlap3D.box(L, 1, k, r:r+3, c:c+3, z:z+2);
                L = testColocObjectOverlap3D.box(L, 2, k, r+1:r+2, c+1:c+2, z+1);
                L2 = testColocObjectOverlap3D.box(L2, 1, k, r:r+3, c:c+3, z:z+2);
                L2(5 + 5*(k-1):6 + 5*(k-1), 30:31, 2, 6) = k;    % far from A
            end
            [~, ~, S]  = testColocObjectOverlap3D.runPair(L,  99, [], struct('rngSeed', 1));
            [~, ~, S2] = testColocObjectOverlap3D.runPair(L2, 99, [], struct('rngSeed', 1));
            tc.verifyLessThan(S.pMoreOverlap, 0.05);
            tc.verifyLessThan(S.pMoreContact, 0.05);
            tc.verifyGreaterThan(S2.pMoreOverlap, 0.5);
            tc.verifyEqual(S.nShuffles, 99);
        end

        function shuffle3DKeepsShapes(tc)
            win = false(30, 30, 10);
            win(5:25, 5:25, 2:9) = true;
            [R, C, Z] = ndgrid(10:12, 10:13, 4:5);
            obj = {sub2ind(size(win), R(:), C(:), Z(:)); sub2ind(size(win), R(:)+5, C(:)+5, Z(:)+2)};
            rng(2)
            [placed, fb] = shuffleObjectsInWindow(obj, win);
            tc.verifyFalse(any(fb));
            tc.verifyTrue(all(win(cat(1, placed{:}))));
            tc.verifyEmpty(intersect(placed{1}, placed{2}));
            for k = 1:2
                [r0, c0, z0] = ind2sub(size(win), obj{k});
                [r1, c1, z1] = ind2sub(size(win), placed{k});
                d = [r1 c1 z1] - [r0 c0 z0];
                tc.verifyEqual(size(unique(d, 'rows'), 1), 1);   % pure translation
            end
        end

        function cellsExtrudedThroughZ(tc)
            L = zeros([40 40 2 12]);
            L = testColocObjectOverlap3D.box(L, 1, 1, 5:8,  5:8, 3:5);    % cell 1
            L = testColocObjectOverlap3D.box(L, 2, 1, 5:8, 22:25, 3:5);   % cell 2, touching in XY? no: far
            L = testColocObjectOverlap3D.box(L, 2, 2, 6:7,  6:7, 4);      % cell 1, inside A1
            cellID = ones(40, 40);  cellID(:, 21:end) = 2;
            [O, P, S] = testColocObjectOverlap3D.runPair(L, 0, cellID, struct(), 1, [2; 1]);
            tc.verifyEqual(height(S), 2);
            tc.verifyEqual(height(P), 1);
            tc.verifyEqual(P.organelleIDB, 2);
            rowB1 = O(strcmp(O.population, 'B') & O.organelleID == 1, :);
            tc.verifyTrue(isnan(rowB1.nearestPartnerDistance));  % no A in cell 2
        end

        function sharedLabelsError(tc)
            L = zeros([40 40 1 12]);
            L(5:8, 5:8, 1, 3:5) = 1;
            ms = {table(0, 1, 'VariableNames', {'cellID','objectID'}); table(0, 1, 'VariableNames', {'cellID','objectID'})};
            tc.verifyError(@() colocObjectOverlap3D(L, ms, [], 1, 2, [0.1 0.3], 0.2, 0, 'f'), ...
                'colocObjectOverlap3D:sharedLabels');
        end

        function classImage(tc)
            L = zeros([40 40 2 12]);
            L = testColocObjectOverlap3D.box(L, 1, 1, 5:14, 5:14, 3:6);
            L = testColocObjectOverlap3D.box(L, 2, 1, 10:19, 5:14, 5:8);
            L = testColocObjectOverlap3D.box(L, 2, 2, 30:32, 30:32, 10);   % lone B
            O = colocObjectOverlap3D(L, {testColocObjectOverlap3D.stats(L,1); testColocObjectOverlap3D.stats(L,2)}, ...
                [], 1, 2, testColocObjectOverlap3D.Sp, 0.2, 0, 'f');
            rgb = colocObjectClassImage3D(O, L, 1, 2);
            tc.verifyEqual(size(rgb), [40 40 1 12 1 3]);
            cmap = uint8(round(255 * colocObjectClassColours()));
            tc.verifyEqual(squeeze(rgb(11, 10, 1, 5, 1, :))', cmap(7, :));   % overlap voxel
            tc.verifyEqual(squeeze(rgb(31, 31, 1, 10, 1, :))', cmap(3, :));  % lone B: no contact
            tc.verifyEqual(squeeze(rgb(5, 9, 1, 3, 1, :))', cmap(1, :));     % A outline
        end
    end
end
