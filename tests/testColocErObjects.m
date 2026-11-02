classdef testColocErObjects < matlab.unittest.TestCase
%TESTCOLOCEROBJECTS  Unit tests for colocErObjects (whole ER, tubules +
% cisternae, as one object population) and colocErOverlapBreakdown (the
% tubule / cisterna / stream split of each organelle's overlap).
%
% COVERAGE
%   ER objects: per-cell connected pieces, cisternae added to the mask   — 1 test
%   breakdown: tubule / cisterna / stream areas and partner counts,
%   per-cell fractions, re-run at a new threshold, no cisternae -> NaN   — 3 tests

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
            idx = sub2ind(testColocErObjects.Sz, R(:), C(:));
        end
        function [erMask, cellID, cs, ms] = scene()
            % cell 1 = columns 1:60, cell 2 = columns 61:100. ER: a 2-px
            % tubule across both cells, a branch up to cisterna 1 (slow);
            % cisterna 2 (fast) is detached and missing from the mask.
            r = @testColocErObjects.rect;
            sz = testColocErObjects.Sz;
            cellID = ones(sz);
            cellID(:, 61:100) = 2;
            erMask = false(sz);
            erMask(r(30:31, 10:90)) = true;
            erMask(r(20:30, 15:16)) = true;
            erMask(r(11:20, 11:20)) = true;
            cs = {table({r(11:20,11:20); r(11:20,41:50)}, [1;2], [1;1], single([0.1; 0.8]), ...
                'VariableNames', {'cisternaePixelIdxList','cisternaeID','cellID','cisternaeSpeedMax'})};
            % organelles (channel 2): 1 on a tubule, 2 on the fast cisterna,
            % 3 on the slow cisterna, 4 on nothing (cell 2)
            T = table({r(28:33,20:25); r(15:24,46:55); r(5:12,12:14); r(70:75,70:75)}, ...
                (1:4)', [1;1;1;2], 'VariableNames', {'organellePixelIdxList','organelleID','cellID'});
            ms = {table.empty; T};
        end
        function [O, S] = runScene(thr, withCisternae)
            [erMask, cellID, cs, ms] = testColocErObjects.scene();
            A = colocErObjects(erMask, cellID, 1, cs);
            opts = struct('kNearest', 0);
            opts.statsB = ms;
            [O, ~, S] = colocObjectOverlap(A, cellID, 1, 2, testColocErObjects.Cal, 0, 0, 'x', opts);
            cisObj = {};
            if withCisternae
                cisObj = colocCisternaeObjects(cs, 'all', thr);
            end
            [O, S] = colocErOverlapBreakdown(O, S, ms, 2, cisObj, 1, cellID, testColocErObjects.Cal);
            O = O{1}; S = S{1};
        end
    end

    methods (Test)
        function testErObjects(tc)
            [erMask, cellID, cs] = testColocErObjects.scene();
            A = colocErObjects(erMask, cellID, 1, cs);
            T = A{1};
            tc.verifyEqual(sort(T.cellID), [1;1;2], 'tubule split at the cell border; detached cisterna is its own object');
            tc.verifyEqual(T.organelleID, (1:3)');
            allPix = cat(1, T.organellePixelIdxList{:});
            tc.verifyTrue(all(ismember(cs{1}.cisternaePixelIdxList{2}, allPix)), 'cisterna missing from the mask is added');
            tc.verifyEqual(numel(allPix), nnz(erMask) + 100);
            A0 = colocErObjects(erMask, cellID, 1);
            tc.verifyEqual(height(A0{1}), 2, 'without cisternae stats the mask alone is used');
        end

        function testBreakdown(tc)
            [O, S] = testColocErObjects.runScene(0.5, true);
            B = O(strcmp(O.population, 'B'), :);
            B = sortrows(B, 'organelleID');
            tol = {'AbsTol', 1e-12};
            tc.verifyEqual(B.overlapArea,          [0.12; 0.30; 0.06; 0], tol{:});
            tc.verifyEqual(B.overlapAreaTubules,   [0.12; 0; 0; 0], tol{:});
            tc.verifyEqual(B.overlapAreaCisternae, [0; 0.30; 0.06; 0], tol{:});
            tc.verifyEqual(B.overlapAreaStream,    [0; 0.30; 0; 0], tol{:});
            tc.verifyEqual(B.nPartnersCisternae,   [0; 1; 1; 0]);
            tc.verifyEqual(B.nPartnersStream,      [0; 1; 0; 0]);
            tc.verifyEqual(B.maxOverlapPartnerCisternaeSpeedMax, [NaN; 0.8; 0.1; NaN], 'AbsTol', 1e-6);
            A = O(strcmp(O.population, 'A'), :);
            tc.verifyTrue(all(isnan(A.overlapAreaTubules)), 'ER rows get NaN');
            S = sortrows(S, 'cellID');
            tc.verifyEqual(S.nCisternae, [2; 0]);
            tc.verifyEqual(S.nStreams,   [1; 0]);
            tc.verifyEqual(S.fracBOverlappingTubules,   [1/3; 0], tol{:});
            tc.verifyEqual(S.fracBOverlappingCisternae, [2/3; 0], tol{:});
            tc.verifyEqual(S.fracBOverlappingStream,    [1/3; 0], tol{:});
        end

        function testNewThresholdReapplied(tc)
            [erMask, cellID, cs, ms] = testColocErObjects.scene();
            [O, S] = testColocErObjects.runScene(0.5, true);
            % re-apply to the stored outputs at a higher threshold
            cisObj = colocCisternaeObjects(cs, 'all', 1.0);
            [O2, S2] = colocErOverlapBreakdown({O}, {S}, ms, 2, cisObj, 1, cellID, testColocErObjects.Cal);
            O2 = O2{1}; S2 = S2{1};
            tc.verifyEqual(width(O2), width(O), 'columns overwritten, not duplicated');
            B = sortrows(O2(strcmp(O2.population, 'B'), :), 'organelleID');
            tc.verifyEqual(B.overlapAreaStream, zeros(4,1));
            tc.verifyEqual(B.overlapAreaCisternae, [0; 0.30; 0.06; 0], 'AbsTol', 1e-12, 'threshold-free columns unchanged');
            tc.verifyEqual(sortrows(S2, 'cellID').nStreams, [0; 0]);
            tc.verifyNotEmpty(erMask);
        end

        function testNoCisternaeGivesNaN(tc)
            [O, S] = testColocErObjects.runScene(0.5, false);
            B = O(strcmp(O.population, 'B'), :);
            tc.verifyTrue(all(isnan(B.overlapAreaTubules)));
            tc.verifyTrue(all(isnan(B.overlapAreaStream)));
            tc.verifyTrue(all(isnan(S.fracBOverlappingCisternae)));
        end
    end
end
