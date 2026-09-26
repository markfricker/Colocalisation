classdef testColocCisternaeObjects < matlab.unittest.TestCase
%TESTCOLOCCISTERNAEOBJECTS  Unit tests for colocCisternaeObjects and for
% colocObjectOverlap with a separate population-B array (opts.statsB):
% ER cisternae (A) vs organelles (B).
%
% COVERAGE
%   column mapping, all/ordinary/streams split, threshold 0 = no streams,
%   missing speed column, invalid class                                — 4 tests
%   colocObjectOverlap with opts.statsB (same channel index allowed),
%   cell present in only one population                                — 2 tests

    properties (Constant)
        Sz = [100 100]
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
            idx = sub2ind(testColocCisternaeObjects.Sz, R(:), C(:));
        end
        function S = cisternae(withSpeed)
            r = @testColocCisternaeObjects.rect;
            S = table({r(11:20,11:20); r(11:20,41:50); r(61:70,61:70)}, [1;2;3], [1;1;1], ...
                'VariableNames', {'cisternaePixelIdxList','cisternaeID','cellID'});
            if withSpeed
                S.cisternaeSpeedMax = single([0.1; 0.8; 1.5]);
            end
        end
    end

    methods (Test)
        function testMappingAndSplit(tc)
            cs = {testColocCisternaeObjects.cisternae(true)};
            [all_, iA] = colocCisternaeObjects(cs, 'all', 0.5);
            tc.verifyEqual(all_{1}.organelleID, [1;2;3]);
            tc.verifyEqual(all_{1}.isStream, [false; true; true]);
            tc.verifyEqual([iA.nTotal iA.nKept], [3 3]);
            ord = colocCisternaeObjects(cs, 'ordinary', 0.5);
            tc.verifyEqual(ord{1}.organelleID, 1);
            str = colocCisternaeObjects(cs, 'streams', 0.5);
            tc.verifyEqual(str{1}.organelleID, [2;3]);
            tc.verifyEqual(str{1}.organellePixelIdxList{1}, cs{1}.cisternaePixelIdxList{2});
        end

        function testThresholdZeroMeansNoStreams(tc)
            cs = {testColocCisternaeObjects.cisternae(true)};
            ord = colocCisternaeObjects(cs, 'ordinary', 0);
            tc.verifyEqual(height(ord{1}), 3);
            [str, info] = colocCisternaeObjects(cs, 'streams', 0);
            tc.verifyEqual(height(str{1}), 0);
            tc.verifyNotEmpty(info.message);
        end

        function testMissingSpeed(tc)
            cs = {testColocCisternaeObjects.cisternae(false)};
            [str, info] = colocCisternaeObjects(cs, 'streams', 0.5);
            tc.verifyFalse(info.hasSpeed);
            tc.verifyEqual(height(str{1}), 0);
            tc.verifyTrue(contains(info.message, 'optical flow'));
            allc = colocCisternaeObjects(cs, 'all', 0.5);
            tc.verifyEqual(height(allc{1}), 3);
        end

        function testInvalidClass(tc)
            tc.verifyError(@() colocCisternaeObjects({}, 'tubules', 0), 'colocCisternaeObjects:class');
        end

        function testOverlapWithSeparateStatsB(tc)
            r = @testColocCisternaeObjects.rect;
            cis = colocCisternaeObjects({testColocCisternaeObjects.cisternae(true)}, 'all', 0.5);
            org = {table({r(15:17,15:17); r(80:82,80:82)}, [1;2], [1;1], ...
                'VariableNames', {'organellePixelIdxList','organelleID','cellID'})};
            % ER channel 1 vs organelle "channel 1" of a different array: allowed
            [O, P, S, N] = colocObjectOverlap(cis, ones(testColocCisternaeObjects.Sz), 1, 1, ...
                0.1, 0.2, 19, 'f', struct('statsB', {org}, 'rngSeed', 1));
            tc.verifyEqual(height(P{1}), 1);                         % organelle 1 inside cisterna 1
            tc.verifyEqual([P{1}.organelleIDA P{1}.organelleIDB], [1 1]);
            tc.verifyEqual(P{1}.fracOfB, 1);                         % whole organelle inside
            tc.verifyEqual(S{1}.nA, 3);
            tc.verifyEqual(S{1}.nB, 2);
            tc.verifyFalse(isnan(S{1}.pMoreOverlap));
            tc.verifyEqual(height(O{1}), 5);
            tc.verifyGreaterThan(height(N{1}), 0);
        end

        function testSavedFlagTakesPrecedence(tc)
            cs = {testColocCisternaeObjects.cisternae(true)};
            cs{1}.cisternaeIsStream = [1; 0; 0];         % deliberately disagrees with speed
            str = colocCisternaeObjects(cs, 'streams', 0.5);
            tc.verifyEqual(str{1}.organelleID, 1);        % flag wins over speed >= 0.5
            ord = colocCisternaeObjects(cs, 'ordinary', 0.5);
            tc.verifyEqual(ord{1}.organelleID, [2; 3]);
        end

        function testIsStreamCarriedAndLabelsStamped(tc)
            r = @testColocCisternaeObjects.rect;
            cellID = ones(testColocCisternaeObjects.Sz); cellID(:, 51:end) = 2;
            cs = testColocCisternaeObjects.cisternae(true);
            cs.cellID = [1; 1; 2];                       % cisterna 3 in cell 2
            cis = colocCisternaeObjects({cs}, 'all', 0.5);
            org = {table({r(15:17,15:17); r(64:66,64:66); r(90:92,90:92)}, [1;2;3], [1;2;2], ...
                'VariableNames', {'organellePixelIdxList','organelleID','cellID'})};
            lab = struct('cisternaeClass', 'all', 'streamsThreshold', 0.5);
            [O, P, S, N] = colocObjectOverlap(cis, cellID, 1, 1, 0.1, 0.2, 0, 'f', ...
                struct('statsB', {org}, 'labelColumns', lab));
            O = O{1}; P = P{1}; S = S{1}; N = N{1};
            OA = O(O.partnerChannel == 1 & ~isnan(O.isStream), :);
            tc.verifyEqual(sortrows([OA.organelleID OA.isStream]), [1 0; 2 1; 3 1]);
            tc.verifyEqual(nnz(isnan(O.isStream)), 3);                 % the 3 organelle rows
            for T = {O, P, S, N}
                tc.verifyTrue(all(strcmp(T{1}.cisternaeClass, 'all')));
                tc.verifyTrue(all(T{1}.streamsThreshold == 0.5));
            end
            tc.verifyEqual(height(S), 2);                % both cells, one plane: concatenated fine
        end

        function testCellInOnlyOnePopulation(tc)
            r = @testColocCisternaeObjects.rect;
            cellID = ones(testColocCisternaeObjects.Sz); cellID(:, 51:end) = 2;
            A = {table({r(11:20,11:20)}, 1, 1, 'VariableNames', {'organellePixelIdxList','organelleID','cellID'})};
            B = {table({r(11:13,61:63)}, 1, 2, 'VariableNames', {'organellePixelIdxList','organelleID','cellID'})};
            [~, ~, S] = colocObjectOverlap(A, cellID, 1, 1, 0.1, 0.2, 0, 'f', struct('statsB', {B}));
            tc.verifyEqual(sortrows(S{1}(:, {'cellID','nA','nB'})).Variables, [1 1 0; 2 0 1]);
        end
    end
end
