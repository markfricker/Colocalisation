classdef testColocObjectVisualisation < matlab.unittest.TestCase
%TESTCOLOCOBJECTVISUALISATION  Tests for the object-colocalisation
% visualisation helpers: link coordinates in colocObjectOverlap,
% colocObjectLinks, colocObjectClassImage (+ colocObjectClassColours) and
% colocObjectWriteBack.
%
% Layout (100x100 px, 0.1 um/px, one cell):
%   cisternae (A, channel 1): C1 ordinary r11:20 c11:20, C2 stream r11:20
%                             c41:50, C3 stream r61:70 c61:70
%   organelles (B, channel 2): O1 inside C1, O2 inside C2, O3 isolated
%                             r85:87 c85:87, O4 r22:24 c15:17 (2 px below
%                             C1 = 0.2 um -> in contact, no overlap)

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
            idx = sub2ind(testColocObjectVisualisation.Sz, R(:), C(:));
        end
        function [cis, org, O, P, S, N] = runLayout()
            r = @testColocObjectVisualisation.rect;
            cs = table({r(11:20,11:20); r(11:20,41:50); r(61:70,61:70)}, [1;2;3], [1;1;1], ...
                single([0.1; 0.8; 1.5]), 'VariableNames', {'cisternaePixelIdxList','cisternaeID','cellID','cisternaeSpeedMax'});
            cis = colocCisternaeObjects({cs}, 'all', 0.5);
            org = cell(2, 1);
            org{2} = table({r(15:17,15:17); r(15:17,45:47); r(85:87,85:87); r(22:24,15:17)}, (1:4)', ones(4,1), ...
                'VariableNames', {'organellePixelIdxList','organelleID','cellID'});
            [O, P, S, N] = colocObjectOverlap(cis, ones(testColocObjectVisualisation.Sz), 1, 2, 0.1, 0.2, 0, 'f', ...
                struct('statsB', {org}, 'kNearest', 1, 'partnerAttributes', {{'isStream','cisternaeSpeedMax'}}));
        end
        function c = px(rgb, row, col)
            c = double(squeeze(rgb(row, col, 1, 1, 1, :)))' / 255;
        end
    end

    methods (Test)
        function testPopulationAndLinkColumns(tc)
            [~, ~, O] = testColocObjectVisualisation.runLayout();
            O = O{1};
            tc.verifyEqual(sort(unique(O.population))', {'A','B'});
            OB = sortrows(O(strcmp(O.population,'B'), :), 'organelleID');
            % O3: own closest pixel (85,85) -> C3 corner (70,70), x = col
            tc.verifyEqual([OB.linkX1(3) OB.linkY1(3) OB.linkX2(3) OB.linkY2(3)], [85 85 70 70]);
            tc.verifyEqual(OB.nearestPartnerDistance(3), hypot(15,15)*0.1, 'AbsTol', 1e-9);
            % O4: top row 22 -> C1 bottom row 20, same column
            tc.verifyEqual(OB.linkY1(4) - OB.linkY2(4), 2);
            tc.verifyEqual(OB.linkX1(4), OB.linkX2(4));
        end

        function testLinksSkipOverlapping(tc)
            [~, ~, O] = testColocObjectVisualisation.runLayout();
            L = colocObjectLinks(O, 'B');
            L = L{1};
            tc.verifyTrue(all(isnan(L(1,:))));
            segs = reshape(L(2:end,:)', 3, 3, []);        % [x y d] x 3 rows x nSeg
            tc.verifyEqual(size(segs, 3), 2);             % O3, O4 (O1/O2 overlap -> skipped)
            tc.verifyEqual(sort(squeeze(segs(3,1,:)))', sort([0.2 hypot(15,15)*0.1]), 'AbsTol', 1e-9);
        end

        function testClassImageColours(tc)
            [cis, org, O] = testColocObjectVisualisation.runLayout();
            rgb = colocObjectClassImage(O, cis, org, 1, 2, testColocObjectVisualisation.Sz);
            cmap = colocObjectClassColours();
            tc.verifySize(rgb, [100 100 1 1 1 3]);
            tc.verifyClass(rgb, 'uint8');
            cmap = round(255 * cmap) / 255;         % uint8 quantisation (grey 0.5 -> 128/255)
            px = @(r,c) testColocObjectVisualisation.px(rgb, r, c);
            tc.verifyEqual(px(16,16), cmap(7,:));   % O1 inside C1: overlap pixel (white)
            tc.verifyEqual(px(86,86), cmap(3,:));   % O3 isolated: grey
            tc.verifyEqual(px(23,16), cmap(4,:));   % O4 in contact: yellow
            tc.verifyEqual(px(11,11), cmap(1,:));   % C1 edge: ordinary outline (cyan)
            tc.verifyEqual(px(11,41), cmap(2,:));   % C2 edge: stream outline (red)
            tc.verifyEqual(px(15,15), cmap(7,:));   % O2/O1 pixels overlap -> white
            tc.verifyEqual(px(50,5), [0 0 0]);      % background
        end

        function testClassImageStreamPrecedence(tc)
            % an organelle straddling an ordinary and a stream cisterna is magenta
            r = @testColocObjectVisualisation.rect;
            cs = table({r(11:20,11:20); r(11:20,23:32)}, [1;2], [1;1], single([0.1; 0.9]), ...
                'VariableNames', {'cisternaePixelIdxList','cisternaeID','cellID','cisternaeSpeedMax'});
            cis = colocCisternaeObjects({cs}, 'all', 0.5);
            org = cell(2,1);
            org{2} = table({r(24:26,18:27)}, 1, 1, 'VariableNames', {'organellePixelIdxList','organelleID','cellID'});
            org{2}.organellePixelIdxList{1} = [r(18:20,19:26); r(21:23,19:26)];   % overlaps both cisternae
            O = colocObjectOverlap(cis, ones(100), 1, 2, 0.1, 0.2, 0, 'f', ...
                struct('statsB', {org}, 'kNearest', 0, 'partnerAttributes', {{'isStream'}}));
            rgb = colocObjectClassImage(O, cis, org, 1, 2, [100 100]);
            cmap = colocObjectClassColours();
            tc.verifyEqual(testColocObjectVisualisation.px(rgb, 22, 21), cmap(6,:));  % B-only pixel: magenta
        end

        function testWriteBack(tc)
            [cis, org, O] = testColocObjectVisualisation.runLayout();
            target = org;
            target{1} = table((1:2)', [1;1], 'VariableNames', {'organelleID','cellID'});   % untouched slot
            cols = {'overlapAreaStream','nearestPartnerIsStream','inContact'};
            out = colocObjectWriteBack(target, O, 'B', 2, 'organelleID', 'erObj', cols);
            T = out{2};
            tc.verifyEqual(T.erObjNearestPartnerIsStream, [0; 1; 1; 0]);
            tc.verifyEqual(T.erObjInContact, [1; 1; 0; 1]);
            tc.verifyEqual(T.erObjOverlapAreaStream, [0; 9*0.01; 0; 0], 'AbsTol', 1e-12);
            tc.verifyTrue(all(isnan(out{1}.erObjInContact)), 'other slots get NaN columns');
            % idempotent: re-running replaces rather than duplicating
            out2 = colocObjectWriteBack(out, O, 'B', 2, 'organelleID', 'erObj', cols);
            tc.verifyEqual(out2{2}.Properties.VariableNames, T.Properties.VariableNames);
            tc.verifyEqual(out2{2}.erObjInContact, T.erObjInContact);
            % cisterna side (A) onto a cisternaeStats-style table
            cs = {cis{1}(:, {'organelleID','cellID'})};
            cs{1}.Properties.VariableNames{'organelleID'} = 'cisternaeID';
            outA = colocObjectWriteBack(cs, O, 'A', 1, 'cisternaeID', 'erObj', {'nPartners'});
            tc.verifyEqual(outA{1}.erObjNPartners, [1; 1; 0]);   % overlapping partners: C1-O1, C2-O2 (O4 only in contact)
            % both populations into ONE array in one call: neither wipes the other
            both = {cis{1}(:, {'organelleID','cellID'}); org{2}(:, {'organelleID','cellID'})};
            outAB = colocObjectWriteBack(both, O, {'A','B'}, [1 2], 'organelleID', 'orgObj', {'nPartners'});
            tc.verifyEqual(outAB{1}.orgObjNPartners, [1; 1; 0]);
            tc.verifyEqual(outAB{2}.orgObjNPartners, [1; 1; 0; 0]);
        end
    end
end
