classdef testColocConformance < matlab.unittest.TestCase
%TESTCOLOCCONFORMANCE  Pixel colocalisation against known answers and Fiji Coloc 2.
%
% Conformance tests for the method-provenance records coloc.pearson,
% coloc.costes, coloc.manders (AnalyzERproject_sandbox core/utils/
% funcMethodDocs.m, field .provenance). Three kinds of test:
%
%   analytical  synthetic images with an exact expected answer
%               (colocConformanceCases: identical, anticorrelated,
%               partial overlap, Costes' signal-plus-background)
%   golden      the same images run through Fiji Coloc 2 3.1.0
%               (fixtures/coloc2_golden.mat, made by fixtures/makeColocGolden
%               + coloc2Golden.groovy; provenance stored in the .mat)
%   deviation   one test per documented deviation from Coloc 2
%               (testDev*), pinning our behaviour AND Coloc 2's, so a change
%               on either side fails here and the record must be revisited
%
% Run from this folder:
%   results = run(matlab.unittest.TestSuite.fromClass(?testColocConformance))
%
% Golden tolerances: Pearson's r and the regression line are the same
% closed forms (Coloc 2 accumulates in float -> AbsTol 1e-6 on r, RelTol
% 1e-3 on the line). Thresholds: Coloc 2 steps one grey level at a time and
% stores thresholds in the pixel type; we walk pixel by pixel on a
% continuous line -> within 1 grey level.

    properties
        cases
        golden
    end

    methods (TestClassSetup)
        function setup(tc)
            here = fileparts(mfilename('fullpath'));
            root = fullfile(here, '..', '..');
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'src')));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'utils')));
            tc.applyFixture(matlab.unittest.fixtures.PathFixture(here));
            tc.cases = colocConformanceCases();
            S = load(fullfile(here, 'fixtures', 'coloc2_golden.mat'), 'golden');
            tc.golden = S.golden;
        end
    end

    methods (Access = private)
        function c = case_(tc, name)
            c = tc.cases(strcmp({tc.cases.name}, name));
            c.ch1 = double(c.ch1);
            c.ch2 = double(c.ch2);
        end

        function v = g(tc, caseName, quantity, impl)
            % Coloc 2 value; impl defaults to 'Costes' (the dialog default).
            if nargin < 4
                impl = "Costes";
            end
            V = tc.golden.values;
            v = V.num(V.caseName == caseName & V.implementation == impl & V.name == quantity);
            tc.assertNumElements(v, 1, sprintf('golden %s / %s / %s', caseName, impl, quantity));
        end
    end

    % =====================================================================
    % Fixture provenance
    % =====================================================================
    methods (Test)
        function testGoldenProvenanceRecorded(tc)
            P = tc.golden.provenance;
            tc.verifyTrue(contains(P.versions, 'Colocalisation_Analysis-3.1.0'));
            tc.verifyNotEmpty(P.date);
            tc.verifyNotEmpty(P.settings);
            tc.verifyTrue(isfile(fullfile(fileparts(mfilename('fullpath')), 'fixtures', 'coloc2Golden.groovy')));
        end

        function testFixtureImagesMatchCaseGenerator(tc)
            % The TIFFs Coloc 2 read are the images the tests use.
            imDir = fullfile(fileparts(mfilename('fullpath')), 'fixtures', 'images');
            for k = 1:numel(tc.cases)
                c = tc.cases(k);
                tc.verifyEqual(imread(fullfile(imDir, [c.name '_ch1.tif'])), c.ch1, c.name);
                tc.verifyEqual(imread(fullfile(imDir, [c.name '_ch2.tif'])), c.ch2, c.name);
            end
        end
    end

    % =====================================================================
    % Analytical cases
    % =====================================================================
    methods (Test)
        function testIdenticalChannels(tc)
            c = tc.case_('identical');
            tc.verifyEqual(colocPearson(c.ch1, c.ch2, []), 1, 'AbsTol', 1e-12);
            [M1, M2] = colocManders(c.ch1, c.ch2, [], 0, 0);
            tc.verifyEqual([M1 M2], [1 1], 'AbsTol', 1e-12);
            [T1, T2, ~, conv, full, a, b] = colocCostesThreshold(c.ch1, c.ch2, []);
            tc.verifyEqual([a b], [1 0], 'AbsTol', 1e-9);
            tc.verifyTrue(conv);
            tc.verifyFalse(full);              % r of the dim pixels never reaches 0
            tc.verifyEqual(T1, T2, 'AbsTol', 1e-9);
        end

        function testAnticorrelatedChannels(tc)
            c = tc.case_('anticorrelated');
            tc.verifyEqual(colocPearson(c.ch1, c.ch2, []), -1, 'AbsTol', 1e-12);
            [T1, ~, ~, conv, ~, a, b] = colocCostesThreshold(c.ch1, c.ch2, []);
            tc.verifyFalse(conv);
            tc.verifyTrue(isnan(T1));
            % Major axis of ch2 = 255 - ch1 (slope fix: the minor-axis root
            % gave +1 before 2026-10-19).
            tc.verifyEqual([a b], [-1 255], 'AbsTol', 1e-9);
        end

        function testIndependentNoise(tc)
            rng(7);
            ch1 = randi(255, 200);
            ch2 = randi(255, 200);
            r = colocPearson(ch1, ch2, []);
            tc.verifyLessThan(abs(r), 4 / sqrt(numel(ch1)));   % 4 SD of r under H0
            p = colocCostesRandomization(ch1, ch2, [], 100, 3, 8);
            tc.verifyGreaterThan(p, 0.05);
        end

        function testPartialOverlapExactManders(tc)
            % Squares A, B of equal area overlapping by half, zero
            % background: M1 = M2 = |A n B| / |A| = 0.5 at any thresholds
            % between 0 and the signal, and r of the two indicator images
            % = (pAB - pA pB) / (pA (1 - pA)) = (1/8 - 1/16) / (3/16) = 1/3.
            c = tc.case_('partialOverlap');
            tc.verifyEqual(colocPearson(c.ch1, c.ch2, []), 1/3, 'AbsTol', 1e-12);
            [M1, M2] = colocManders(c.ch1, c.ch2, [], 0, 0);
            tc.verifyEqual([M1 M2], [0.5 0.5], 'AbsTol', 1e-12);
            [T1, T2, ~, ~, full] = colocCostesThreshold(c.ch1, c.ch2, []);
            tc.verifyTrue(full);
            [tM1, tM2] = colocManders(c.ch1, c.ch2, [], T1, T2);
            tc.verifyEqual([tM1 tM2], [0.5 0.5], 'AbsTol', 1e-12);
        end

        function testCostesClassicThresholdsBelowSignalBand(tc)
            % Correlated signal 80-200 on independent background 0-40:
            % the Costes threshold falls below the signal band and the
            % pixels below it have r <= 0.
            c = tc.case_('costesClassic');
            [T1, T2, rAt, ~, full] = colocCostesThreshold(c.ch1, c.ch2, []);
            tc.verifyTrue(full);
            tc.verifyLessThanOrEqual(rAt, 0);
            tc.verifyLessThan(T1, 80);
            tc.verifyLessThan(T2, 0.8 * 80 + 10);
            tc.verifyLessThanOrEqual(colocPearson(c.ch1, c.ch2, c.ch1 < T1 | c.ch2 < T2), 0);
            [tM1, tM2] = colocManders(c.ch1, c.ch2, [], T1, T2);
            blk = c.sigA;
            tc.verifyGreaterThanOrEqual(tM1, sum(c.ch1(blk)) / sum(c.ch1(:)));
            tc.verifyGreaterThanOrEqual(tM2, 0.99 * sum(c.ch2(blk)) / sum(c.ch2(:)));
        end
    end

    % =====================================================================
    % Golden values (Fiji Coloc 2 3.1.0)
    % =====================================================================
    methods (Test)
        function testGoldenPearson(tc)
            for k = 1:numel(tc.cases)
                c = tc.case_(tc.cases(k).name);
                tc.verifyEqual(colocPearson(c.ch1, c.ch2, []), ...
                    tc.g(c.name, "Pearson's R value (no threshold)"), 'AbsTol', 1e-6, c.name);
            end
        end

        function testGoldenRegressionLine(tc)
            for k = 1:numel(tc.cases)
                c = tc.case_(tc.cases(k).name);
                [~, ~, ~, ~, ~, a, b] = colocCostesThreshold(c.ch1, c.ch2, []);
                tc.verifyEqual(a, tc.g(c.name, "m (slope)"), 'RelTol', 1e-3, 'AbsTol', 1e-6, c.name);
                tc.verifyEqual(b, tc.g(c.name, "b (y-intercept)"), 'RelTol', 1e-3, 'AbsTol', 1e-3, c.name);
            end
        end

        function testGoldenUnthresholdedManders(tc)
            for k = 1:numel(tc.cases)
                c = tc.case_(tc.cases(k).name);
                [M1, M2] = colocManders(c.ch1, c.ch2, [], 0, 0);
                tc.verifyEqual(M1, tc.g(c.name, "Manders' M1 (Above zero intensity of Ch2)"), 'AbsTol', 1e-9, c.name);
                tc.verifyEqual(M2, tc.g(c.name, "Manders' M2 (Above zero intensity of Ch1)"), 'AbsTol', 1e-9, c.name);
            end
        end

        function testGoldenMandersAndBelowRAtColoc2Thresholds(tc)
            % Given Coloc 2's own thresholds, our tM1/tM2 and below-threshold
            % r reproduce Coloc 2 exactly: the formulas are the same.
            for name = ["independent" "anticorrelated" "partialOverlap" "partialNoisy" "costesClassic"]
                c = tc.case_(char(name));
                T1 = tc.g(name, "Ch1 Max Threshold");
                T2 = tc.g(name, "Ch2 Max Threshold");
                [tM1, tM2] = colocManders(c.ch1, c.ch2, [], T1, T2);
                tc.verifyEqual(tM1, tc.g(name, "Manders' tM1 (Above autothreshold of Ch2)"), 'AbsTol', 1e-9, char(name));
                tc.verifyEqual(tM2, tc.g(name, "Manders' tM2 (Above autothreshold of Ch1)"), 'AbsTol', 1e-9, char(name));
                tc.verifyEqual(colocPearson(c.ch1, c.ch2, c.ch1 < T1 | c.ch2 < T2), ...
                    tc.g(name, "Pearson's R value (below threshold)"), 'AbsTol', 1e-6, char(name));
            end
        end

        function testGoldenCostesThresholdsWithinOneGreyLevel(tc)
            % Cases where Coloc 2's Costes stepper reaches r <= 0 normally.
            for name = ["independent" "partialOverlap" "partialNoisy" "costesClassic"]
                c = tc.case_(char(name));
                [T1, T2, ~, ~, full] = colocCostesThreshold(c.ch1, c.ch2, []);
                tc.verifyTrue(full, char(name));
                tc.verifyLessThanOrEqual(abs(T1 - tc.g(name, "Ch1 Max Threshold")), 1, char(name));
                tc.verifyLessThanOrEqual(abs(T2 - tc.g(name, "Ch2 Max Threshold")), 1, char(name));
            end
        end
    end

    % =====================================================================
    % Deviations from Coloc 2 (one test per funcMethodDocs deviation)
    % =====================================================================
    methods (Test)
        function testDevMinRFallback(tc)
            % coloc.costes: no zero crossing -> the global-minimum-r
            % threshold. Coloc 2's Costes stepper instead runs to its floor
            % (threshold 1, below-r NaN, 'thresholds too low'); its
            % Bisection stepper lands where we do.
            c = tc.case_('identical');
            [T1, T2, ~, conv, full] = colocCostesThreshold(c.ch1, c.ch2, []);
            tc.verifyTrue(conv);
            tc.verifyFalse(full);
            tc.verifyEqual(tc.g("identical", "Ch1 Max Threshold"), 1);
            tc.verifyTrue(isnan(tc.g("identical", "Pearson's R value (below threshold)")));
            tc.verifyEqual(T1, tc.g("identical", "Ch1 Max Threshold", "Bisection"), 'AbsTol', 1);
            tc.verifyEqual(T2, tc.g("identical", "Ch2 Max Threshold", "Bisection"), 'AbsTol', 1);
        end

        function testDevDerivedThresholdNotRounded(tc)
            % coloc.costes: the threshold of the channel read off the
            % regression line is kept continuous; Coloc 2 stores it in the
            % pixel type (rounded for 8/16-bit images). On integer data the
            % pixels at the rounded level then count for Coloc 2 and not
            % for us, which moves tM1 on this case by about 0.05. Kept on
            % purpose (agreed 2026-10-19): the continuous value is the
            % exact point on the regression line.
            c = tc.case_('independent');
            [T1, T2] = colocCostesThreshold(c.ch1, c.ch2, []);
            T2c = tc.g("independent", "Ch2 Max Threshold");
            tc.verifyNotEqual(T2, round(T2));                 % continuous
            tc.verifyEqual(round(T2), T2c);                   % Coloc 2 = rounded
            tM1 = colocManders(c.ch1, c.ch2, [], T1, T2);
            tM1rounded = colocManders(c.ch1, c.ch2, [], T1, round(T2));
            tM1c = tc.g("independent", "Manders' tM1 (Above autothreshold of Ch2)");
            tc.verifyEqual(tM1rounded, tM1c, 'AbsTol', 1e-9);
            tc.verifyLessThan(tM1, tM1c - 0.02);
        end

        function testDevNegativeSlopeGivesNoThreshold(tc)
            % coloc.costes: a negative regression slope gives no threshold
            % (Costes' method assumes the channels rise together); Coloc 2
            % walks the line anyway and reports thresholds and tM1/tM2.
            c = tc.case_('anticorrelated');
            [T1, T2, ~, conv] = colocCostesThreshold(c.ch1, c.ch2, []);
            tc.verifyFalse(conv);
            tc.verifyTrue(isnan(T1) && isnan(T2));
            tc.verifyFalse(isnan(tc.g("anticorrelated", "Ch1 Max Threshold")));
            tc.verifyFalse(isnan(tc.g("anticorrelated", "Manders' tM1 (Above autothreshold of Ch2)")));
        end

        function testDevPValueIsEmpiricalFraction(tc)
            % coloc.costes (significance): p = fraction of block-shuffled r
            % >= observed r, small = significant. Coloc 2's "Costes
            % P-Value" is Phi((r - mean_null) / sd_null), the normal CDF,
            % so ~1 = significant: the two are not on the same scale.
            c = tc.case_('costesClassic');
            [p, robs, rNull] = colocCostesRandomization(c.ch1, c.ch2, [], 100, 3, 1);
            tc.verifyEqual(p, mean(rNull >= robs));
            tc.verifyEqual(p, 0);
            pc = tc.g("costesClassic", "Costes P-Value");
            tc.verifyEqual(pc, 1, 'AbsTol', 1e-6);
            mu = tc.g("costesClassic", "Costes Shuffled Mean");
            sd = tc.g("costesClassic", "Costes Shuffled Std.D.");
            tc.verifyEqual(pc, min(1, 0.5 * erfc(-(robs - mu) / (sd * sqrt(2)))), 'AbsTol', 1e-6);
        end

        function testDevMandersManualFallback(tc)
            % coloc.manders: with no Costes threshold (negative slope) M1/M2
            % use p.manualThreshold1/2 (AnalyzER passes background SD x
            % backgroundSdMultiplier); Coloc 2 reports tM1/tM2 at its
            % stepper's end point instead.
            c = tc.case_('anticorrelated');
            p = colocParamsDefault();
            p.nIterations = 10;
            p.manualThreshold1 = 50;
            p.manualThreshold2 = 120;
            R = colocPixelBasedRun(c.ch1, c.ch2, [], p);
            tc.verifyEqual([R.mandersThreshold1 R.mandersThreshold2], [50 120]);
            [M1, M2] = colocManders(c.ch1, c.ch2, [], 50, 120);
            tc.verifyEqual([R.manders1 R.manders2], [M1 M2], 'AbsTol', 1e-12);
            tc.verifyNotEqual(R.manders1, tc.g("anticorrelated", "Manders' tM1 (Above autothreshold of Ch2)"));
        end

        function testGoldenNullWidth(tc)
            % coloc.costes (significance): the shuffled channel is
            % Gaussian-smoothed (sigma = blockSize) as in Coloc 2, so the
            % null has Coloc 2's width. Both nulls come from 100 random
            % shuffles (SD of an SD estimate ~7%), hence the 30% band.
            for k = 1:numel(tc.cases)
                c = tc.case_(tc.cases(k).name);
                [~, ~, rNull] = colocCostesRandomization(c.ch1, c.ch2, [], 100, 3, 1);
                sdC = tc.g(c.name, "Costes Shuffled Std.D.");
                tc.verifyEqual(std(rNull), sdC, 'RelTol', 0.3, c.name);
                tc.verifyLessThan(abs(mean(rNull) - tc.g(c.name, "Costes Shuffled Mean")), 3 * sdC / sqrt(10), c.name);
            end
        end

        function testNullCalibratedOnIndependentTexture(tc)
            % The property that decided the smoothing: on independent
            % images with PSF-scale texture (FWHM 5.9 px) p < 0.05 should
            % occur for ~5% of pairs. Unsmoothed block shuffling gave
            % ~25% (2026-10-19); smoothed gives ~4%.
            rng(11);
            nPairs = 120;
            nFP = 0;
            for i = 1:nPairs
                a = imgaussfilt(randn(96), 2.5);
                b = imgaussfilt(randn(96), 2.5);
                p = colocCostesRandomization(a, b, [], 50, 3, i);
                nFP = nFP + (p < 0.05);
            end
            tc.verifyLessThanOrEqual(nFP / nPairs, 0.10);
        end
    end
end
