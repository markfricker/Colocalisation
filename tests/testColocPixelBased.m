classdef testColocPixelBased < matlab.unittest.TestCase
%TESTCOLOCPIXELBASED  Unit tests for the pixel-based colocalisation core.
%
% USAGE
%   Run all tests from the MATLAB command window:
%       results = runtests('tests/testColocPixelBased');
%       table(results)
%
% COVERAGE
%   colocPearson              — 6 tests
%   colocManders               — 4 tests
%   colocCostesThreshold       — 3 tests
%   colocCostesRandomization   — 4 tests
%   colocScoreImage            — 6 tests
%   colocPixelBasedRun         — 5 tests
%
% REQUIREMENTS
%   MATLAB R2019b+ (matlab.unittest framework)
%   Statistics and Machine Learning Toolbox not required (no toolbox
%   functions used beyond base MATLAB).
%   All src/ and utils/ on the path (added automatically by
%   TestClassSetup).
%
% All synthetic test data is generated with a fixed rng seed per test for
% reproducibility.

    methods (TestClassSetup)
        function addPaths(tc) %#ok<MANU>
            rootDir = fullfile(fileparts(mfilename('fullpath')), '..');
            addpath(fullfile(rootDir, 'src'));
            addpath(fullfile(rootDir, 'utils'));
        end
    end

    % =====================================================================
    % Shared synthetic helpers
    % =====================================================================
    methods (Static, Access = private)
        function [ch1, ch2] = correlatedPair(imSize, slope, intercept, noiseStd, seedVal)
            rng(seedVal);
            ch1 = 100 * rand(imSize);
            ch2 = slope .* ch1 + intercept + noiseStd * randn(imSize);
        end

        function [ch1, ch2] = uncorrelatedPair(imSize, seedVal)
            rng(seedVal);
            ch1 = 100 * rand(imSize);
            ch2 = 100 * rand(imSize);
        end

        function [ch1, ch2, sigMask] = signalPlusBackgroundPair(imSize, seedVal)
            % Dim, mutually independent background everywhere (0-20), with
            % a bright, near-perfectly-correlated "signal" block (50-100)
            % overwritten in the centre -- used to check that
            % colocCostesThreshold finds a threshold separating the two.
            rng(seedVal);
            ch1 = 20 * rand(imSize);
            ch2 = 20 * rand(imSize);

            sigMask = false(imSize);
            rows = round(imSize(1)*0.35):round(imSize(1)*0.65);
            cols = round(imSize(2)*0.35):round(imSize(2)*0.65);
            sigMask(rows, cols) = true;

            nSig = nnz(sigMask);
            sigVals1 = 50 + 50 * rand(nSig, 1);
            ch1(sigMask) = sigVals1;
            ch2(sigMask) = sigVals1 + 2 * randn(nSig, 1);
        end

        function [ch1, ch2] = wholeRangeCorrelatedPair(imSize, seedVal)
            % Exact (noise-free) positive linear relationship across the
            % WHOLE intensity range -- no independent background
            % population anywhere. This is the clean idealisation of a
            % tight cell-boundary mask (true background already excluded)
            % or two near-identical adjacent frames: because the
            % relationship holds exactly at every intensity, every
            % below-threshold subpopulation has r=1, so r never crosses
            % zero and colocCostesThreshold's min-r partial-convergence
            % fallback must engage. (Any amount of real per-pixel noise
            % lets a handful of the very dimmest pixels decorrelate by
            % chance, which is enough for the strict criterion to be met
            % right at the noise floor -- exact data is what reliably
            % forces the fallback path for a deterministic unit test.)
            rng(seedVal); %#ok<NASGU> -- kept for interface symmetry/reproducibility
            ch1 = 100 * rand(imSize);
            ch2 = 0.8 * ch1 + 5;
        end
    end

    % =====================================================================
    % colocPearson
    % =====================================================================
    methods (Test)
        function testPearsonPerfectPositiveCorrelation(tc)
            [ch1, ch2] = tc.correlatedPair([40 40], 2, 5, 0, 1);
            r = colocPearson(ch1, ch2, []);
            tc.verifyEqual(r, 1, 'AbsTol', 1e-9);
        end

        function testPearsonPerfectNegativeCorrelation(tc)
            [ch1, ch2] = tc.correlatedPair([40 40], -3, 10, 0, 2);
            r = colocPearson(ch1, ch2, []);
            tc.verifyEqual(r, -1, 'AbsTol', 1e-9);
        end

        function testPearsonNearZeroForIndependentNoise(tc)
            [ch1, ch2] = tc.uncorrelatedPair([120 120], 3);
            r = colocPearson(ch1, ch2, []);
            tc.verifyLessThan(abs(r), 0.05);
        end

        function testPearsonRespectsMask(tc)
            % Whole image is uncorrelated noise except a sub-block that is
            % perfectly correlated; masking to that sub-block should give
            % r close to 1 even though the unmasked image does not.
            [ch1, ch2] = tc.uncorrelatedPair([60 60], 4);
            mask = false(60, 60);
            mask(10:30, 10:30) = true;
            ch2(mask) = ch1(mask);   % force perfect correlation inside mask

            rMasked   = colocPearson(ch1, ch2, mask);
            rUnmasked = colocPearson(ch1, ch2, []);

            tc.verifyEqual(rMasked, 1, 'AbsTol', 1e-9);
            tc.verifyLessThan(abs(rUnmasked), 0.9);
        end

        function testPearsonNaNWhenTooFewPixels(tc)
            ch1 = [1 2; 3 4];
            ch2 = [4 3; 2 1];
            mask = false(2, 2);
            mask(1, 1) = true;   % only 1 pixel selected
            [r, n] = colocPearson(ch1, ch2, mask);
            tc.verifyTrue(isnan(r));
            tc.verifyEqual(n, 1);
        end

        function testPearsonNaNForZeroVarianceChannel(tc)
            ch1 = 5 * ones(10, 10);      % constant -> zero variance
            ch2 = rand(10, 10);
            r = colocPearson(ch1, ch2, []);
            tc.verifyTrue(isnan(r));
        end
    end

    % =====================================================================
    % colocManders
    % =====================================================================
    methods (Test)
        function testMandersFullOverlapGivesOne(tc)
            ch = rand(20, 20) + 0.1;   % strictly positive
            [M1, M2] = colocManders(ch, ch, [], 0, 0);
            tc.verifyEqual(M1, 1, 'AbsTol', 1e-9);
            tc.verifyEqual(M2, 1, 'AbsTol', 1e-9);
        end

        function testMandersDisjointRegionsGiveZero(tc)
            ch1 = zeros(20, 20);
            ch2 = zeros(20, 20);
            ch1(1:10, :)  = rand(10, 20) + 0.1;   % ch1 signal, top half only
            ch2(11:20, :) = rand(10, 20) + 0.1;   % ch2 signal, bottom half only
            [M1, M2] = colocManders(ch1, ch2, [], 0, 0);
            tc.verifyEqual(M1, 0);
            tc.verifyEqual(M2, 0);
        end

        function testMandersPartialOverlapMatchesHandComputed(tc)
            % ch1 signal in cols 1:10, ch2 signal in cols 6:15 (5-col
            % overlap out of 10 ch1-signal cols and 10 ch2-signal cols) --
            % with uniform intensity, M1 and M2 both reduce to the simple
            % fractional-column overlap.
            ch1 = zeros(20, 20);
            ch2 = zeros(20, 20);
            ch1(:, 1:10)  = 1;
            ch2(:, 6:15)  = 1;
            [M1, M2] = colocManders(ch1, ch2, [], 0, 0);
            tc.verifyEqual(M1, 5/10, 'AbsTol', 1e-9);
            tc.verifyEqual(M2, 5/10, 'AbsTol', 1e-9);
        end

        function testMandersNaNWhenNoAboveThresholdSignal(tc)
            ch1 = zeros(10, 10);
            ch2 = rand(10, 10) + 0.1;
            [M1, M2] = colocManders(ch1, ch2, [], 0, 0);
            tc.verifyTrue(isnan(M1));
            tc.verifyFalse(isnan(M2));
        end
    end

    % =====================================================================
    % colocCostesThreshold
    % =====================================================================
    methods (Test)
        function testCostesThresholdSeparatesSignalFromBackground(tc)
            [ch1, ch2] = tc.signalPlusBackgroundPair([100 100], 5);
            [Tr, Tg, rAtThreshold, converged, fullyConverged] = colocCostesThreshold(ch1, ch2, []);

            tc.verifyTrue(converged);
            tc.verifyTrue(fullyConverged);
            tc.verifyLessThan(Tr, 30);          % well below the 50-100 signal band
            tc.verifyLessThanOrEqual(rAtThreshold, 1e-9);

            % Above-threshold population should recover strong overlap
            % (the near-perfectly-correlated signal block).
            [M1, M2] = colocManders(ch1, ch2, [], Tr, Tg);
            tc.verifyGreaterThan(M1, 0.7);
            tc.verifyGreaterThan(M2, 0.7);
        end

        function testCostesThresholdPartialConvergenceUsesMinRFallback(tc)
            % When r never crosses <=0 (no genuinely uncorrelated
            % background population -- e.g. a tight mask, or two very
            % similar-looking frames), the strict criterion isn't met but
            % a usable threshold is still returned via the min-r fallback.
            % No console warning -- this is an expected, handled outcome
            % (converged/fullyConverged already tell the caller), not an
            % exceptional one; firing per-cell now that colocalisation
            % runs per real cell region would flood the console.
            [ch1, ch2] = tc.wholeRangeCorrelatedPair([120 120], 51);

            [Tr, Tg, rAtThreshold, converged, fullyConverged] = colocCostesThreshold(ch1, ch2, []);

            tc.verifyTrue(converged);
            tc.verifyFalse(fullyConverged);
            tc.verifyFalse(isnan(Tr));
            tc.verifyFalse(isnan(Tg));
            tc.verifyGreaterThan(rAtThreshold, 0);   % never actually reached <=0
        end

        function testCostesThresholdNonPositiveSlopeReturnsNaNNoWarning(tc)
            % No console warning -- expected, handled outcome (converged
            % tells the caller), not exceptional.
            [ch1, ch2] = tc.correlatedPair([50 50], -1, 100, 2, 6);

            [Tr, ~, ~, converged] = colocCostesThreshold(ch1, ch2, []);

            tc.verifyFalse(converged);
            tc.verifyTrue(isnan(Tr));
        end

        function testCostesThresholdTooFewPixelsReturnsNotConverged(tc)
            ch1 = rand(2, 2);
            ch2 = rand(2, 2);
            mask = false(2, 2);
            mask(1, 1) = true;   % only 1 pixel -> below the n<3 floor
            [Tr, ~, ~, converged] = colocCostesThreshold(ch1, ch2, mask);
            tc.verifyFalse(converged);
            tc.verifyTrue(isnan(Tr));
        end
    end

    % =====================================================================
    % colocCostesRandomization
    % =====================================================================
    methods (Test)
        function testRandomizationSignificantForRealSignal(tc)
            [ch1, ch2] = tc.correlatedPair([60 60], 2, 5, 3, 10);
            [pValue, robs, rNull] = colocCostesRandomization(ch1, ch2, [], 200, 4, 42);

            tc.verifyLessThan(pValue, 0.05);
            tc.verifyGreaterThan(robs, 0.5);
            tc.verifyEqual(numel(rNull), 200);
        end

        function testRandomizationNotSignificantForIndependentNoise(tc)
            [ch1, ch2] = tc.uncorrelatedPair([80 80], 11);
            [pValue, robs] = colocCostesRandomization(ch1, ch2, [], 200, 4, 43);

            tc.verifyLessThan(abs(robs), 0.2);
            tc.verifyGreaterThan(pValue, 0.05);
        end

        function testRandomizationReproducibleWithSeed(tc)
            [ch1, ch2] = tc.correlatedPair([40 40], 1.5, 0, 5, 20);
            [p1, r1, null1] = colocCostesRandomization(ch1, ch2, [], 50, 4, 7);
            [p2, r2, null2] = colocCostesRandomization(ch1, ch2, [], 50, 4, 7);

            tc.verifyEqual(p1, p2);
            tc.verifyEqual(r1, r2);
            tc.verifyEqual(null1, null2);
        end

        function testRandomizationBlockTooLargeErrors(tc)
            ch1 = rand(10, 10);
            ch2 = rand(10, 10);
            tc.verifyError(@() colocCostesRandomization(ch1, ch2, [], 10, 20), ...
                'colocCostesRandomization:blockTooLarge');
        end
    end

    % =====================================================================
    % colocScoreImage
    % =====================================================================
    methods (Test)
        function testScoreImageGatesBelowThresholdToNaN(tc)
            ch1 = [2 8; 12 20];
            ch2 = [3 9; 11 21];
            [scoreImage, gatedMask] = colocScoreImage(ch1, ch2, [], 10, 10, 'geomean');

            tc.verifyEqual(gatedMask, logical([0 0; 1 1]));
            tc.verifyTrue(isnan(scoreImage(1,1)));
            tc.verifyTrue(isnan(scoreImage(1,2)));
            tc.verifyFalse(isnan(scoreImage(2,1)));
            tc.verifyFalse(isnan(scoreImage(2,2)));
        end

        function testScoreImageGeomeanRangeAndBrightestScoresOne(tc)
            % A pixel that is simultaneously the brightest in both
            % channels (among gated pixels) should score exactly 1; every
            % gated score should fall in [0,1].
            ch1 = [15 25 100];
            ch2 = [16 30 90];
            scoreImage = colocScoreImage(ch1, ch2, true(1,3), 10, 10, 'geomean');

            tc.verifyEqual(scoreImage(3), 1, 'AbsTol', 1e-9);
            tc.verifyGreaterThanOrEqual(min(scoreImage), 0);
            tc.verifyLessThanOrEqual(max(scoreImage), 1);
        end

        function testScoreImagePdmSignConvention(tc)
            % Gated population is idx 3 and 6; idx3 sits below the whole-
            % mask mean in ch1 but above it in ch2 (negative PDM), idx6
            % sits above the mean in both (positive PDM).
            ch1 = [5, 8, 12, 30, 15, 20];
            ch2 = [5, 9, 25, 8, 5, 30];
            [scoreImage, gatedMask] = colocScoreImage(ch1, ch2, true(1,6), 10, 10, 'pdm');

            tc.verifyEqual(find(gatedMask), [3 6]);
            tc.verifyEqual(scoreImage(3), -34, 'AbsTol', 1e-3);
            tc.verifyEqual(scoreImage(6), 81.6667, 'AbsTol', 1e-3);
            tc.verifyTrue(isnan(scoreImage(1)));
        end

        function testScoreImageRespectsMask(tc)
            ch1 = [20 20];
            ch2 = [20 20];
            mask = logical([1 0]);
            [scoreImage, gatedMask] = colocScoreImage(ch1, ch2, mask, 5, 5, 'geomean');

            tc.verifyFalse(gatedMask(2));
            tc.verifyTrue(isnan(scoreImage(2)));
        end

        function testScoreImageAllBelowThresholdGivesAllNaN(tc)
            ch1 = [1 2; 3 4];
            ch2 = [1 2; 3 4];
            [scoreImage, gatedMask] = colocScoreImage(ch1, ch2, [], 100, 100, 'geomean');

            tc.verifyFalse(any(gatedMask, 'all'));
            tc.verifyTrue(all(isnan(scoreImage), 'all'));
        end

        function testScoreImageUnknownTypeErrors(tc)
            tc.verifyError(@() colocScoreImage([1 2], [1 2], [], 0, 0, 'bogus'), ...
                'colocScoreImage:unknownScoreType');
        end
    end

    % =====================================================================
    % colocPixelBasedRun (integration)
    % =====================================================================
    methods (Test)
        function testPixelBasedRunFieldsAndConsistency(tc)
            [ch1, ch2] = tc.correlatedPair([50 50], 2, 3, 2, 30);
            p = colocParamsDefault();
            p.nIterations = 30;   % keep the test fast
            results = colocPixelBasedRun(ch1, ch2, [], p);

            tc.verifyEqual(results.pearsonN, numel(ch1));
            tc.verifyGreaterThan(results.pearsonR, 0.5);
            tc.verifyTrue(isfield(results, 'manders1'));
            tc.verifyTrue(isfield(results, 'manders2'));
            tc.verifyTrue(isfield(results, 'costesFullyConverged'));
            tc.verifyEqual(size(results.pixelValues, 1), numel(ch1));
            tc.verifyEqual(numel(results.randNullR), p.nIterations);
            tc.verifyEqual(size(results.scoreImageGeomean), size(ch1));
            tc.verifyEqual(size(results.scoreImagePdm), size(ch1));
        end

        function testPixelBasedRunPixelValuesMatchMask(tc)
            mask = false(30, 30);
            mask(5:25, 5:25) = true;
            [ch1, ch2] = tc.correlatedPair([30 30], 1, 0, 1, 31);
            p = colocParamsDefault();
            p.nIterations = 20;
            p.blockSize   = 2;
            results = colocPixelBasedRun(ch1, ch2, mask, p);

            tc.verifyEqual(size(results.pixelValues, 1), nnz(mask));
            tc.verifyEqual(results.pixelValues(:,1), double(ch1(results.maskLinearIdx)));
            tc.verifyEqual(results.pixelValues(:,2), double(ch2(results.maskLinearIdx)));
        end

        function testPixelBasedRunDefaultSaturationIsInfNoExclusion(tc)
            p = colocParamsDefault();
            tc.verifyEqual(p.saturationValue1, Inf);
            tc.verifyEqual(p.saturationValue2, Inf);
        end

        function testPixelBasedRunExcludesSaturatedPixels(tc)
            [ch1, ch2] = tc.correlatedPair([40 40], 1, 0, 1, 40);
            ch1(1:5, 1:5) = 500;   % clipped/saturated block
            p = colocParamsDefault();
            p.nIterations = 20;
            p.saturationValue1 = 500;

            results = colocPixelBasedRun(ch1, ch2, [], p);

            tc.verifyEqual(results.pearsonN, numel(ch1) - 25);
            tc.verifyTrue(all(results.pixelValues(:,1) < 500));
        end

        function testPixelBasedRunBothScoreImagesAlwaysComputed(tc)
            [ch1, ch2] = tc.correlatedPair([40 40], 2, 5, 1, 41);
            p = colocParamsDefault();
            p.nIterations = 20;

            results = colocPixelBasedRun(ch1, ch2, [], p);

            tc.verifyTrue(any(~isnan(results.scoreImageGeomean), 'all'));
            tc.verifyTrue(any(~isnan(results.scoreImagePdm), 'all'));
            % Same gate for both -- identical NaN pattern.
            tc.verifyEqual(isnan(results.scoreImageGeomean), isnan(results.scoreImagePdm));
        end
    end
end
