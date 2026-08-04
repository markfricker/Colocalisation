function demoColocPixelBased(seedVal)
% DEMOCOLOCPIXELBASED
% Demo for the pixel-based colocalisation core (colocPixelBasedRun).
%
% Inputs
%   seedVal : optional rng seed for the synthetic image pair (default 1)
%
% This demo:
%   1) builds a synthetic two-channel image: dim, mutually independent
%      background plus a bright, near-perfectly-correlated "signal" block
%   2) runs the full pixel-based pipeline (Pearson, Costes threshold,
%      Manders M1/M2, Costes randomisation significance)
%   3) prints a summary and plots the two channels alongside a ch1/ch2
%      intensity scatter (the same pairing the eventual GUI plot axes
%      will show)

    if nargin < 1 || isempty(seedVal)
        seedVal = 1;
    end

    imSize = [150 150];
    rng(seedVal);
    ch1 = 20 * rand(imSize);
    ch2 = 20 * rand(imSize);

    sigMask = false(imSize);
    rows = 50:100;
    cols = 50:100;
    sigMask(rows, cols) = true;
    nSig = nnz(sigMask);
    sigVals1 = 50 + 50 * rand(nSig, 1);
    ch1(sigMask) = sigVals1;
    ch2(sigMask) = sigVals1 + 3 * randn(nSig, 1);

    p = colocParamsDefault();
    results = colocPixelBasedRun(ch1, ch2, [], p);

    fprintf('Pearson''s r (whole image):    %.3f (n=%d)\n', results.pearsonR, results.pearsonN);
    fprintf('Costes threshold converged:   %d\n', results.costesConverged);
    fprintf('  T1 = %.2f, T2 = %.2f\n', results.costesT1, results.costesT2);
    fprintf('Manders M1 / M2 (at threshold): %.3f / %.3f\n', results.manders1, results.manders2);
    fprintf('Costes randomisation p-value: %.3f\n', results.randPValue);

    figure('Color', 'w');

    subplot(1, 3, 1);
    imshow(ch1, []);
    title('Channel 1');

    subplot(1, 3, 2);
    imshow(ch2, []);
    title('Channel 2');

    subplot(1, 3, 3);
    scatter(results.pixelValues(:,1), results.pixelValues(:,2), 4, 'filled', ...
        'MarkerFaceAlpha', 0.3);
    hold on;
    if results.costesConverged
        xline(results.costesT1, 'r--');
        yline(results.costesT2, 'r--');
    end
    hold off;
    xlabel('Channel 1 intensity');
    ylabel('Channel 2 intensity');
    title('Pixel scatter (ch1 vs ch2)');
    axis square;
end
