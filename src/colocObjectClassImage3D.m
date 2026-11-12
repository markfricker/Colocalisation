function rgb = colocObjectClassImage3D(objectOut, labelsIn, chA, chB)
%COLOCOBJECTCLASSIMAGE3D  RGB class volume of a 3D object-colocalisation run
% (colocObjectOverlap3D), for checking the results section by section.
%
%   rgb = colocObjectClassImage3D(objectOut, labelsIn, chA, chB)
%
% The 3D counterpart of colocObjectClassImage, with the same colours and
% meaning (colocObjectClassColours): B objects filled by status (no
% contact / in contact / colocalised), A objects outlined in every section
% (4-neighbour outline within the section, so each slice shows them), and
% voxels where A and B coincide white. Only objects that took part in the
% run (appear in objectOut) are drawn.
%
% INPUTS
%   objectOut - {1 x 1 x nT} per-object tables from colocObjectOverlap3D
%               (population, organelleID, colocalised, inContact).
%   labelsIn  - [nY nX lC nZ nT] label volume the run used; channel c is
%               label channel min(lC, c), organelleID = label value.
%   chA, chB  - channels of populations A and B.
%
% OUTPUT
%   rgb - uint8 [nY x nX x 1 x nZ x nT x 3] (the app's RGB image layout).

cmap = uint8(round(255 * colocObjectClassColours()));
[nY, nX, lC, nZ, nT] = size(labelsIn);
rgb = zeros(nY, nX, 1, nZ, nT, 3, 'uint8');
for iT = 1:min(nT, size(objectOut, 3))
    O = objectOut{1, 1, iT};
    if ~istable(O) || isempty(O)
        continue
    end
    LA = double(reshape(labelsIn(:,:,min(lC,chA),:,iT), nY, nX, nZ));
    LB = double(reshape(labelsIn(:,:,min(lC,chB),:,iT), nY, nX, nZ));

    % --- B objects, filled by status -----------------------------------------
    OB = O(strcmp(O.population, 'B'), :);
    lutB = zeros(max(LB(:)) + 1, 1, 'uint8');       % label + 1 -> class
    cls = zeros(size(OB, 1), 1, 'uint8') + 3;
    cls(OB.inContact > 0)   = 4;
    cls(OB.colocalised > 0) = 5;
    ok = OB.organelleID >= 1 & OB.organelleID <= numel(lutB) - 1;
    lutB(OB.organelleID(ok) + 1) = cls(ok);
    C = lutB(LB + 1);
    Bmask = C > 0;

    % --- A objects, outlined in each section ----------------------------------
    OA = O(strcmp(O.population, 'A'), :);
    inRun = false(max(LA(:)) + 1, 1);
    okA = OA.organelleID >= 1 & OA.organelleID <= numel(inRun) - 1;
    inRun(OA.organelleID(okA) + 1) = true;
    LA(~inRun(LA + 1)) = 0;
    Amask = LA > 0;
    C(outline2D(LA)) = 1;

    % --- overlap voxels --------------------------------------------------------
    C(Amask & Bmask) = 7;

    for ch = 1:3
        lut = [0; cmap(:, ch)];
        rgb(:, :, 1, :, iT, ch) = reshape(lut(double(C) + 1), nY, nX, 1, nZ);
    end
end
end % colocObjectClassImage3D


function E = outline2D(L)
% labelled voxels with an in-section 4-neighbour of a different label
P = padarray(L, [1 1 0], 0);
C = P(2:end-1, 2:end-1, :);
E = C > 0 & (P(1:end-2,2:end-1,:) ~= C | P(3:end,2:end-1,:) ~= C | ...
             P(2:end-1,1:end-2,:) ~= C | P(2:end-1,3:end,:) ~= C);
end
