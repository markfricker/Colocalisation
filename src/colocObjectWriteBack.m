function targetStats = colocObjectWriteBack(targetStats, objectOut, population, channel, idColumn, prefix, columns)
%COLOCOBJECTWRITEBACK  Copy per-object colocalisation results onto the
% object tables the app already displays/exports (organelle morphology,
% ER cisternae), as prefixed columns.
%
%   targetStats = colocObjectWriteBack(targetStats, objectOut, population, ...
%       channel, idColumn, prefix, columns)
%
% Once written back, the results can be shown with the app's existing
% colour-by-metric object display and appear on the morphology / cisternae
% Save Data sheets, with no new drawing code.
%
% INPUTS
%   targetStats - {nC x nZ x nT} object tables to update (e.g.
%                 morphologyStats or cisternaeStats).
%   objectOut   - {1 x nZ x nT} per-object tables from colocObjectOverlap.
%   population  - 'A' or 'B' (or a cellstr, e.g. {'A','B'}): which objectOut
%                 rows to copy. Pass BOTH populations in one call when they
%                 share the same targetStats and prefix -- separate calls
%                 would reset each other's columns.
%   channel     - channel slot(s) of targetStats those objects live in, one
%                 per population.
%   idColumn    - object ID column of targetStats matched against
%                 objectOut.organelleID (e.g. 'organelleID', 'cisternaeID').
%   prefix      - column prefix, e.g. 'erObj' -> erObjOverlapArea.
%   columns     - cellstr of objectOut columns to copy (missing ones skipped).
%
% Every table in targetStats (all channels/planes) gets every prefixed
% column, NaN where no result applies: the app keeps only the columns common
% to all tables when it concatenates them, so a column missing from one
% table would vanish from the Results tab altogether. Existing prefixed
% columns are reset to NaN first, so a re-run fully replaces the previous
% values (idempotent).

newNames = cellfun(@(c) [prefix upper(c(1)) c(2:end)], columns, 'UniformOutput', false);

% 1) every table gets every column, reset to NaN
for k = 1:numel(targetStats)
    T = targetStats{k};
    if ~istable(T) || width(T) == 0
        continue
    end
    for j = 1:numel(newNames)
        T.(newNames{j}) = nan(height(T), 1);
    end
    targetStats{k} = T;
end

% 2) fill matched rows, population by population
population = cellstr(population);
[~, nZ, nT] = size(objectOut);
for iPop = 1:numel(population)
    targetStats = fillPopulation(targetStats, objectOut, population{iPop}, channel(iPop), ...
        idColumn, columns, newNames, nZ, nT);
end
end


% =========================================================================
function targetStats = fillPopulation(targetStats, objectOut, population, channel, idColumn, columns, newNames, nZ, nT)
if channel > size(targetStats, 1)
    return
end
for iT = 1:nT
    for iZ = 1:nZ
        O = objectOut{1, iZ, iT};
        if ~istable(O) || isempty(O) || ~ismember('population', O.Properties.VariableNames)
            continue
        end
        O = O(strcmp(O.population, population), :);
        tz = min(iZ, size(targetStats, 2));
        tt = min(iT, size(targetStats, 3));
        T = targetStats{channel, tz, tt};
        if isempty(O) || ~istable(T) || isempty(T) ...
                || ~all(ismember({'cellID', idColumn}, T.Properties.VariableNames))
            continue
        end
        [tf, loc] = ismember([double(T.cellID) double(T.(idColumn))], ...
            [double(O.cellID) double(O.organelleID)], 'rows');
        for j = 1:numel(columns)
            if ismember(columns{j}, O.Properties.VariableNames)
                v = T.(newNames{j});
                v(tf) = double(O.(columns{j})(loc(tf)));
                T.(newNames{j}) = v;
            end
        end
        targetStats{channel, tz, tt} = T;
    end
end
end
