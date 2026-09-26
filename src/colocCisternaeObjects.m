function [objStats, info] = colocCisternaeObjects(cisternaeStats, cisternaeClass, streamsThreshold)
%COLOCCISTERNAEOBJECTS  ER cisternae as object tables for colocObjectOverlap,
% optionally split into ordinary cisternae and fast-moving "streams".
%
%   [objStats, info] = colocCisternaeObjects(cisternaeStats, cisternaeClass, streamsThreshold)
%
% INPUTS
%   cisternaeStats   - {nC x nZ x nT} cisternae tables from
%                      analyzerFeatureAnalysis(..., 'cisternae', ...): needs
%                      .cisternaePixelIdxList, .cisternaeID and .cellID;
%                      .cisternaeSpeedMax (added when optical flow was run)
%                      for the streams split.
%   cisternaeClass   - 'all' | 'ordinary' | 'streams'.
%   streamsThreshold - cisternaeSpeedMax at/above which a cisterna counts as
%                      a stream (same units and meaning as the Results/Plot
%                      streams filter). <= 0 means "no streams
%                      classification": every cisterna is ordinary.
%
% OUTPUTS
%   objStats - same-size cell array of tables with the columns
%              colocObjectOverlap expects (organellePixelIdxList,
%              organelleID = cisternaeID, cellID) plus isStream and, when
%              available, cisternaeSpeedMax. Filtered to the requested class.
%   info     - struct: .hasSpeed (any table had cisternaeSpeedMax),
%              .nTotal, .nKept (cisternae before/after the class filter),
%              .message (non-empty when the request could not be met, e.g.
%              'streams' without optical flow or with threshold <= 0).

cisternaeClass = lower(char(cisternaeClass));
if ~ismember(cisternaeClass, {'all','ordinary','streams'})
    error('colocCisternaeObjects:class', 'cisternaeClass must be ''all'', ''ordinary'' or ''streams''.');
end

objStats = cell(size(cisternaeStats));
info = struct('hasSpeed', false, 'nTotal', 0, 'nKept', 0, 'message', '');

for k = 1:numel(cisternaeStats)
    S = cisternaeStats{k};
    if ~istable(S) || isempty(S) || ~ismember('cisternaePixelIdxList', S.Properties.VariableNames)
        continue
    end
    n = height(S);
    T = table();
    T.organellePixelIdxList = S.cisternaePixelIdxList;
    if ismember('cisternaeID', S.Properties.VariableNames)
        T.organelleID = double(S.cisternaeID);
    else
        T.organelleID = (1:n)';
    end
    if ismember('cellID', S.Properties.VariableNames)
        T.cellID = double(S.cellID);
    else
        T.cellID = ones(n, 1);
    end
    hasSpeed = ismember('cisternaeSpeedMax', S.Properties.VariableNames);
    if hasSpeed
        T.cisternaeSpeedMax = double(S.cisternaeSpeedMax);
        info.hasSpeed = true;
        T.isStream = streamsThreshold > 0 & T.cisternaeSpeedMax >= streamsThreshold;
    else
        T.isStream = false(n, 1);
    end

    switch cisternaeClass
        case 'ordinary', T = T(~T.isStream, :);
        case 'streams',  T = T(T.isStream, :);
    end
    info.nTotal = info.nTotal + n;
    info.nKept  = info.nKept + height(T);
    objStats{k} = T;
end

if strcmp(cisternaeClass, 'streams')
    if ~info.hasSpeed
        info.message = 'No cisternaeSpeedMax -- run ER optical flow before cisternae analysis to classify streams.';
    elseif streamsThreshold <= 0
        info.message = 'Streams threshold is 0 (no streams classification) -- set it in the ER Analysis cisternae panel.';
    end
end
end
