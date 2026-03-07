function changed = Changed(Problem,Population,varargin)
% Detect stage change by schedule only (no extra evaluations)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2025 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    changed = false;

    % For scheduled dynamic problems (e.g., D_MIXED_SEQ), avoid probing via
    % re-evaluation, which may alter internal state and FE.
    if isprop(Problem,'taut') && isprop(Problem,'lastT')
        if nargin >= 3
            FirstStageGens = max(0,round(varargin{1}));
        else
            FirstStageGens = 300;
        end
        % Existing solutions were evaluated before the current FE counter was
        % incremented, so derive the effective generation from FE-1.
        CurrentGen = floor(max(0,Problem.FE-1) / Problem.N);
        if CurrentGen < FirstStageGens
            tNow = 0;
        else
            tNow = 1 + floor((CurrentGen - FirstStageGens) / Problem.taut);
        end

        tLast = Problem.lastT;

        % Compare effective objective-stage index
        if isprop(Problem,'ObjGroups') && ~isempty(Problem.ObjGroups)
            objNow  = min(tNow  + 1, length(Problem.ObjGroups));
            objLast = min(tLast + 1, length(Problem.ObjGroups));
        else
            objNow  = tNow;
            objLast = tLast;
        end

        % Compare effective frozen-stage index
        if isprop(Problem,'FrozenGroups') && ~isempty(Problem.FrozenGroups)
            frzNow  = min(tNow  + 1, length(Problem.FrozenGroups));
            frzLast = min(tLast + 1, length(Problem.FrozenGroups));
        else
            frzNow  = tNow;
            frzLast = tLast;
        end

        changed = (objNow ~= objLast) || (frzNow ~= frzLast);
    end
end
