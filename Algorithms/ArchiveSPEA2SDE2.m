classdef ArchiveSPEA2SDE2 < SPEA2SDE
% <2025> <multi/many> <real/integer/label/binary/permutation> <dynamic>
% Archive-based SPEA2 with shift-based density estimation

    methods
        function main(Algorithm,Problem)
            %% 1. 初始化
            Population = Problem.Initialization();
            Fitness    = CalFitness(Population.objs);
            
            % 初始化外部存档
            ExternalArchive = [];
            StageBestHistory = {};
            CurrentStage = getCurrentStageSignature(Problem,300);
            CurrentStageInfo = getStageInfo(Problem,Population,CurrentStage);
            StageHistorySaved = false;

            %% 2. 主循环
            try
                while Algorithm.NotTerminated(Population)
                    MatingPool = TournamentSelection(2,Problem.N,Fitness);
                    Offspring  = OperatorGA(Problem,Population(MatingPool));
                    
                    % 按 D_MIXED_SEQ.Evaluation 中的 t 逻辑判断是否进入新阶段
                    NextStage = getCurrentStageSignature(Problem,300);
                    EnvironmentChanged = ~isequal(NextStage,CurrentStage);
                    
                    if EnvironmentChanged
                        % ==========================================================
                        %                环境变化处理 (基于存档的热启动)
                        % ==========================================================
                        
                        BestSolutions = getBestSolutions(ExternalArchive,Population);
                        NextStageInfo = getStageInfo(Problem,Offspring,NextStage);
                        if ~isempty(BestSolutions)
                            StageBestHistory{end+1} = packStageBest(BestSolutions,CurrentStageInfo,Problem.FE); %#ok<AGROW>
                        end

                        % 3. 按变化类型继承上一阶段的 BestSolutions
                        Population = rebuildPopulationAtChange(Problem,BestSolutions,CurrentStageInfo,NextStageInfo);
                        Fitness    = CalFitness(Population.objs);
                        
                        % 5. 清空存档
                        ExternalArchive = [];
                        CurrentStageInfo = NextStageInfo;
                        CurrentStage = NextStage;
                        
                    else
                        % ==========================================================
                        %                稳态进化
                        % ==========================================================
                        
                        [Population,Fitness] = EnvironmentalSelection([Population,Offspring],Problem.N);
                        
                        % 将子代累积到存档中
                        ExternalArchive = [ExternalArchive, Offspring];
                    end
                end
            catch err
                if ~strcmp(err.identifier,'PlatEMO:Termination')
                    rethrow(err);
                end
            end

            if ~StageHistorySaved
                FinalBestSolutions = getBestSolutions(ExternalArchive,Population);
                if ~isempty(FinalBestSolutions)
                    StageBestHistory{end+1} = packStageBest(FinalBestSolutions,CurrentStageInfo,Problem.FE); %#ok<AGROW>
                end
                saveStageBestHistory(StageBestHistory,Algorithm,Problem);
                StageHistorySaved = true;
            end
        end
    end

end

function BestSolutions = getBestSolutions(ExternalArchive,Population)
    if ~isempty(ExternalArchive)
        Source = ExternalArchive;
    else
        Source = Population;
    end

    if isempty(Source)
        BestSolutions = SOLUTION.empty();
        return;
    end

    BestSolutions = Source.best;
end

function StageInfo = getStageInfo(Problem,Population,StageSignature)
    StageInfo = struct();
    StageInfo.index = [];
    StageInfo.M     = size(Population.objs,2);
    StageInfo.active = [];
    StageInfo.frozen = [];

    if nargin < 3 || isempty(StageSignature)
        return;
    end

    if isfield(StageSignature,'objStage') && ~isempty(StageSignature.objStage)
        StageInfo.index = StageSignature.objStage;
        if isprop(Problem,'ObjGroups') && ~isempty(Problem.ObjGroups) ...
                && StageSignature.objStage >= 1 && StageSignature.objStage <= length(Problem.ObjGroups)
            StageInfo.active = Problem.ObjGroups{StageSignature.objStage};
            StageInfo.M      = length(StageInfo.active);
        end
    end
    if isfield(StageSignature,'frozenStage') && ~isempty(StageSignature.frozenStage) ...
            && isprop(Problem,'FrozenGroups') && ~isempty(Problem.FrozenGroups) ...
            && StageSignature.frozenStage >= 1 && StageSignature.frozenStage <= length(Problem.FrozenGroups)
        StageInfo.frozen = Problem.FrozenGroups{StageSignature.frozenStage};
    end
end

function Packed = packStageBest(BestSolutions,StageInfo,FE)
    Packed = struct();
    Packed.stageIndex = StageInfo.index;
    Packed.M          = StageInfo.M;
    Packed.active     = StageInfo.active;
    Packed.frozen     = StageInfo.frozen;
    Packed.FE         = FE;
    Packed.decs       = BestSolutions.decs;
    Packed.objs       = BestSolutions.objs;
    Packed.cons       = BestSolutions.cons;
end

function saveStageBestHistory(StageBestHistory,Algorithm,Problem)
    folder = fullfile('Data',class(Algorithm));
    [~,~]  = mkdir(folder);
    timestamp = round(posixtime(datetime('now'))*1e6);
    file   = fullfile(folder,sprintf('StageBestHistory_%s_M%d_D%d_%d.mat',...
        class(Problem),Problem.M,Problem.D,timestamp));
    save(file,'StageBestHistory');
end

function Population = rebuildPopulationAtChange(Problem,BestSolutions,OldStageInfo,NewStageInfo)
    Inherited = inheritBestSolutions(Problem,BestSolutions,OldStageInfo,NewStageInfo);
    Population = fillPopulationAfterChange(Problem,Inherited);
end

function Inherited = inheritBestSolutions(Problem,BestSolutions,OldStageInfo,NewStageInfo)
    if isempty(BestSolutions)
        Inherited = SOLUTION.empty();
        return;
    end

    objIncrease   = NewStageInfo.M > OldStageInfo.M;
    objDecrease   = NewStageInfo.M < OldStageInfo.M;
    freezeAdded   = ~isempty(setdiff(NewStageInfo.frozen,OldStageInfo.frozen));
    needReeval    = objIncrease || freezeAdded;

    if objDecrease && canProjectObjectives(OldStageInfo,NewStageInfo)
        Inherited = projectBestSolutions(BestSolutions,OldStageInfo,NewStageInfo);
    else
        Inherited = BestSolutions;
        if ~isequal(OldStageInfo.active,NewStageInfo.active)
            needReeval = true;
        end
    end

    Inherited = downsampleSolutions(Inherited,Problem.N);

    if needReeval
        Inherited = Problem.Evaluation(Inherited.decs);
    end
end

function tf = canProjectObjectives(OldStageInfo,NewStageInfo)
    if isempty(OldStageInfo.active) || isempty(NewStageInfo.active)
        tf = false;
        return;
    end
    tf = all(ismember(NewStageInfo.active,OldStageInfo.active));
end

function Inherited = projectBestSolutions(BestSolutions,OldStageInfo,NewStageInfo)
    [tf,loc] = ismember(NewStageInfo.active,OldStageInfo.active);
    if ~all(tf)
        Inherited = SOLUTION.empty();
        return;
    end

    BestDec      = BestSolutions.decs;
    BestObj      = BestSolutions.objs;
    BestCon      = BestSolutions.cons;
    ProjectedObj = BestObj(:,loc);
    Inherited    = SOLUTION(BestDec,ProjectedObj,BestCon);
    Inherited    = Inherited.best;
end

function Population = fillPopulationAfterChange(Problem,Inherited)
    if isempty(Inherited)
        Population = Problem.Initialization();
        return;
    end

    if length(Inherited) >= Problem.N
        SelectIdx  = randperm(length(Inherited),Problem.N);
        Population = Inherited(SelectIdx);
        return;
    end

    NumToRandom = Problem.N - length(Inherited);
    RandomPop   = Problem.Initialization(NumToRandom);
    Population  = [Inherited,RandomPop];
end

function Solutions = downsampleSolutions(Solutions,N)
    if isempty(Solutions) || length(Solutions) <= N
        return;
    end
    SelectIdx = randperm(length(Solutions),N);
    Solutions = Solutions(SelectIdx);
end

function signature = getCurrentStageSignature(Problem,FirstStageGens)
    if ~(isprop(Problem,'taut') && isprop(Problem,'lastT'))
        signature = [];
        return;
    end

    CurrentGen = floor(Problem.FE / Problem.N);
    if CurrentGen < FirstStageGens
        t = 0;
    else
        t = 1 + floor((CurrentGen - FirstStageGens) / Problem.taut);
    end

    signature = struct('objStage',t,'frozenStage',t);
    if isprop(Problem,'ObjGroups') && ~isempty(Problem.ObjGroups)
        signature.objStage = min(t + 1,length(Problem.ObjGroups));
    end
    if isprop(Problem,'FrozenGroups') && ~isempty(Problem.FrozenGroups)
        signature.frozenStage = min(t + 1,length(Problem.FrozenGroups));
    end
end
