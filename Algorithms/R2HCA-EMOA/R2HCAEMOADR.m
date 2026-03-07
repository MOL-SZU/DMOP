classdef R2HCAEMOADR < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation> <dynamic>
% Dynamic R2HCA-EMOA with stage-pool transfer for D_MIXED_SEQ
% numVec       --- 100  --- Direction vector number
% pDirect      --- 0.4  --- Ratio of direct inherited solutions
% alphaBase    --- 0.15 --- Base ratio of FE budget used for re-evaluation at change
% firstStageGens --- 300 --- Warm-up generations before dynamic switching
% poolChunk    --- 200  --- StagePool lightweight pruning interval (offspring count)
% poolCapFac   --- 8    --- StagePool cap factor, cap = poolCapFac * N
% divIncRatio  --- 0.2  --- Minimum random diversity ratio when objective increases

%--------------------------------------------------------------------------
% Core idea:
% 1) Maintain one stagePool for the current stage.
% 2) When the environment changes, finalize the previous stagePool.
% 3) Use only that finalized stagePool to rebuild the next-stage population.
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            params = parseParameters(Algorithm,Problem);
            Algorithm.save = sign(Algorithm.save)*inf;

            Population   = Problem.Initialization();
            currentState = snapshotState(Problem,params.firstStageGens);
            stagePool    = initializeStagePool(Population);
            allStagePools = {};
            stageSaved   = false;
            forceChange  = false;
            r2State      = buildR2State(Problem,Population,params.numVec);

            while Algorithm.NotTerminated(Population)
                if forceChange || Changed(Problem,Population,params.firstStageGens)
                    inheritPool = finalizeStagePool(stagePool);
                    allStagePools{end+1} = inheritPool; %#ok<AGROW>

                    nextState  = inferNextState(Problem,currentState,params.firstStageGens);
                    Population = rebuildFromPreviousStage(Problem,inheritPool,currentState,nextState,params);
                    currentState = snapshotState(Problem,params.firstStageGens);
                    stagePool  = initializeStagePool(Population);
                    r2State    = buildR2State(Problem,Population,params.numVec);
                    forceChange = false;

                    if Problem.FE >= Problem.maxFE
                        if ~stageSaved
                            allStagePools{end+1} = finalizeStagePool(stagePool); %#ok<AGROW>
                            saveStagePools(allStagePools,Algorithm,Problem);
                            stageSaved = true;
                        end
                        continue;
                    end
                end

                [Population,stagePool,r2State,forceChange] = ...
                    evolveOneGeneration(Problem,Population,stagePool,r2State,params);

                if Problem.FE >= Problem.maxFE && ~stageSaved
                    allStagePools{end+1} = finalizeStagePool(stagePool); %#ok<AGROW>
                    saveStagePools(allStagePools,Algorithm,Problem);
                    stageSaved = true;
                end
            end
        end
    end
end

function params = parseParameters(Algorithm,Problem)
    [params.numVec,params.pDirect,params.alphaBase,params.firstStageGens,...
        params.poolChunk,params.poolCapFac,params.divIncRatio] = ...
        Algorithm.ParameterSet(100,0.4,0.15,300,200,8,0.2);

    params.pDirect        = min(max(params.pDirect,0),1);
    params.alphaBase      = max(params.alphaBase,0);
    params.firstStageGens = max(0,round(params.firstStageGens));
    params.poolChunk      = max(1,round(params.poolChunk));
    params.poolCap        = max(Problem.N,round(params.poolCapFac*Problem.N));
    params.divIncRatio    = min(max(params.divIncRatio,0),0.8);
end

function stagePool = initializeStagePool(Population)
    stagePool.decs    = Population.decs;
    stagePool.objs    = Population.objs;
    stagePool.cons    = Population.cons;
    stagePool.pending = 0;
end

function stagePool = appendToStagePool(stagePool,Offspring,poolChunk,poolCap)
    stagePool.decs    = [stagePool.decs;Offspring.decs];
    stagePool.objs    = [stagePool.objs;Offspring.objs];
    stagePool.cons    = [stagePool.cons;Offspring.cons];
    stagePool.pending = stagePool.pending + size(Offspring.decs,1);
    if stagePool.pending >= poolChunk
        stagePool = compactStagePool(stagePool,poolCap);
    end
end

function inheritPool = finalizeStagePool(stagePool)
    inheritPool = pruneStagePool(stagePool);
end

function state = snapshotState(Problem,firstStageGens)
    state = struct('stage',-1,'active',1:Problem.M,'frozen',[],'M',Problem.M);

    if isprop(Problem,'lastT') && Problem.lastT >= 0
        state.stage = Problem.lastT;
    end
    if isprop(Problem,'CurrentIndices') && ~isempty(Problem.CurrentIndices)
        state.active = Problem.CurrentIndices(:)';
        state.M      = length(state.active);
    end
    if isprop(Problem,'FrozenGroups') && isprop(Problem,'lastT') && ~isempty(Problem.FrozenGroups)
        idx = min(max(Problem.lastT+1,1),length(Problem.FrozenGroups));
        state.frozen = unique(Problem.FrozenGroups{idx});
    end

    if state.stage >= 0 && ~isempty(state.active)
        return;
    end

    state = inferStateBySchedule(Problem,firstStageGens);
end

function state = inferNextState(Problem,currentState,firstStageGens)
    state = snapshotState(Problem,firstStageGens);
    if state.stage > currentState.stage || ~isequal(state.active,currentState.active) || ~isequal(state.frozen,currentState.frozen)
        return;
    end

    if isprop(Problem,'ObjGroups') && ~isempty(Problem.ObjGroups)
        stageIdx = min(currentState.stage+2,length(Problem.ObjGroups));
        state.active = Problem.ObjGroups{stageIdx};
        state.M      = length(state.active);
        state.stage  = currentState.stage + 1;
    end
    if isprop(Problem,'FrozenGroups') && ~isempty(Problem.FrozenGroups)
        frozenIdx = min(currentState.stage+2,length(Problem.FrozenGroups));
        state.frozen = unique(Problem.FrozenGroups{frozenIdx});
    end
end

function state = inferStateBySchedule(Problem,firstStageGens)
    state = struct('stage',-1,'active',1:Problem.M,'frozen',[],'M',Problem.M);

    if ~(isprop(Problem,'taut') && isprop(Problem,'ObjGroups'))
        return;
    end

    currentGen = floor(max(0,Problem.FE-1) / Problem.N);
    if currentGen < firstStageGens
        tNow = 0;
    else
        tNow = 1 + floor((currentGen - firstStageGens) / Problem.taut);
    end
    state.stage = tNow;

    if ~isempty(Problem.ObjGroups)
        objIdx = min(max(tNow+1,1),length(Problem.ObjGroups));
        state.active = Problem.ObjGroups{objIdx};
        state.M      = length(state.active);
    end
    if isprop(Problem,'FrozenGroups') && ~isempty(Problem.FrozenGroups)
        frozenIdx = min(max(tNow+1,1),length(Problem.FrozenGroups));
        state.frozen = unique(Problem.FrozenGroups{frozenIdx});
    end
end

function Population = rebuildFromPreviousStage(Problem,inheritPool,oldState,nextState,params)
    alpha      = adaptiveAlpha(params.alphaBase,oldState,nextState);
    Population = inheritCandidates(Problem,inheritPool,oldState,nextState,params.pDirect,alpha);

    if nextState.M > oldState.M
        reserveRandom = min(Problem.N,max(1,round(params.divIncRatio*Problem.N)));
        keepN         = max(0,Problem.N-reserveRandom);
        if length(Population) > keepN
            Population = Population(1:keepN);
        end
    end

    Population = fillPopulation(Problem,Population,Problem.N);
end

function Population = inheritCandidates(Problem,stagePool,oldState,newState,pDirect,alpha)
    Population = SOLUTION.empty();
    if isemptyStagePool(stagePool)
        return;
    end

    stageDec = stagePool.decs;
    stageObj = stagePool.objs;
    stageCon = stagePool.cons;
    candN    = size(stageDec,1);
    oldOrder = qualityOrder(stageObj,stageCon);

    [projectedObj,projectedOrder,canProject] = projectObjectives(stageObj,stageCon,oldState,newState,oldOrder);

    objIncrease = newState.M > oldState.M;
    objDecrease = newState.M < oldState.M;
    frozenAdded = ~isempty(setdiff(newState.frozen,oldState.frozen));
    needReeval  = objIncrease || frozenAdded || ~canProject;

    directIDs = [];
    if ~needReeval && ~isempty(projectedObj)
        if objDecrease
            frontNo  = NDSort(projectedObj,1);
            directIDs = find(frontNo == 1);
            directIDs = intersect(projectedOrder,directIDs,'stable');
            if isempty(directIDs)
                directIDs = projectedOrder;
            end
            targetDirect = min(Problem.N,length(directIDs));
        else
            directIDs    = projectedOrder;
            targetDirect = min(round(pDirect*Problem.N),candN);
        end
        directIDs = directIDs(1:min(targetDirect,length(directIDs)));
        Population = buildInheritedSolutions(Problem,stageDec,projectedObj,stageCon,directIDs);
    end

    needN = Problem.N - length(Population);
    if needN <= 0
        Population = Population(1:Problem.N);
        return;
    end

    remainFE   = max(0,Problem.maxFE-Problem.FE);
    evalBudget = min([needN,remainFE,max(1,round(alpha*Problem.N)),candN]);
    if evalBudget <= 0
        return;
    end

    if objDecrease && ~isempty(projectedObj)
        reeOrder = projectedOrder;
    else
        reeOrder = oldOrder;
    end
    if ~isempty(directIDs)
        reeOrder = setdiff(reeOrder,directIDs,'stable');
    end
    if isempty(reeOrder)
        reeOrder = oldOrder;
    end

    reeIDs = reeOrder(1:min(evalBudget,length(reeOrder)));
    reePop = Problem.Evaluation(stageDec(reeIDs,:));
    if isempty(Population)
        Population = reePop;
    else
        Population = [Population,reePop];
    end
end

function [projectedObj,projectedOrder,canProject] = projectObjectives(stageObj,stageCon,oldState,newState,defaultOrder)
    projectedObj   = [];
    projectedOrder = defaultOrder;
    canProject     = all(ismember(newState.active,oldState.active));
    if ~canProject
        return;
    end

    [tf,loc] = ismember(newState.active,oldState.active);
    if ~all(tf)
        canProject = false;
        return;
    end

    projectedObj   = stageObj(:,loc);
    projectedOrder = qualityOrder(projectedObj,stageCon);
end

function Population = buildInheritedSolutions(Problem,stageDec,projectedObj,stageCon,ids)
    Population = SOLUTION.empty();
    if isempty(ids)
        return;
    end
    ids    = ids(:)';
    PopDec = stageDec(ids,:);
    PopObj = projectedObj(ids,:);
    PopCon = stageCon(ids,:);
    PopAdd = repmat(Problem.FE,size(PopDec,1),1);
    Population = SOLUTION(PopDec,PopObj,PopCon,PopAdd);
end

function Population = fillPopulation(Problem,Population,targetN)
    if length(Population) > targetN
        Population = Population(1:targetN);
        return;
    end

    remainFE = max(0,Problem.maxFE-Problem.FE);
    addN     = min(targetN-length(Population),remainFE);
    if addN > 0
        Population = [Population,Problem.Initialization(addN)];
    end
    if length(Population) > targetN
        Population = Population(1:targetN);
    end
end

function [Population,stagePool,r2State,forceChange] = evolveOneGeneration(Problem,Population,stagePool,r2State,params)
    forceChange = false;
    for i = 1 : Problem.N
        Offspring = OperatorGAhalf(Problem,Population(randperm(length(Population),2)));
        if size(Offspring.objs,2) ~= size(Population.objs,2)
            forceChange = true;
            return;
        end

        stagePool = appendToStagePool(stagePool,Offspring,params.poolChunk,params.poolCap);
        merged    = spliceOffspring(Population,Offspring,r2State.worst,Problem.N);
        [Population,r2State.worst,r2State.tensor] = ...
            Reduce(Problem,merged,r2State.W,r2State.worst,r2State.tensor,r2State.r,r2State.row);
    end
end

function Population1 = spliceOffspring(Population,Offspring,worst,N)
    if worst == 1
        Population1 = [Offspring,Population];
    elseif worst == N + 1
        Population1 = [Population,Offspring];
    else
        Population1 = [Population(1:worst-1),Offspring,Population(worst:end)];
    end
end

function stagePoolND = pruneStagePool(stagePool)
    stagePoolND = stagePool;
    if isemptyStagePool(stagePool)
        return;
    end

    frontNo = constrainedFront(stagePool.objs,stagePool.cons);
    ndMask  = frontNo == 1;
    if any(ndMask)
        stagePoolND = subsetStagePool(stagePool,find(ndMask));
    end
end

function r2State = buildR2State(Problem,Population,numVec)
    r2State.r = 1 + 1/getH(Problem.M,Problem.N);
    [r2State.W,r2State.row] = UniformVector(numVec,Problem.M);

    PopObj = alignObjectives(Population.objs,Problem.M);
    if isempty(PopObj)
        PopObj = zeros(0,Problem.M);
    end
    fmax  = max(PopObj,[],1);
    fmin  = min(PopObj,[],1);
    span  = max(fmax-fmin,1e-12);
    PopObj = (PopObj-repmat(fmin,size(PopObj,1),1))./repmat(span,size(PopObj,1),1);

    if size(PopObj,1) < Problem.N
        PopObj = [PopObj;zeros(Problem.N-size(PopObj,1),Problem.M)];
    end
    PopObj = [PopObj;zeros(1,Problem.M)];

    r2State.tensor = InitializeUtilityTensor(Problem,PopObj,r2State.W,r2State.r,r2State.row);
    r2State.worst  = Problem.N + 1;
end

function PopObj = alignObjectives(PopObj,M)
    if isempty(PopObj)
        PopObj = zeros(0,M);
    elseif size(PopObj,2) > M
        PopObj = PopObj(:,1:M);
    elseif size(PopObj,2) < M
        PopObj = [PopObj,zeros(size(PopObj,1),M-size(PopObj,2))];
    end
end

function alpha = adaptiveAlpha(alphaBase,oldState,newState)
    alpha = alphaBase;
    if newState.M > oldState.M
        alpha = alpha + 0.10;
    end
    if ~isequal(oldState.active,newState.active)
        alpha = alpha + 0.05;
    end
    if ~isempty(setdiff(newState.frozen,oldState.frozen))
        alpha = alpha + 0.10;
    end
    alpha = min(max(alpha,0.05),0.50);
end

function order = qualityOrder(PopObj,PopCon)
    if isempty(PopObj)
        order = [];
        return;
    end
    if isempty(PopCon)
        PopCon = zeros(size(PopObj,1),1);
    end

    CV       = sum(max(PopCon,0),2);
    feasible = CV <= 0;
    order    = [];

    if any(feasible)
        ObjF   = PopObj(feasible,:);
        FrontF = NDSort(ObjF,inf);
        CrowF  = CrowdingDistance(ObjF,FrontF);
        idxF   = find(feasible);
        [~,rk] = sortrows([FrontF(:),-CrowF(:)],[1,2]);
        order  = [order;idxF(rk)]; %#ok<AGROW>
    end

    if any(~feasible)
        idxI     = find(~feasible);
        [~,rkI]  = sort(CV(idxI),'ascend');
        order    = [order;idxI(rkI)]; %#ok<AGROW>
    end

    order = order(:)';
end

function stagePool = compactStagePool(stagePool,poolCap)
    if isemptyStagePool(stagePool)
        return;
    end

    decKey = round(stagePool.decs*1e10)/1e10;
    [~,ia] = unique(decKey,'rows','stable');
    stagePool = subsetStagePool(stagePool,ia);
    if size(stagePool.decs,1) > poolCap
        order = qualityOrder(stagePool.objs,stagePool.cons);
        stagePool = subsetStagePool(stagePool,order(1:poolCap));
    end
    stagePool.pending = 0;
end

function frontNo = constrainedFront(PopObj,PopCon)
    if isempty(PopObj)
        frontNo = [];
        return;
    end

    if isempty(PopCon)
        PopCon = zeros(size(PopObj,1),1);
    end
    CV       = sum(max(PopCon,0),2);
    feasible = CV <= 0;
    frontNo  = inf(size(PopObj,1),1);
    if any(feasible)
        frontNo(feasible) = NDSort(PopObj(feasible,:),inf);
    end
end

function stagePool = subsetStagePool(stagePool,idx)
    idx = idx(:);
    stagePool.decs = stagePool.decs(idx,:);
    stagePool.objs = stagePool.objs(idx,:);
    if isempty(stagePool.cons)
        stagePool.cons = zeros(length(idx),0);
    else
        stagePool.cons = stagePool.cons(idx,:);
    end
end

function tf = isemptyStagePool(stagePool)
    tf = isempty(stagePool) || ~isfield(stagePool,'decs') || isempty(stagePool.decs);
end

function saveStagePools(allStagePools,Algorithm,Problem)
    folder = fullfile('Data',class(Algorithm));
    [~,~]  = mkdir(folder);
    file   = fullfile(folder,sprintf('StagePools_%s_M%d_D%d_%d.mat',...
        class(Problem),Problem.M,Problem.D,round(now*1e8)));
    StagePools = allStagePools; %#ok<NASGU>
    save(file,'StagePools');
end
