classdef LEC < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation> <dynamic>
% Learning to Expand and Contract Pareto Sets
% lambda --- 20 --- Parameter for contraction sampling (popsize/20)
%
% 修正版：
%   1. 使用 stage-signature 先验环境检测，与保守版 DTAEA 一致；
%   2. 不再只依赖 Problem.M 或 Offspring.objs 维度变化；
%   3. 不允许 UniformPoint 改写 Problem.N；
%   4. 固定 PopSize = Problem.N，避免阶段判断 FE / Problem.N 错乱；
%   5. 目标数增加时保留 LEC 原 PSExpansion；
%   6. 目标数减少时保留 LEC 原 PSContraction；
%   7. 目标数不变但 stage 变化时，只重评 Population；
%   8. 增加 Population 规模修复、邻域索引保护、DensitySelection 越界保护。

%------------------------------- Reference --------------------------------
% G. Ruan, L. L. Minku, S. Menzel, B. Sendhoff, and X. Yao,
% "Learning to Expand/Contract Pareto Sets in Dynamic Multiobjective
% Optimization With a Changing Number of Objectives,"
% IEEE Transactions on Evolutionary Computation, 2025.
%--------------------------------------------------------------------------

    properties
        lambda = 20;    % Contraction sampling parameter
        W;              % Weight vectors
        Neighbors;      % Neighborhood index
        T = 20;         % Neighborhood size
        delta = 0.9;    % Probability of choosing parents from neighborhood
    end

    methods
        function main(Algorithm, Problem)

            %% Parameter setting
            Algorithm.lambda = Algorithm.ParameterSet(20);
            % Paper setting: LEC uses the same random mating optimizer as KTDMOEA.
            Algorithm.delta = 0;

            %% 固定种群规模
            % 关键修复：
            %   不要让 UniformPoint 改写 Problem.N。
            %   Problem.N 应保持实验设置中的 popsize，例如 300。
            PopSize = Problem.N;

            %% Initialization: decomposition strategy setup
            [Algorithm.W, ~] = Algorithm.UniformPoint(PopSize, Problem.M);

            Algorithm.T = min(size(Algorithm.W, 1), 20);
            Algorithm.Neighbors = pdist2(Algorithm.W, Algorithm.W, 'euclidean');
            [~, Algorithm.Neighbors] = sort(Algorithm.Neighbors, 2);
            Algorithm.Neighbors = Algorithm.Neighbors(:, 1:Algorithm.T);

            %% Initialize Population
            Population = Problem.Initialization();
            Population = EnsurePopulationSize_LEC(Algorithm, Problem, Population, PopSize);

            %% Stage-signature prior change detection
            if isprop(Problem, 'FirstStageGens') && ~isempty(Problem.FirstStageGens)
                firstGens = Problem.FirstStageGens;
            elseif isprop(Problem, 'taut') && ~isempty(Problem.taut)
                firstGens = Problem.taut;
            else
                error('LEC requires Problem.FirstStageGens or Problem.taut to determine the first change generation.');
            end

            CurrentStage     = getCurrentStageSignature_LEC(Problem, firstGens);
            CurrentStageInfo = getStageInfo_LEC(Problem, Population, CurrentStage);

            lastM = CurrentStageInfo.M;

            %% Optimization Loop
            while Algorithm.NotTerminated(Population)

                % =========================================================
                % Step 0: Prior environment detection by stage signature
                % =========================================================
                NextStage = getCurrentStageSignature_LEC(Problem, firstGens);
                EnvironmentChanged = ~isequal(NextStage, CurrentStage);

                if EnvironmentChanged

                    OldStageInfo = CurrentStageInfo;
                    NewStageInfo = getStageInfo_LEC(Problem, Population, NextStage);

                    fprintf('[LEC Change Detected] FE=%d | OldStage=%s | NewStage=%s | OldM=%d | NewM=%d\n', ...
                        Problem.FE, ...
                        stageSigToString_LEC(CurrentStage), ...
                        stageSigToString_LEC(NextStage), ...
                        OldStageInfo.M, ...
                        NewStageInfo.M);

                    % -----------------------------------------------------
                    % 更新权重向量，但不要修改 Problem.N
                    % -----------------------------------------------------
                    [Algorithm.W, ~] = Algorithm.UniformPoint(PopSize, NewStageInfo.M);

                    Algorithm.T = min(size(Algorithm.W, 1), 20);
                    Algorithm.Neighbors = pdist2(Algorithm.W, Algorithm.W, 'euclidean');
                    [~, Algorithm.Neighbors] = sort(Algorithm.Neighbors, 2);
                    Algorithm.Neighbors = Algorithm.Neighbors(:, 1:Algorithm.T);

                    % -----------------------------------------------------
                    % 环境响应
                    % -----------------------------------------------------
                    if NewStageInfo.M > OldStageInfo.M

                        % 目标数增加：
                        % 保留 LEC 原论文 PS expansion 逻辑
                        Population = Algorithm.PSExpansion(Population, Problem, OldStageInfo.M);

                    elseif NewStageInfo.M < OldStageInfo.M

                        % 目标数减少：
                        % 保留 LEC 原论文 PS contraction 逻辑
                        Population = Algorithm.PSContraction(Population, Problem, OldStageInfo.M);

                    else

                        % 目标数不变但 stage 变化：
                        % 保守处理，只重评，不额外扩张/收缩
                        Population = Problem.Evaluation(Population.decs);

                    end

                    % 保证环境响应后的种群规模稳定为 PopSize
                    Population = EnsurePopulationSize_LEC(Algorithm, Problem, Population, PopSize);

                    % 更新当前目标维度，以实际种群为准
                    if ~isempty(Population)
                        NewStageInfo.M = size(Population(1).objs, 2);
                    end

                    CurrentStage     = NextStage;
                    CurrentStageInfo = NewStageInfo;
                    lastM = NewStageInfo.M; %#ok<NASGU>

                    % 环境响应后跳过本轮 reproduction
                    continue;
                end

                % =========================================================
                % Step 1: Mating selection
                % =========================================================
                Population = EnsurePopulationSize_LEC(Algorithm, Problem, Population, PopSize);

                currentM = size(Population(1).objs, 2);
                PopLen   = length(Population);

                OffspringSize = PopSize;
                RemainingFE = RemainingEvaluationsInCurrentStage_LEC(Problem, CurrentStageInfo, firstGens, PopSize);

                if isfinite(RemainingFE) && RemainingFE < PopSize
                    OffspringSize = 2 * floor(RemainingFE / 2);

                    if OffspringSize < 2
                        if RemainingFE > 0
                            Problem.Initialization(RemainingFE);
                        end
                        continue;
                    end
                end

                MatingPool = zeros(1, OffspringSize);

                for i = 1 : OffspringSize

                    if PopLen <= 0
                        error('LEC population is empty before mating selection.');
                    end

                    if i <= size(Algorithm.Neighbors, 1) && rand < Algorithm.delta

                        P = Algorithm.Neighbors(i, randi(Algorithm.T));

                        % 邻域编号来自权重向量，不一定小于 Population 长度
                        if P > PopLen
                            P = randi(PopLen);
                        end

                    else
                        P = randi(PopLen);
                    end

                    MatingPool(i) = P;
                end

                %% Generate Offspring
                Offspring = OperatorGA(Problem, Population(MatingPool));

                %% Population update
                Population = Algorithm.UpdateCA([Population, Offspring], PopSize);

                if ~isequal(getCurrentStageSignature_LEC(Problem, firstGens), CurrentStage)
                    continue;
                end

                Population = EnsurePopulationSize_LEC(Algorithm, Problem, Population, PopSize);
            end
        end


        %% Part 1: PS Expansion
        function NewPop = PSExpansion(Algorithm, OldPop, Problem, oldM)

            % Extract Old Pareto Set
            [FrontNo, ~] = NDSort(OldPop.objs, 1);
            PSt = OldPop(FrontNo == 1);

            if isempty(PSt)
                PSt = OldPop;
            end

            % 1. Learn Candidate Expansion Directions
            Dirs_Cand = [];

            if length(PSt) > 2

                try
                    [Coeff, ~, ~] = pca(PSt.decs);
                    nEig = min(size(Coeff, 2), oldM - 1);

                    if nEig > 0

                        EVs = Coeff(:, 1:nEig);
                        NullBasis = null(EVs');

                        N = Problem.N;
                        nNull = size(NullBasis, 2);

                        if nNull > 0

                            Coefs = 2 * lhsdesign(N, nNull) - 1;
                            Dirs_Cand = Coefs * NullBasis';

                            len = sqrt(sum(Dirs_Cand.^2, 2));
                            Dirs_Cand = Dirs_Cand ./ (len + 1e-10);
                        end
                    end

                catch
                    Dirs_Cand = [];
                end
            end

            % 2. Select Most Promising Directions
            D_exp = [];

            % 重评旧 PS 到新环境
            OldP_NewEnv = Problem.Evaluation(PSt.decs);

            if ~isempty(Dirs_Cand)

                idx = randi(length(PSt));
                x_base = PSt(idx);

                for i = 1 : size(Dirs_Cand, 1)

                    D_vec = Dirs_Cand(i, :);
                    y_dec = Algorithm.GenerateSolution(x_base.decs, D_vec, Problem.lower, Problem.upper);
                    y = Problem.Evaluation(y_dec);

                    dominated = false;

                    for k = 1 : length(OldP_NewEnv)

                        if all(OldP_NewEnv(k).objs <= y.objs) && ...
                                any(OldP_NewEnv(k).objs < y.objs)

                            dominated = true;
                            break;
                        end
                    end

                    if ~dominated
                        D_exp = [D_exp; D_vec];
                    end
                end
            end

            % 3. Expand PS
            if isempty(D_exp)
                NewPop = OldP_NewEnv;
                return;
            end

            N = Problem.N;
            nSelect = min(length(PSt), N);

            [W_old, ~] = Algorithm.UniformPoint(nSelect, oldM);

            Objs = PSt.objs;

            if size(Objs, 2) ~= size(W_old, 2)
                [W_old, ~] = Algorithm.UniformPoint(nSelect, size(Objs, 2));
            end

            [~, RegionIdx] = max(1 - pdist2(Objs, W_old, 'cosine'), [], 2);

            BaseSols = [];
            usedRegions = unique(RegionIdx);

            for i = 1 : length(usedRegions)

                regionMembers = PSt(RegionIdx == usedRegions(i));
                randIdx = randi(length(regionMembers));

                BaseSols = [BaseSols, regionMembers(randIdx)];
            end

            if isempty(BaseSols)
                BaseSols = PSt;
            end

            TransferredDecs = [];
            nDirs = size(D_exp, 1);
            counter = 0;

            while size(TransferredDecs, 1) < N

                counter = counter + 1;

                idxBase = mod(counter - 1, length(BaseSols)) + 1;
                x_i = BaseSols(idxBase).decs;

                idxDir = mod(counter - 1, nDirs) + 1;
                D_vec = D_exp(idxDir, :);

                x_new = Algorithm.GenerateSolution(x_i, D_vec, Problem.lower, Problem.upper);
                TransferredDecs = [TransferredDecs; x_new];
            end

            NewPop = Problem.Evaluation(TransferredDecs(1:min(size(TransferredDecs, 1), N), :));
        end


        %% Part 2: PS Contraction
        function NewPop = PSContraction(Algorithm, OldPop, Problem, oldM)

            [FrontNo, ~] = NDSort(OldPop.objs, 1);
            PSt = OldPop(FrontNo == 1);

            if isempty(PSt)
                PSt = OldPop;
            end

            OldP_NewEnv = Problem.Evaluation(PSt.decs);

            % 1. Learn Candidate Contraction Directions
            C_con = [];

            if length(PSt) >= oldM

                try
                    [Coeff, ~, ~] = pca(PSt.decs);
                    nEig = min(size(Coeff, 2), oldM - 1);

                    if nEig > 0
                        C_con = Coeff(:, 1:nEig)';
                    end

                catch
                    C_con = [];
                end
            end

            % 2. Select Most Promising Directions
            D_con = [];

            if ~isempty(C_con) && ~isempty(OldP_NewEnv)

                idx = randi(length(OldP_NewEnv));
                x_base = OldP_NewEnv(idx);

                for i = 1 : size(C_con, 1)

                    D_vec = C_con(i, :);
                    nGen = floor(Problem.N / max(1, (oldM - 1)));

                    Y_decs = [];

                    for k = 1 : nGen

                        y_dec = Algorithm.GenerateSolution(x_base.decs, D_vec, Problem.lower, Problem.upper);
                        Y_decs = [Y_decs; y_dec];
                    end

                    if isempty(Y_decs)
                        continue;
                    end

                    Y = Problem.Evaluation(Y_decs);

                    promising = false;

                    for k = 1 : length(Y)

                        if all(Y(k).objs <= x_base.objs) && ...
                                any(Y(k).objs < x_base.objs)

                            promising = true;
                            break;
                        end
                    end

                    if promising
                        D_con = [D_con; D_vec];
                    end
                end
            end

            % 3. Contract PS
            if isempty(D_con)
                NewPop = OldP_NewEnv;
                return;
            end

            nBase = floor(Problem.N / Algorithm.lambda);

            if nBase < 1
                nBase = 1;
            end

            perm = randperm(length(PSt));
            P_base = PSt(perm(1:min(length(PSt), nBase)));

            P_g_Decs = [];
            N_con = size(D_con, 1);

            for i = 1 : length(P_base)

                x_i = P_base(i).decs;

                for j = 1 : N_con

                    D_vec = D_con(j, :);

                    Ng = floor(Problem.N / (N_con * length(P_base)));

                    if Ng < 1
                        Ng = 1;
                    end

                    for k = 1 : Ng

                        x_new = Algorithm.GenerateSolution(x_i, D_vec, Problem.lower, Problem.upper);
                        P_g_Decs = [P_g_Decs; x_new];
                    end
                end
            end

            if isempty(P_g_Decs)
                NewPop = OldP_NewEnv;
                return;
            end

            P_g = Problem.Evaluation(P_g_Decs);

            Combined = [OldP_NewEnv, P_g];

            if isempty(Combined)
                NewPop = Problem.Initialization();
                return;
            end

            currentM = size(Combined(1).objs, 2);

            NewPop = Algorithm.DensitySelection(Combined, Problem.N, currentM);
        end


        %% Helper: Generate Solution with Step Size
        function x_new = GenerateSolution(~, x, D_vec, lower, upper)

            D_vec(abs(D_vec) < 1e-10) = 1e-10;

            para = zeros(size(x));

            maskPos = D_vec > 0;
            maskNeg = ~maskPos;

            if isscalar(lower)
                lower = repmat(lower, size(x));
            end

            if isscalar(upper)
                upper = repmat(upper, size(x));
            end

            if size(lower, 1) > 1
                lower = lower(:)';
            end

            if size(upper, 1) > 1
                upper = upper(:)';
            end

            lower = lower(1:length(x));
            upper = upper(1:length(x));

            para(maskPos) = (upper(maskPos) - x(maskPos)) ./ D_vec(maskPos);
            para(maskNeg) = (lower(maskNeg) - x(maskNeg)) ./ D_vec(maskNeg);

            ss = min(para);

            if isempty(ss) || ~isfinite(ss) || ss < 0
                ss = 0;
            end

            x_new = x + ss * rand() * D_vec;
            x_new = max(min(x_new, upper), lower);
        end


        %% Helper: CA Update Mechanism used by KTDMOEA/LEC
        function Population = UpdateCA(~, MixedPop, N)

            if isempty(MixedPop)
                Population = MixedPop;
                return;
            end

            [~, uniqueIdx] = unique(MixedPop.objs, 'rows');
            MixedPop = MixedPop(uniqueIdx);

            [FrontNo, ~] = NDSort(MixedPop.objs, 1);
            Population = MixedPop(FrontNo == 1);

            if isempty(Population)
                Population = MixedPop(1:min(N, length(MixedPop)));
                return;
            end

            if length(Population) > N
                [~, Rank] = sort(CrowdingDistance(Population.objs), 'descend');
                Population = Population(Rank(1:N));
            end
        end


        %% Helper: Density Selection
        function Population = DensitySelection(Algorithm, Population, N, M)

            if isempty(Population)
                return;
            end

            if length(Population) <= N
                return;
            end

            % 使用实际目标维度，防止 M 与 Population.objs 不一致
            M_actual = size(Population(1).objs, 2);

            if nargin < 4 || isempty(M) || M ~= M_actual
                M = M_actual;
            end

            % Non-dominated sorting
            [FrontNo, MaxFNo] = NDSort(Population.objs, Population.cons, N);

            Next = FrontNo < MaxFNo;
            Last = find(FrontNo == MaxFNo);

            nNext = sum(Next);
            nNeed = N - nNext;

            if nNeed <= 0
                idx = find(Next);

                if length(idx) > N
                    idx = idx(1:N);
                end

                Population = Population(idx);
                return;
            end

            if isempty(Last)
                Population = Population(1:min(N, length(Population)));
                return;
            end

            LastSols = Population(Last);
            nLast = length(LastSols);

            if nLast <= nNeed
                Next(Last) = true;
                Population = Population(Next);

                if length(Population) > N
                    Population = Population(1:N);
                end

                return;
            end

            if ~isempty(Algorithm.W) && size(Algorithm.W, 2) == M
                W_current = Algorithm.W;
            else
                [W_current, ~] = Algorithm.UniformPoint(N, M);
            end

            Objs = LastSols.objs;

            Zmin = min(Objs, [], 1);
            Zmax = max(Objs, [], 1);

            ObjsNorm = (Objs - Zmin) ./ (Zmax - Zmin + 1e-10);

            RegionIdx = zeros(1, nLast);

            for i = 1 : nLast

                sol = ObjsNorm(i, :);
                norm_sol = norm(sol);

                if norm_sol == 0
                    norm_sol = 1e-6;
                end

                if size(W_current, 2) ~= length(sol)
                    [W_current, ~] = Algorithm.UniformPoint(N, length(sol));
                end

                cosine = (sol * W_current') / norm_sol;
                [~, RegionIdx(i)] = max(cosine);
            end

            perm = randperm(nLast);
            RegionIdxPerm = RegionIdx(perm);

            [~, unique_idx] = unique(RegionIdxPerm, 'stable');
            unique_idx = unique_idx(:)';

            if length(unique_idx) >= nNeed

                SelectedPermIdx = unique_idx(1:nNeed);

            else

                SelectedPermIdx = unique_idx;

                remaining = setdiff(1:nLast, unique_idx);

                needExtra = nNeed - length(SelectedPermIdx);
                takeNum = min(needExtra, length(remaining));

                if takeNum > 0
                    SelectedPermIdx = [SelectedPermIdx, remaining(1:takeNum)];
                end

                if length(SelectedPermIdx) < nNeed

                    allIdx = 1:nLast;
                    rest = setdiff(allIdx, SelectedPermIdx);

                    needExtra = nNeed - length(SelectedPermIdx);
                    takeNum = min(needExtra, length(rest));

                    if takeNum > 0
                        SelectedPermIdx = [SelectedPermIdx, rest(1:takeNum)];
                    end
                end
            end

            if length(SelectedPermIdx) > nNeed
                SelectedPermIdx = SelectedPermIdx(1:nNeed);
            end

            RealSelectedIdx = perm(SelectedPermIdx);
            Next(Last(RealSelectedIdx)) = true;

            Population = Population(Next);

            if length(Population) > N
                Population = Population(1:N);
            end
        end


        %% Helper: UniformPoint
        function [W, N] = UniformPoint(~, N, M)

            H1 = 1;

            while nchoosek(H1 + M, M - 1) <= N
                H1 = H1 + 1;
            end

            W = nchoosek(1:H1 + M - 1, M - 1) - ...
                repmat(0:M - 2, nchoosek(H1 + M - 1, M - 1), 1) - 1;

            W = ([W, zeros(size(W, 1), 1) + H1] - ...
                [zeros(size(W, 1), 1), W]) / H1;

            if H1 < M

                H2 = 0;

                while nchoosek(H1 + M - 1, M - 1) + ...
                        nchoosek(H2 + M, M - 1) <= N
                    H2 = H2 + 1;
                end

                if H2 > 0

                    W2 = nchoosek(1:H2 + M - 1, M - 1) - ...
                        repmat(0:M - 2, nchoosek(H2 + M - 1, M - 1), 1) - 1;

                    W2 = ([W2, zeros(size(W2, 1), 1) + H2] - ...
                        [zeros(size(W2, 1), 1), W2]) / H2;

                    W = [W; W2 / 2 + 1 / (2 * M)];
                end
            end

            W(W < 1e-6) = 1e-6;
            N = size(W, 1);
        end
    end
end


% =========================================================================
% Stage-signature prior change detection helpers for LEC
% =========================================================================
function signature = getCurrentStageSignature_LEC(Problem, FirstStageGens)

    if ~isprop(Problem, 'taut')
        signature = [];
        return;
    end

    CurrentGen = floor(Problem.FE / Problem.N);

    if CurrentGen < FirstStageGens
        t = 0;
    else
        t = 1 + floor((CurrentGen - FirstStageGens) / Problem.taut);
    end

    signature = struct();
    signature.objStage = t;
    signature.frozenStage = t;

    if isprop(Problem, 'ObjGroups') && ~isempty(Problem.ObjGroups)
        signature.objStage = min(t + 1, length(Problem.ObjGroups));
    end

    if isprop(Problem, 'FrozenGroups') && ~isempty(Problem.FrozenGroups)
        signature.frozenStage = min(t + 1, length(Problem.FrozenGroups));
    end
end


function StageInfo = getStageInfo_LEC(Problem, Population, StageSignature)

    StageInfo = struct();
    StageInfo.index  = [];
    StageInfo.M      = size(Population.objs, 2);
    StageInfo.active = [];
    StageInfo.frozen = [];

    if nargin < 3 || isempty(StageSignature)
        return;
    end

    if isfield(StageSignature, 'objStage') && ~isempty(StageSignature.objStage)

        StageInfo.index = StageSignature.objStage;

        if isprop(Problem, 'ObjGroups') && ~isempty(Problem.ObjGroups) && ...
                StageSignature.objStage >= 1 && ...
                StageSignature.objStage <= length(Problem.ObjGroups)

            StageInfo.active = Problem.ObjGroups{StageSignature.objStage};
            StageInfo.M = length(StageInfo.active);
        end
    end

    if isfield(StageSignature, 'frozenStage') && ~isempty(StageSignature.frozenStage) && ...
            isprop(Problem, 'FrozenGroups') && ~isempty(Problem.FrozenGroups) && ...
            StageSignature.frozenStage >= 1 && ...
            StageSignature.frozenStage <= length(Problem.FrozenGroups)

        StageInfo.frozen = Problem.FrozenGroups{StageSignature.frozenStage};
    end
end


function str = stageSigToString_LEC(sig)

    if isempty(sig)
        str = '[]';
        return;
    end

    if isstruct(sig)

        if isfield(sig, 'objStage') && isfield(sig, 'frozenStage')
            str = sprintf('objStage=%d,frozenStage=%d', ...
                sig.objStage, sig.frozenStage);

        elseif isfield(sig, 'objStage')
            str = sprintf('objStage=%d', sig.objStage);

        else
            str = 'struct';
        end

    else
        str = mat2str(sig);
    end
end


% =========================================================================
% Ensure Population Size
% =========================================================================
function Population = EnsurePopulationSize_LEC(Algorithm, Problem, Population, PopSize)

    if isempty(Population)
        Population = Problem.Initialization(PopSize);
        return;
    end

    if length(Population) > PopSize

        M = size(Population(1).objs, 2);
        Population = Algorithm.DensitySelection(Population, PopSize, M);
        return;
    end

    if length(Population) < PopSize

        Need = PopSize - length(Population);

        try
            RandomPop = Problem.Initialization(Need);
        catch
            Decs = repmat(Problem.lower, Need, 1) + ...
                   rand(Need, Problem.D) .* ...
                   (repmat(Problem.upper, Need, 1) - repmat(Problem.lower, Need, 1));

            RandomPop = Problem.Evaluation(Decs);
        end

        Population = [Population, RandomPop];
    end
end


function RemainingFE = RemainingEvaluationsInCurrentStage_LEC(Problem, CurrentStageInfo, FirstStageGens, PopSize)

    if isempty(CurrentStageInfo.index)
        RemainingFE = inf;
        return;
    end

    if isprop(Problem, 'ObjGroups') && ~isempty(Problem.ObjGroups)
        totalStages = length(Problem.ObjGroups);
    else
        RemainingFE = inf;
        return;
    end

    if CurrentStageInfo.index >= totalStages
        RemainingFE = inf;
        return;
    end

    nextChangeGen = FirstStageGens + (CurrentStageInfo.index - 1) * Problem.taut;
    nextChangeFE  = nextChangeGen * PopSize;
    RemainingFE   = nextChangeFE - Problem.FE;

    if RemainingFE < 0
        RemainingFE = 0;
    end
end
