classdef DTAEA < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation> <dynamic>
% Dynamic Two-Archive Evolutionary Algorithm
%
% 修复版：
%   1. 使用 stage-signature 先验环境检测；
%   2. 环境变化时，在重评 CA / DA 之前保存旧环境最后 CA 的 MHV；
%   3. 算法结束后额外保存最后一个阶段的末代 CA 的 MHV；
%   4. 目标数不变但 stage 变化时，只重评 CA / DA，不额外刷新 DA；
%   5. 保存结果写入全局变量 MHV_TEMP_HISTORY，格式为：
%          [TargetGen, TriggerGen, MHV]

    methods
        function main(Algorithm, Problem)

            %% 1. 初始化
            N = Problem.N;

            Population = Problem.Initialization();
            CA = Population;
            DA = Problem.Initialization();

            % -------------------------------------------------------------
            % 读取初始阶段长度
            % 优先使用 Problem.FirstStageGens，否则使用 Problem.taut
            % -------------------------------------------------------------
            if isprop(Problem, 'FirstStageGens') && ~isempty(Problem.FirstStageGens)
                firstGens = Problem.FirstStageGens;
            elseif isprop(Problem, 'taut') && ~isempty(Problem.taut)
                firstGens = Problem.taut;
            else
                error('DTAEA requires Problem.FirstStageGens or Problem.taut.');
            end

            % -------------------------------------------------------------
            % 初始化阶段签名
            % -------------------------------------------------------------
            CurrentStage = getCurrentStageSignature_DTAEA(Problem, firstGens);
            CurrentStageInfo = getStageInfo_DTAEA(Problem, CA, CurrentStage);

            LastM = CurrentStageInfo.M;
            W = UniformPoint(N, LastM);

            % 全局记录
            global GLOBAL_HISTORY;
            GLOBAL_HISTORY = {};

            global DTAEA_CURRENT_CA;
            DTAEA_CURRENT_CA = CA;
            
            % 退出保护：
            % 即使 while 后面的显式保存因为 PlatEMO 终止机制没有正常执行，
            % 也会在 main 函数退出时尝试保存最后阶段。
            cleanupObj = onCleanup(@() SaveDTAEA_FinalStage_MHV_FromGlobal(Problem));

            %% 2. 优化循环
            while Algorithm.NotTerminated(CA)

                % =========================================================
                % Step 0: 先验环境变化检测
                % =========================================================
                NextStage = getCurrentStageSignature_DTAEA(Problem, firstGens);
                EnvironmentChanged = ~isequal(NextStage, CurrentStage);

                if EnvironmentChanged

                    OldStageInfo = CurrentStageInfo;
                    NewStageInfo = getStageInfo_DTAEA(Problem, CA, NextStage);

                    fprintf('[DTAEA Change Detected] FE=%d | OldStage=%s | NewStage=%s | OldM=%d | NewM=%d\n', ...
                        Problem.FE, ...
                        stageSigToString_DTAEA(CurrentStage), ...
                        stageSigToString_DTAEA(NextStage), ...
                        OldStageInfo.M, ...
                        NewStageInfo.M);

                    % =====================================================
                    % 关键修复：
                    % 在环境响应前，保存旧环境最后 CA 的 MHV
                    % =====================================================
                    SaveDTAEA_MHV_BeforeChange(CA, Problem, OldStageInfo, CurrentStage);

                    % -----------------------------------------------------
                    % 1. 更新权重向量到新目标数
                    % -------------------------
                    % ----------------------------
                    W = UniformPoint(N, NewStageInfo.M);

                    % -----------------------------------------------------
                    % 2. 重新评估 CA / DA 到新环境
                    % -----------------------------------------------------
                    CA = Problem.Evaluation(CA.decs);
                    DA = Problem.Evaluation(DA.decs);

                    % -----------------------------------------------------
                    % 3. 按原 DTAEA 逻辑处理目标数增加 / 减少
                    % -----------------------------------------------------
                    if NewStageInfo.M > OldStageInfo.M

                        % =================================================
                        % Case 1: 目标增加
                        % 原论文逻辑：CA 保留，DA 用 LHS 重置
                        % =================================================
                        DA = LatinHypercubeSampling(Problem, N);

                    elseif NewStageInfo.M < OldStageInfo.M

                        % =================================================
                        % Case 2: 目标减少
                        % 原论文逻辑：
                        %   非支配解进入 CA；
                        %   其余解进入 DA；
                        %   CA 不足用变异填充；
                        %   DA 不足用 LHS 填充。
                        % =================================================

                        [FrontNo, ~] = NDSort(CA.objs, inf);

                        NextCA = CA(FrontNo == 1);
                        MovedToDA = CA(FrontNo > 1);

                        NextDA = MovedToDA;

                        % DA 不足，用 LHS 填充
                        if length(NextDA) < N
                            NextDA = [NextDA, LatinHypercubeSampling(Problem, N - length(NextDA))];
                        end

                        % DA 超出，截断
                        if length(NextDA) > N
                            NextDA = NextDA(1:N);
                        end

                        % CA 极端为空时，直接初始化保护
                        if isempty(NextCA)

                            NextCA = LatinHypercubeSampling(Problem, N);

                        elseif length(NextCA) < N

                            Needed = N - length(NextCA);
                            Density = EstimateDensity(NextCA.objs, W);

                            ParentsToMutate = [];

                            for k = 1 : Needed

                                p1 = randi(length(NextCA));
                                p2 = randi(length(NextCA));

                                if Density(p1) < Density(p2)
                                    BestIdx = p1;
                                elseif Density(p1) > Density(p2)
                                    BestIdx = p2;
                                else
                                    if rand() < 0.5
                                        BestIdx = p1;
                                    else
                                        BestIdx = p2;
                                    end
                                end

                                ParentsToMutate = [ParentsToMutate, NextCA(BestIdx)];
                            end

                            if ~isempty(ParentsToMutate)
                                NewMutants = OperatorGA(Problem, ParentsToMutate, {0, 20, 1, 20});
                                NextCA = [NextCA, NewMutants];
                            end
                        end

                        % CA 超出，截断
                        if length(NextCA) > N
                            NextCA = NextCA(1:N);
                        end

                        CA = NextCA;
                        DA = NextDA;

                    else

                        % =================================================
                        % Case 3: 目标数不变，但 stage 变化
                        %
                        % 保守处理：
                        %   前面已经重评 CA / DA；
                        %   这里不刷新 DA，不变异 CA，不额外重构。
                        % =================================================

                    end

                    % -----------------------------------------------------
                    % 4. 更新阶段状态
                    % -----------------------------------------------------
                    CurrentStage = NextStage;
                    CurrentStageInfo = NewStageInfo;
                    LastM = NewStageInfo.M;

                    DTAEA_CURRENT_CA = CA;

                    % 环境变化后，跳过本轮常规繁殖和更新
                    continue;
                end

                % =========================================================
                % Step 1: 正常记录
                % =========================================================
                GLOBAL_HISTORY{end+1} = CA.objs;
                DTAEA_CURRENT_CA = CA;

                % =========================================================
                % Step 2: 繁殖
                % =========================================================

                % 确保 W 维度匹配当前目标数
                if size(CA.objs, 2) ~= size(W, 2)
                    W = UniformPoint(N, size(CA.objs, 2));
                end

                Objs = CA.objs;
                Zmin = min(Objs, [], 1);
                Zmax = max(Objs, [], 1);
                NormObjs = (Objs - Zmin) ./ (Zmax - Zmin + 1e-6);

                Region = GetAssociation(NormObjs, W);
                Occupied = length(unique(Region));
                I_CA = Occupied / size(W, 1);

                Parents = [];

                for i = 1 : N

                    P1 = CA(randi(length(CA)));

                    if rand() < I_CA
                        P2 = CA(randi(length(CA)));
                    else
                        P2 = DA(randi(length(DA)));
                    end

                    Parents = [Parents, P1, P2];
                end

                Offspring = OperatorGA(Problem, Parents);

                % =========================================================
                % Step 3: 更新 CA / DA
                % =========================================================
                CA = UpdateCA(CA, Offspring, W, N);
                DA = UpdateDA(CA, DA, Offspring, W, N);

                DTAEA_CURRENT_CA = CA;
            end

            % =============================================================
            % 关键修复：
            % 保存最后一个阶段的末代 CA 的 MHV
            %
            % 因为最后阶段后面没有下一次环境变化，
            % 不会再触发 EnvironmentChanged。
            % =============================================================
            SaveDTAEA_FinalStage_MHV(CA, Problem);
            % 结束时同步一次
            DTAEA_CURRENT_CA = CA;

            % 可选保存
            try
                Data = GLOBAL_HISTORY;
                save('DTAEA_Strict_Result.mat', 'Data');
            catch
            end
        end
    end
end


% =========================================================================
% 阶段签名检测函数
% =========================================================================
function signature = getCurrentStageSignature_DTAEA(Problem, FirstStageGens)

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


% =========================================================================
% 根据阶段签名读取阶段信息
% =========================================================================
function StageInfo = getStageInfo_DTAEA(Problem, Population, StageSignature)

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


% =========================================================================
% 阶段签名转字符串，仅用于打印
% =========================================================================
function str = stageSigToString_DTAEA(sig)

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
% DTAEA 内部保存：环境变化前最后 CA 的 MHV
% =========================================================================
% =========================================================================
% DTAEA 内部保存：环境变化前最后 CA 的 MHV
% =========================================================================
% =========================================================================
% DTAEA 内部保存：环境变化前最后 CA 的 MHV
% =========================================================================
function SaveDTAEA_MHV_BeforeChange(CA, Problem, OldStageInfo, CurrentStage)
% SaveDTAEA_MHV_BeforeChange
%
% 在 DTAEA 检测到环境变化后、执行环境响应前，
% 保存旧环境最后 CA 的 MHV。

    global MHV_TEMP_HISTORY
    global DTAEA_RUN_IDX

    if isempty(DTAEA_RUN_IDX)
        RunIdx = -1;
    else
        RunIdx = DTAEA_RUN_IDX;
    end

    if isempty(MHV_TEMP_HISTORY)
        MHV_TEMP_HISTORY = [];
    end

    if isempty(CA)
        return;
    end

    % ---------------------------------------------------------
    % 1. 获取旧阶段编号
    % ---------------------------------------------------------
    StageToSave = [];

    try
        if isfield(OldStageInfo, 'index') && ~isempty(OldStageInfo.index)
            StageToSave = OldStageInfo.index;
        end
    catch
        StageToSave = [];
    end

    if isempty(StageToSave)
        try
            if isstruct(CurrentStage) && isfield(CurrentStage, 'objStage')
                StageToSave = CurrentStage.objStage;
            end
        catch
            StageToSave = [];
        end
    end

    if isempty(StageToSave)
        return;
    end

    % ---------------------------------------------------------
    % 2. 计算 TargetGen
    % ---------------------------------------------------------
    FirstStageGens = 300;

    if isprop(Problem, 'FirstStageGens') && ~isempty(Problem.FirstStageGens)
        FirstStageGens = Problem.FirstStageGens;
    end

    if ~isprop(Problem, 'taut') || isempty(Problem.taut)
        return;
    end

    if StageToSave == 1
        TargetGen = FirstStageGens - 1;
    else
        TargetGen = FirstStageGens + (StageToSave - 1) * Problem.taut - 1;
    end

    % 已保存则跳过
    if ~isempty(MHV_TEMP_HISTORY)
        if any(MHV_TEMP_HISTORY(:, 1) == TargetGen)
            return;
        end
    end

    % ---------------------------------------------------------
    % 3. 暂存真实 FE，并临时切到旧阶段末代
    % ---------------------------------------------------------
    oldFE = Problem.FE;
    tempFE = TargetGen * Problem.N + 1;

    try
        Problem.FE = tempFE;
    catch
    end

    % ---------------------------------------------------------
    % 4. 计算并保存 MHV
    % ---------------------------------------------------------
    try
        [score, HV_Pop, HV_PF] = MHV_Strict(CA, Problem);

        TriggerGen = floor(max(0, oldFE - 1) / Problem.N);

        fprintf(['   [Run %d] [DTAEA Internal Lock] Stage %d | ' ...
                 'TargetGen %d | TriggerGen %d | FE_before %d | ' ...
                 'HV_Pop: %.4f | HV_PF: %.4f | MHV: %.4f\n'], ...
                RunIdx, StageToSave, TargetGen, TriggerGen, oldFE, ...
                HV_Pop, HV_PF, score);

        % 保存格式：
        %   第1列 = TargetGen
        %   第2列 = TriggerGen
        %   第3列 = MHV
        MHV_TEMP_HISTORY = [MHV_TEMP_HISTORY; TargetGen, TriggerGen, score];

    catch ME
        fprintf('   [Run %d] [DTAEA Internal MHV Error] Stage %d | TargetGen %d | %s\n', ...
                RunIdx, StageToSave, TargetGen, ME.message);
    end

    % ---------------------------------------------------------
    % 5. 恢复真实 FE
    % ---------------------------------------------------------
    try
        Problem.FE = oldFE;
    catch
    end
end

% =========================================================================
% DTAEA 退出保护：从全局当前 CA 保存最后阶段 MHV
% =========================================================================
function SaveDTAEA_FinalStage_MHV_FromGlobal(Problem)

    global DTAEA_CURRENT_CA

    try
        if ~isempty(DTAEA_CURRENT_CA)
            SaveDTAEA_FinalStage_MHV(DTAEA_CURRENT_CA, Problem);
        end
    catch ME
        fprintf('[DTAEA Final Cleanup Error] %s\n', ME.message);
    end
end
% =========================================================================
% DTAEA 内部保存：最后一个阶段末代 CA 的 MHV
% =========================================================================
% =========================================================================
% DTAEA 内部保存：最后一个阶段末代 CA 的 MHV
% =========================================================================
% =========================================================================
% DTAEA 内部保存：最后一个阶段末代 CA 的 MHV
% =========================================================================
% =========================================================================
% DTAEA 内部保存：最后一个阶段末代 CA 的 MHV
% =========================================================================
function SaveDTAEA_FinalStage_MHV(CA, Problem)
% SaveDTAEA_FinalStage_MHV
%
% 算法结束时保存最后一个阶段的末代 CA 的 MHV。
%
% 关键修复：
%   1. 优先使用 RunExperiments_WithMHV 传入的 DTAEA_TOTAL_STAGES；
%   2. 不再只依赖 Problem.ObjGroups；
%   3. 根据 MHV_TEMP_HISTORY 已保存的行数，保存下一个缺失阶段；
%   4. 如果前面已经保存了 Stage 1~8，则这里保存 Stage 9。

    global MHV_TEMP_HISTORY
    global DTAEA_RUN_IDX
    global DTAEA_TOTAL_STAGES

    if isempty(DTAEA_RUN_IDX)
        RunIdx = -1;
    else
        RunIdx = DTAEA_RUN_IDX;
    end

    if isempty(MHV_TEMP_HISTORY)
        MHV_TEMP_HISTORY = [];
    end

    if isempty(CA)
        fprintf('   [Run %d] [DTAEA Final Skip] CA is empty.\n', RunIdx);
        return;
    end

    % ---------------------------------------------------------
    % 1. 读取总阶段数
    % ---------------------------------------------------------
    TotalStages = [];

    % 优先使用实验脚本传入的总阶段数
    if ~isempty(DTAEA_TOTAL_STAGES)
        TotalStages = DTAEA_TOTAL_STAGES;
    end

    % 兜底：从 Problem.ObjGroups 读取
    if isempty(TotalStages)
        try
            if isprop(Problem, 'ObjGroups') && ~isempty(Problem.ObjGroups)
                TotalStages = length(Problem.ObjGroups);
            end
        catch
            TotalStages = [];
        end
    end

    if isempty(TotalStages)
        fprintf('   [Run %d] [DTAEA Final Skip] Cannot determine TotalStages.\n', RunIdx);
        return;
    end

    % ---------------------------------------------------------
    % 2. 根据已保存数量，确定下一个缺失阶段
    % ---------------------------------------------------------
    SavedCount = size(MHV_TEMP_HISTORY, 1);
    StageToSave = SavedCount + 1;

    if StageToSave > TotalStages
        return;
    end

    if StageToSave < TotalStages
        fprintf(['   [Run %d] [DTAEA Final Warning] ' ...
                 'Only %d stages saved before final lock, next missing stage is %d, total stages is %d.\n'], ...
                RunIdx, SavedCount, StageToSave, TotalStages);
    end

    % ---------------------------------------------------------
    % 3. 读取 FirstStageGens 和 tau
    % ---------------------------------------------------------
    FirstStageGens = 300;

    if isprop(Problem, 'FirstStageGens') && ~isempty(Problem.FirstStageGens)
        FirstStageGens = Problem.FirstStageGens;
    end

    if ~isprop(Problem, 'taut') || isempty(Problem.taut)
        fprintf('   [Run %d] [DTAEA Final Skip] Problem.taut is missing.\n', RunIdx);
        return;
    end

    % ---------------------------------------------------------
    % 4. 计算该阶段理论末代 TargetGen
    % ---------------------------------------------------------
    if StageToSave == 1
        TargetGen = FirstStageGens - 1;
    else
        TargetGen = FirstStageGens + (StageToSave - 1) * Problem.taut - 1;
    end

    % 如果已经保存过这个 TargetGen，则不重复保存
    if ~isempty(MHV_TEMP_HISTORY)
        if any(MHV_TEMP_HISTORY(:, 1) == TargetGen)
            return;
        end
    end

    % ---------------------------------------------------------
    % 5. 暂存真实 FE，并临时切到该阶段末代
    % ---------------------------------------------------------
    oldFE = Problem.FE;
    tempFE = TargetGen * Problem.N + 1;

    try
        Problem.FE = tempFE;
    catch
    end

    % ---------------------------------------------------------
    % 6. 计算并保存 MHV
    % ---------------------------------------------------------
    try
        [score, HV_Pop, HV_PF] = MHV_Strict(CA, Problem);

        TriggerGen = floor(max(0, oldFE - 1) / Problem.N);

        fprintf(['   [Run %d] [DTAEA Final Lock] Stage %d | ' ...
                 'TargetGen %d | TriggerGen %d | FE_end %d | ' ...
                 'HV_Pop: %.4f | HV_PF: %.4f | MHV: %.4f\n'], ...
                RunIdx, StageToSave, TargetGen, TriggerGen, oldFE, ...
                HV_Pop, HV_PF, score);

        % 保存格式：
        %   第1列 = TargetGen
        %   第2列 = TriggerGen
        %   第3列 = MHV
        MHV_TEMP_HISTORY = [MHV_TEMP_HISTORY; TargetGen, TriggerGen, score];

    catch ME
        fprintf('   [Run %d] [DTAEA Final MHV Error] Stage %d | TargetGen %d | %s\n', ...
                RunIdx, StageToSave, TargetGen, ME.message);
    end

    % ---------------------------------------------------------
    % 7. 恢复真实 FE
    % ---------------------------------------------------------
    try
        Problem.FE = oldFE;
    catch
    end
end

% =========================================================================
% 垂直距离关联
% =========================================================================
function Region = GetAssociation(NormObjs, W)

    [N, ~] = size(NormObjs);
    [NW, ~] = size(W);

    Region = zeros(N, 1);

    for i = 1 : N

        x = NormObjs(i, :);
        minDist = inf;
        bestK = 1;

        for k = 1 : NW

            w = W(k, :);
            w_norm_sq = sum(w.^2);

            if w_norm_sq < 1e-10
                proj = zeros(size(w));
            else
                scalar_proj = (x * w') / w_norm_sq;
                proj = scalar_proj * w;
            end

            d_perp = norm(x - proj);

            if d_perp < minDist
                minDist = d_perp;
                bestK = k;
            end
        end

        Region(i) = bestK;
    end
end


% =========================================================================
% 密度估计
% =========================================================================
function Density = EstimateDensity(Objs, W)

    if isempty(Objs)
        Density = [];
        return;
    end

    Zmin = min(Objs, [], 1);
    Zmax = max(Objs, [], 1);

    NormObjs = (Objs - Zmin) ./ (Zmax - Zmin + 1e-6);

    Region = GetAssociation(NormObjs, W);

    N_W = size(W, 1);
    Counts = histcounts(Region, 1:N_W+1);

    Density = Counts(Region);
end


% =========================================================================
% LHS 采样
% =========================================================================
function Pop = LatinHypercubeSampling(Problem, N)

    if N <= 0
        Pop = SOLUTION.empty();
        return;
    end

    try
        X = lhsdesign(N, Problem.D);
    catch
        X = rand(N, Problem.D);
    end

    X = X .* (Problem.upper - Problem.lower) + Problem.lower;

    Pop = Problem.Evaluation(X);
end


% =========================================================================
% UpdateCA
% =========================================================================
function NewCA = UpdateCA(CA, Q, W, N)

    R = [CA, Q];

    [FrontNo, ~] = NDSort(R.objs, inf);

    NextCA = [];
    i = 1;

    while length(NextCA) + length(find(FrontNo == i)) <= N
        NextCA = [NextCA, R(FrontNo == i)];
        i = i + 1;
    end

    LastFront = R(FrontNo == i);
    Candidates = [NextCA, LastFront];

    if length(Candidates) <= N
        NewCA = Candidates;
        return;
    end

    Objs = Candidates.objs;

    Zmin = min(Objs, [], 1);
    Zmax = max(Objs, [], 1);

    NormObjs = (Objs - Zmin) ./ (Zmax - Zmin + 1e-6);

    Region = GetAssociation(NormObjs, W);

    CurrentSet = Candidates;
    CurrentObj = NormObjs;
    CurrentRegion = Region;

    while length(CurrentSet) > N

        [Counts, Edges] = histcounts(CurrentRegion, 1:size(W,1)+1);
        [~, CrowdedIdx] = max(Counts);
        CrowdedWIdx = Edges(CrowdedIdx);

        InWIdx = find(CurrentRegion == CrowdedWIdx);

        SubW = W(CrowdedWIdx, :);
        SubW(SubW < 1e-6) = 1e-6;

        TchDist = max(abs(CurrentObj(InWIdx, :)) ./ ...
            repmat(SubW, length(InWIdx), 1), [], 2);

        [~, WorstLocalIdx] = max(TchDist);
        RemoveIdx = InWIdx(WorstLocalIdx);

        CurrentSet(RemoveIdx) = [];
        CurrentObj(RemoveIdx, :) = [];
        CurrentRegion(RemoveIdx) = [];
    end

    NewCA = CurrentSet;
end


% =========================================================================
% UpdateDA
% =========================================================================
function NewDA = UpdateDA(CA, DA, Q, W, N)

    R = [DA, Q];

    AllObjs = [CA.objs; R.objs];

    Zmin = min(AllObjs, [], 1);
    Zmax = max(AllObjs, [], 1);

    NormR  = (R.objs  - Zmin) ./ (Zmax - Zmin + 1e-6);
    NormCA = (CA.objs - Zmin) ./ (Zmax - Zmin + 1e-6);

    RegionR  = GetAssociation(NormR, W);
    RegionCA = GetAssociation(NormCA, W);

    S = [];
    RemainingIdx = 1:length(R);

    itr = 1;

    while length(S) < N

        for i = 1 : size(W, 1)

            if length(S) >= N
                break;
            end

            CountCA = sum(RegionCA == i);

            if CountCA < itr

                CurrentCandidatesIdx = RemainingIdx(RegionR(RemainingIdx) == i);

                if ~isempty(CurrentCandidatesIdx)

                    CandObjs = NormR(CurrentCandidatesIdx, :);

                    [FrontNo, ~] = NDSort(CandObjs, 1);

                    BestCandidatesIdx = CurrentCandidatesIdx(FrontNo == 1);

                    SubW = W(i, :);
                    SubW(SubW < 1e-6) = 1e-6;

                    BestCandObjs = NormR(BestCandidatesIdx, :);

                    TchDist = max(abs(BestCandObjs) ./ ...
                        repmat(SubW, length(BestCandidatesIdx), 1), [], 2);

                    [~, MinLocalIdx] = min(TchDist);

                    TargetIdx = BestCandidatesIdx(MinLocalIdx);

                    S = [S, TargetIdx];

                    RemainingIdx(RemainingIdx == TargetIdx) = [];
                end
            end
        end

        itr = itr + 1;

        % 安全跳出，防止极端情况下死循环
        if itr > N + 10 && length(S) < N

            Rest = setdiff(1:length(R), S);
            Need = N - length(S);

            if length(Rest) >= Need
                S = [S, Rest(1:Need)];
            else
                S = [S, Rest];
            end

            break;
        end
    end

    NewDA = R(S);
end