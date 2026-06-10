function RunExperiments_WithMHV
% RunExperiments_WithMHV
%
% 记录"每个环境变化前最后一个种群"的 MHV 比值。
%
% 输出:
%   MHV_Data{r} 为一个 [Stages x 3] 矩阵:
%
%   第1列 = TargetGen
%           理论锁定代数，即每个环境变化前最后一代
%
%   第2列 = TriggerGen
%           实际触发记录时，由 floor((FE-1)/N) 得到的代数
%
%   第3列 = MHV
%           该阶段变化前种群的 MHV 比值
%
% 注意:
%   DTAEA 使用算法内部保存 MHV。
%   其他算法仍然使用 outputFcn 保存。
%
% 关键修复:
%   1. 新增全局变量 DTAEA_TOTAL_STAGES；
%   2. 每次 run 前把当前 Stages 传给 DTAEA；
%   3. DTAEA 内部最终保存 Stage 9 时不再只依赖 Problem.ObjGroups。

    % =========================================================================
    % 1. 基础参数
    % =========================================================================
    N = 300;
    Runs = 31;
    Taus = [25, 50, 100];

    AlgoList  = {@DTAEA, @KTDMOEA, @LEC, @STA};
    AlgoNames = {'DTAEA', 'KTDMOEA', 'LEC', 'STA'};
    % AlgoList  = {@ArchiveSPEA2SDE};
    % AlgoNames = {'ArchiveSPEA2SDE'};

    % =========================================================================
    % 2. 实验设置
    % =========================================================================
    Settings = struct();

    % --- Setting I ---
    Settings(1).Name    = 'Setting I (Mild)';
    Settings(1).ProbTag = 'DTLZ4';
    Settings(1).Config  = '1-6:DTLZ4';
    Settings(1).Seq     = ['2,4; ' ...
                            '2,4,5; ' ...
                            '1,2,4,5; ' ...
                            '1,2,4,5,6; ' ...
                            '1,2,3,4,5,6; ' ...
                            '2,3,4,5,6; ' ...
                            '2,3,4,5; ' ...
                            '2,3,5; ' ...
                            '3,5'];
    Settings(1).Frozen  = '';

    % --- Setting II ---
    Settings(2).Name    = 'Setting II (Moderate)';
    Settings(2).ProbTag = 'DTLZ4';
    Settings(2).Config  = '1-10:DTLZ4';
    Settings(2).Seq     = ['2,7; ' ...
                            '2,5,7,10; ' ...
                            '1,2,5,6,7,10; ' ...
                            '1,2,4,5,6,7,9,10; ' ...
                            '1,2,3,4,5,6,7,8,9,10; ' ...
                            '1,2,3,5,6,8,9,10; ' ...
                            '2,3,5,6,9,10; ' ...
                            '2,5,6,9; ' ...
                            '5,6'];
    Settings(2).Frozen  = '';

    % --- Setting III ---
    Settings(3).Name    = 'Setting III (Severe)';
    Settings(3).ProbTag = 'DTLZ4';
    Settings(3).Config  = '1-10:DTLZ4';
    Settings(3).Seq     = ['3,8; ' ...
                            '2,3,6,7,8; ' ...
                            '1,2,3,4,5,6,7,8,9,10; ' ...
                            '1,3,5,6,7,10; ' ...
                            '3,7,8; ' ...
                            '1,3,4,5,6,7,8,9; ' ...
                            '2,5,7,10; ' ...
                            '1,2,4,5,6,9,10; ' ...
                            '1,2,3,4,5,6,7,8,10'];
    Settings(3).Frozen  = '';

    % =========================================================================
    % 3. 输出目录
    % =========================================================================
    if ~exist('Results_StepMHV', 'dir')
        mkdir('Results_StepMHV');
    end

    global MHV_TEMP_HISTORY
    global DTAEA_CURRENT_CA
    global DTAEA_RUN_IDX
    global DTAEA_TOTAL_STAGES

    fprintf('=== 开始执行"环境变化前最后种群"MHV实验 ===\n');

    % =========================================================================
    % 4. 主循环
    % =========================================================================
    for s = 1:length(Settings)

        CurrentSetting = Settings(s);

        for tau = Taus

            % 总阶段数
            Stages = length(strsplit(CurrentSetting.Seq, ';'));

            % 默认第一阶段长度。
            % 如果 D_MIXED_SEQ 内部定义 FirstStageGens=300，这里保持一致。
            FirstStageGens_Default = 300;

            % 最后一个阶段末代 TargetGen:
            %   FirstStageGens + (Stages-1)*tau - 1
            %
            % 为了让算法运行跨过最后一个 TargetGen，运行到:
            %   FirstStageGens + (Stages-1)*tau
            MaxGen = FirstStageGens_Default + (Stages - 1) * tau;
            MaxFE  = MaxGen * N;

            fprintf('\n=== %s [%s], Tau: %d (Stages: %d) ===\n', ...
                    CurrentSetting.Name, CurrentSetting.ProbTag, tau, Stages);

            for a = 1:length(AlgoList)

                AlgorithmFunc = AlgoList{a};
                AlgoName = AlgoNames{a};

                resFileName = fullfile('Results_StepMHV', ...
                    sprintf('MHV_Step_%s_%s_Set%d_Tau%d.mat', ...
                    AlgoName, CurrentSetting.ProbTag, s, tau));

                if exist(resFileName, 'file')
                    fprintf('Skipping %s (Exists)\n', AlgoName);
                    continue;
                end

                fprintf('Running %s ...\n', AlgoName);

                MHV_Data = cell(1, Runs);

                for r = 1:Runs

                    % 每次运行前清空临时历史
                    MHV_TEMP_HISTORY = [];

                    % 清空 DTAEA 当前 CA，防止上一轮残留
                    DTAEA_CURRENT_CA = [];

                    % 把当前 run 编号和总阶段数传给 DTAEA 内部保存函数
                    DTAEA_RUN_IDX = r;
                    DTAEA_TOTAL_STAGES = Stages;

                    CurrentProblem = {@D_MIXED_SEQ, tau, CurrentSetting.Seq, ...
                                      CurrentSetting.Config, CurrentSetting.Frozen};

                    try
                        platemo('algorithm', AlgorithmFunc, ...
                                'problem', CurrentProblem, ...
                                'N', N, ...
                                'maxFE', MaxFE, ...
                                'save', 0, ...
                                'outputFcn', @(alg, pro) ...
                                    RecordMHV_Adaptive_Callback_Final( ...
                                        alg, pro, tau, r, Stages, AlgoName));

                        MHV_Data{r} = MHV_TEMP_HISTORY;

                    catch ME
                        fprintf('   [Run %d Error]: %s\n', r, ME.message);
                        MHV_Data{r} = [];
                    end
                end

                save(resFileName, ...
                     'MHV_Data', ...
                     'CurrentSetting', ...
                     'tau', ...
                     'AlgoName');

                fprintf('=> %s Done.\n', AlgoName);
            end
        end
    end

    clear global MHV_TEMP_HISTORY
    clear global DTAEA_CURRENT_CA
    clear global DTAEA_RUN_IDX
    clear global DTAEA_TOTAL_STAGES

    fprintf('\n所有实验结束。\n');
end


% =========================================================================
% 记录"每个环境变化前最后种群"MHV 的回调
% =========================================================================
function RecordMHV_Adaptive_Callback_Final(Algorithm, Problem, tau, runIdx, TotalStages, AlgoName)
% RecordMHV_Adaptive_Callback_Final
%
% 非 DTAEA 算法使用该回调保存 MHV。
%
% DTAEA 已经在 DTAEA.m 内部保存：
%   Stage 1 ~ Stage 倒数第二阶段：环境变化前保存
%   最后阶段：算法结束后保存
%
% 因此 DTAEA 在这里直接 return，避免重复保存或漏记。

    global MHV_TEMP_HISTORY

    % DTAEA 使用算法内部保存
    if strcmp(AlgoName, 'DTAEA')
        return;
    end

    if isempty(MHV_TEMP_HISTORY)
        MHV_TEMP_HISTORY = [];
    end

    % 已经保存了多少个阶段
    SavedCount = size(MHV_TEMP_HISTORY, 1);

    % 已经记录完所有阶段
    if SavedCount >= TotalStages
        return;
    end

    % 当前需要保存的阶段编号
    StageToSave = SavedCount + 1;

    % 当前种群对应的已完成代数
    CurrentGen = floor(max(0, Problem.FE - 1) / Problem.N);

    % 从 Problem 读取 FirstStageGens
    FirstStageGens = GetFirstStageGens(Problem);

    % 当前阶段的理论锁定代数
    if StageToSave == 1
        TargetGen = FirstStageGens - 1;
    else
        TargetGen = FirstStageGens + (StageToSave - 1) * tau - 1;
    end

    % 阶段保护
    CurrentStage = GetStageByFEminus1(Problem);

    if CurrentStage ~= StageToSave
        return;
    end

    % 使用 >=，防止跳过目标代数时漏记
    if CurrentGen < TargetGen
        return;
    end

    % 防止同一阶段重复记录
    if SavedCount > 0 && MHV_TEMP_HISTORY(end, 1) == TargetGen
        return;
    end

    % 取种群
    Population = [];

    try
        Population = Algorithm.result;
    catch
        Population = [];
    end

    if isempty(Population)
        return;
    end

    % 计算 MHV
    try
        [score, HV_Pop, HV_PF] = MHV_Strict(Population, Problem);

        fprintf(['   [Run %d] Stage %d Locked. ' ...
                 'TargetGen %d | TriggerGen %d | ' ...
                 'HV_Pop: %.4f | HV_PF: %.4f | MHV: %.4f\n'], ...
                runIdx, StageToSave, TargetGen, CurrentGen, ...
                HV_Pop, HV_PF, score);

        % 保存格式:
        %   第1列 TargetGen
        %   第2列 TriggerGen
        %   第3列 MHV
        MHV_TEMP_HISTORY = [MHV_TEMP_HISTORY; TargetGen, CurrentGen, score];

    catch ME
        fprintf('   [Run %d] Error at Stage %d TargetGen %d TriggerGen %d: %s\n', ...
                runIdx, StageToSave, TargetGen, CurrentGen, ME.message);
    end
end


% =========================================================================
% 用 FE-1 判断当前种群所属阶段
% =========================================================================
function Stage = GetStageByFEminus1(Problem)

    CurrentGen = floor(max(0, Problem.FE - 1) / Problem.N);

    FirstStageGens = GetFirstStageGens(Problem);

    if isprop(Problem, 'taut') && ~isempty(Problem.taut)
        tau = Problem.taut;
    else
        tau = 50;
    end

    if CurrentGen < FirstStageGens
        t = 0;
    else
        t = 1 + floor((CurrentGen - FirstStageGens) / tau);
    end

    if isprop(Problem, 'ObjGroups') && ~isempty(Problem.ObjGroups)
        Stage = min(t + 1, length(Problem.ObjGroups));
    else
        Stage = t + 1;
    end
end


% =========================================================================
% 读取 FirstStageGens
% =========================================================================
function FirstStageGens = GetFirstStageGens(Problem)

    FirstStageGens = 300;

    if isprop(Problem, 'FirstStageGens') && ~isempty(Problem.FirstStageGens)
        FirstStageGens = Problem.FirstStageGens;
    end
end