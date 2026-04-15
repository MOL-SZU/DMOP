function RunExperiments_WithMHV
% RunExperiments_WithMHV
% 记录"每个阶段最后一代"的 MHV 比值
%
% 输出:
%   MHV_Data{r} 为一个 [Stages x 2] 矩阵:
%   第1列 = 锁定代数
%   第2列 = 该代数对应的 MHV 比值

    % =========================================================================
    % 1. 基础参数
    % =========================================================================
    N = 300;
    Runs = 31;
    Taus = [25, 50, 100];

    AlgoList  = {@DTAEA, @KTDMOEA, @LEC, @STA, @ArchiveSPEA2SDE};
    AlgoNames = {'DTAEA', 'KTDMOEA', 'LEC', 'STA', 'ArchiveSPEA2SDE'};

    % =========================================================================
    % 2. 实验设置
    % =========================================================================
    Settings = struct();

    % --- Setting I ---
    Settings(1).Name    = 'Setting I (Mild)';
    Settings(1).ProbTag = 'DTLZ2';
    Settings(1).Config  = '1-6:DTLZ2';
    Settings(1).Seq     = '2,4; 2,4,5; 1,2,4,5; 1,2,4,5,6; 1,2,3,4,5,6; 2,3,4,5,6; 2,3,4,5; 2,3,5; 3,5';
    Settings(1).Frozen  = '';

    % --- Setting II ---
    Settings(2).Name    = 'Setting II (Moderate)';
    Settings(2).ProbTag = 'DTLZ2';
    Settings(2).Config  = '1-10:DTLZ2';
    Settings(2).Seq     = ['2,7; 2,5,7,10; 1,2,5,6,7,10; 1,2,4,5,6,7,9,10; ' ...
                           '1,2,3,4,5,6,7,8,9,10; 1,2,3,5,6,8,9,10; 2,3,5,6,9,10; 2,5,6,9; 5,6'];
    Settings(2).Frozen  = '';

    % --- Setting III ---
    Settings(3).Name    = 'Setting III (Severe)';
    Settings(3).ProbTag = 'DTLZ2';
    Settings(3).Config  = '1-10:DTLZ2';
    Settings(3).Seq     = ['3,8; 2,3,6,7,8; 1,2,3,4,5,6,7,8,9,10; ' ...
                           '1,3,5,6,7,10; 3,7,8; 1,3,4,5,6,7,8,9; 2,5,7,10; ' ...
                           '1,2,4,5,6,9,10; 1,2,3,4,5,6,7,8,10'];
    Settings(3).Frozen  = '';

    if ~exist('Results_StepMHV', 'dir')
        mkdir('Results_StepMHV');
    end

    global MHV_TEMP_HISTORY

    fprintf('=== 开始执行"每阶段最后一代"MHV实验 ===\n');

    for s = 1:length(Settings)
        CurrentSetting = Settings(s);

        for tau = Taus
            % 总阶段数
            Stages = length(strsplit(CurrentSetting.Seq, ';'));

            % 最后一个阶段最后一代 = 300 + (Stages-1)*tau - 1
            % 所以运行到 300 + (Stages-1)*tau 代即可
            MaxGen = 300 + (Stages - 1) * tau;
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
                    MHV_TEMP_HISTORY = [];

                    CurrentProblem = {@D_MIXED_SEQ, tau, CurrentSetting.Seq, ...
                                      CurrentSetting.Config, CurrentSetting.Frozen};

                    try
                        platemo('algorithm', AlgorithmFunc, ...
                                'problem', CurrentProblem, ...
                                'N', N, ...
                                'maxFE', MaxFE, ...
                                'save', 0, ...
                                'outputFcn', @(alg, pro) ...
                                    RecordMHV_Adaptive_Callback_Final(alg, pro, tau, r, Stages));

                        MHV_Data{r} = MHV_TEMP_HISTORY;

                    catch ME
                        fprintf('   [Run %d Error]: %s\n', r, ME.message);
                        MHV_Data{r} = [];
                    end
                end

                save(resFileName, 'MHV_Data', 'CurrentSetting', 'tau', 'AlgoName');
                fprintf('=> %s Done.\n', AlgoName);
            end
        end
    end

    clear global MHV_TEMP_HISTORY;
    fprintf('\n所有实验结束。\n');
end

% =========================================================================
% 记录"每个阶段最后一代"MHV 的回调
% =========================================================================
function RecordMHV_Adaptive_Callback_Final(Algorithm, Problem, tau, runIdx, TotalStages)
    global MHV_TEMP_HISTORY

    if isempty(MHV_TEMP_HISTORY)
        MHV_TEMP_HISTORY = [];
    end

    CurrentGen = ceil(Problem.FE / Problem.N);
    SavedCount = size(MHV_TEMP_HISTORY, 1);

    % 已经记录完所有阶段
    if SavedCount >= TotalStages
        return;
    end

    % 每个阶段最后一代：
    % Stage 1 -> 299
    % Stage 2 -> 300 + tau - 1
    % Stage 3 -> 300 + 2*tau - 1
    % ...
    if SavedCount == 0
        TargetGen = 299;
    else
        TargetGen = 300 + SavedCount * tau - 1;
    end

    if CurrentGen == TargetGen
        % 防止同一代重复记录
        if SavedCount == 0 || MHV_TEMP_HISTORY(end,1) ~= TargetGen
            Population = Algorithm.result;
            if ~isempty(Population)
                try
                    [score, HV_Pop, HV_PF] = MHV_Strict(Population, Problem);

                    fprintf('   [Run %d] Stage %d Locked. Gen %d | HV_Pop: %.4f | HV_PF: %.4f | MHV: %.4f\n', ...
                            runIdx, SavedCount + 1, TargetGen, HV_Pop, HV_PF, score);

                    % 记录格式：[锁定代数, MHV]
                    MHV_TEMP_HISTORY = [MHV_TEMP_HISTORY; TargetGen, score];
                catch ME
                    fprintf('   [Run %d] Error at Gen %d: %s\n', ...
                            runIdx, TargetGen, ME.message);
                end
            end
        end
    end
end