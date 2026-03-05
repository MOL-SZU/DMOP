classdef ArchiveSPEA2SDE < SPEA2SDE
% <2025> <multi/many> <real/integer/label/binary/permutation> <dynamic>
% Archive-based SPEA2 with shift-based density estimation

    methods
        function main(Algorithm,Problem)
            %% 1. 初始化
            Population = Problem.Initialization();
            Fitness    = CalFitness(Population.objs);
            
            % 初始化外部存档
            ExternalArchive = [];

            %% 2. 主循环
            while Algorithm.NotTerminated(Population)
                
                MatingPool = TournamentSelection(2,Problem.N,Fitness);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                
                % 检测环境是否发生变化 (通过维度对比)
                EnvironmentChanged = size(Offspring.objs, 2) ~= size(Population.objs, 2);
                
                if EnvironmentChanged
                    % ==========================================================
                    %                环境变化处理 (基于存档的热启动)
                    % ==========================================================
                    
                    ArchiveDecs = [];
                    if ~isempty(ExternalArchive)
                        % 1. 存档去重
                        [~, DistinctIdx] = unique(ExternalArchive.objs, 'rows');
                        UniqueArchive    = ExternalArchive(DistinctIdx);
                        
                        % 2. 在旧环境下进行非支配排序，选 Rank 1
                        [FrontNo, ~]  = NDSort(UniqueArchive.objs, 1);
                        BestSolutions = UniqueArchive(FrontNo == 1);
                        
                        ArchiveDecs = BestSolutions.decs;
                    end
                    
                    % 3. 生成新种群决策变量
                    NumInArchive = size(ArchiveDecs, 1);
                    NewDecs = [];
                    
                    if NumInArchive >= Problem.N
                        % 存档优良解足够多，随机选 N 个
                        SelectIdx = randperm(NumInArchive, Problem.N);
                        NewDecs   = ArchiveDecs(SelectIdx, :);
                    elseif NumInArchive > 0
                        % 存档不够，随机补足
                        NumToRandom = Problem.N - NumInArchive;
                        RandomPop   = Problem.Initialization(NumToRandom);
                        NewDecs     = [ArchiveDecs; RandomPop.decs];
                    else
                        % 存档为空，全随机
                        RandomPop   = Problem.Initialization(Problem.N);
                        NewDecs     = RandomPop.decs;
                    end
                    
                    % 4. 在新环境下评估
                    Population = Problem.Evaluation(NewDecs);
                    Fitness    = CalFitness(Population.objs);
                    
                    % 5. 清空存档
                    ExternalArchive = [];
                    
                else
                    % ==========================================================
                    %                稳态进化
                    % ==========================================================
                    
                    [Population,Fitness] = EnvironmentalSelection([Population,Offspring],Problem.N);
                    
                    % 将子代累积到存档中
                    ExternalArchive = [ExternalArchive, Offspring];
                end
            end
        end
    end

end
