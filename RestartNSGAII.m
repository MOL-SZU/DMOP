classdef RestartNSGAII < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation> <dynamic>
% Nondominated sorting genetic algorithm II
% Baseline 版本：当环境变化时，完全丢弃旧种群，重新初始化
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% 1. 初始化种群
            Population = Problem.Initialization();
            [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Problem.N);
            

            %% 2. 优化循环
            while Algorithm.NotTerminated(Population) 
                
                % 生成子代
                MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                
                % 检测环境是否变化
                EnvironmentChanged = size(Offspring.objs,2) ~= size(Population.objs,2);
                
                if EnvironmentChanged
                   % 环境发生变化，直接重新初始化
                   Population = Problem.Initialization();
                   
                   % 重新计算新种群的非支配排序和拥挤距离
                   [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Problem.N);
                else
                   % 环境未变化，执行标准的 NSGA-II 环境选择
                   [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Problem.N);
                end
                % ====================================================
            end
            

        end
    end
end