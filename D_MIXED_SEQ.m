classdef D_MIXED_SEQ < PROBLEM
% <dynamic> <multi/many> <real> <large/none> <expensive/none>
% Dynamic Mixed-Benchmark (Standardized Version)
%
% Parameters:
%   taut        --- 25     --- Generations per stage
%   SeqString   --- '2,4;2,4,5;1,2,4,5;1,2,4,5,6;1,2,3,4,5,6;2,3,4,5,6;2,3,4,5;2,3,5;3,5' --- Active objectives
%   ObjConfig   --- '1-6:WFG9'  --- Problem types
%   FrozenStr   --- '0'             --- Frozen variables per stage

    properties(SetAccess = private)
        ObjGroups      % Active indices per stage
        ObjMap         % Map global index to Problem Name
        FrozenGroups   % Frozen variable indices per stage
        FrozenVal      % Fixed value for frozen variables
        taut           % Generations per stage
        lastT          % Last stage index
        GlobalMaxM     % Max dimension across all stages
        CurrentIndices % Current active indices
        WFG_k, WFG_l, WFG_S % WFG parameters
    end

    methods
        %% Default settings
        function Setting(obj)
            defaultSeqString = '2,4;2,4,5;1,2,4,5;1,2,4,5,6;1,2,3,4,5,6;2,3,4,5,6;2,3,4,5;2,3,5;3,5';
            defaultObjConfig = '1-6:WFG9'; 
            defaultFrozenStr = ''; 
            
            % 参数设置 (GUI 会读取上方的注释)
            [obj.taut, seqString, configString, frozenStr] = obj.ParameterSet(25, defaultSeqString, defaultObjConfig, defaultFrozenStr);
            
            obj.ObjGroups = obj.parseSequenceString(seqString);
            allIndices = [obj.ObjGroups{:}];
            obj.GlobalMaxM = max(allIndices);
            obj.ObjMap = obj.parseConfigString(configString, obj.GlobalMaxM);
            
            if isempty(obj.D), obj.D = obj.GlobalMaxM + 9; end
            
            obj.FrozenGroups = obj.parseSequenceString(frozenStr);
            obj.FrozenVal = 0.5; 
            
            % WFG Parameters
            obj.WFG_k = obj.GlobalMaxM - 1;
            obj.WFG_l = obj.D - obj.WFG_k;
            obj.WFG_S = 2 : 2 : 2*obj.GlobalMaxM; 
            
            obj.lower = zeros(1, obj.D);
            obj.upper = ones(1, obj.D);
            obj.encoding = ones(1, obj.D);
            
            obj.lastT = -1;
            obj.M = length(obj.ObjGroups{1});
        end

        %% Evaluation
        function Population = Evaluation(obj, varargin)
            CurrentGen = floor(obj.FE / obj.N);
            FirstStageGens = 300; % 保持300代延迟
            
            if CurrentGen < FirstStageGens
                t = 0;
            else
                t = 1 + floor((CurrentGen - FirstStageGens) / obj.taut);
            end

            if t > obj.lastT
                idx = min(t + 1, length(obj.ObjGroups));
                obj.CurrentIndices = obj.ObjGroups{idx};
                obj.M = length(obj.CurrentIndices); 
                obj.lastT = t;
            end

            PopDec = obj.CalDec(varargin{1});
            
            % Variable Freezing
            fIdx = min(t + 1, length(obj.FrozenGroups)); 
            currentFrozenVars = obj.FrozenGroups{fIdx};
            if ~isempty(currentFrozenVars)
                validVars = currentFrozenVars(currentFrozenVars <= size(PopDec, 2) & currentFrozenVars > 0);
                if ~isempty(validVars)
                    PopDec(:, validVars) = obj.FrozenVal;
                end
            end
            
            PopObj = zeros(size(PopDec,1), obj.M);
            activeIndices = obj.CurrentIndices;
            involvedProbs = unique(obj.ObjMap(activeIndices));
            
            for i = 1:length(involvedProbs)
                pName = involvedProbs{i};
                isThisProb = strcmp(obj.ObjMap(activeIndices), pName);
                targetCols = find(isThisProb);     
                globalCols = activeIndices(targetCols); 
                
                if isempty(globalCols), continue; end
                
                if startsWith(pName, 'DTLZ')
                    FullObj = obj.CalDTLZ(PopDec, pName);
                elseif startsWith(pName, 'WFG')
                    FullObj = obj.CalWFG(PopDec, pName);
                end
                PopObj(:, targetCols) = FullObj(:, globalCols);
            end

            % 保持最大化逻辑 (-1)
            PopObj = PopObj * (-1); 

            PopCon = obj.CalCon(PopDec);
            Population = SOLUTION(PopDec, PopObj, PopCon, zeros(size(PopDec,1), 1) + obj.FE);
            obj.FE = obj.FE + length(Population);
        end

        %% Calculation Engines (DTLZ Standard Checked)
        function PopObj = CalDTLZ(obj, PopDec, pName)
            M = obj.GlobalMaxM; [N, D] = size(PopDec);
            probIndex = str2double(pName(5:end));
            XM = PopDec(:, M:end);
            
            switch probIndex
                case 1
                    g = 100 * (D - M + 1 + sum((XM - 0.5).^2 - cos(20 * pi * (XM - 0.5)), 2));
                    PopObj = 0.5 * repmat(1 + g, 1, M) .* fliplr(cumprod([ones(N, 1), PopDec(:, 1:M-1)], 2)) .* [ones(N, 1), 1 - PopDec(:, M-1:-1:1)];
                case 2
                    g = sum((XM - 0.5).^2, 2);
                    PopObj = repmat(1 + g, 1, M) .* fliplr(cumprod([ones(N, 1), cos(PopDec(:, 1:M-1) * pi / 2)], 2)) .* [ones(N, 1), sin(PopDec(:, M-1:-1:1) * pi / 2)];
                case 3
                    g = 100 * (D - M + 1 + sum((XM - 0.5).^2 - cos(20 * pi * (XM - 0.5)), 2));
                    PopObj = repmat(1 + g, 1, M) .* fliplr(cumprod([ones(N, 1), cos(PopDec(:, 1:M-1) * pi / 2)], 2)) .* [ones(N, 1), sin(PopDec(:, M-1:-1:1) * pi / 2)];
                case 4
                    g = sum((XM - 0.5).^2, 2);
                    PopDec(:, 1:M-1) = PopDec(:, 1:M-1).^100;
                    PopObj = repmat(1 + g, 1, M) .* fliplr(cumprod([ones(N, 1), cos(PopDec(:, 1:M-1) * pi / 2)], 2)) .* [ones(N, 1), sin(PopDec(:, M-1:-1:1) * pi / 2)];
                case 5
                    g = sum((XM - 0.5).^2, 2);
                    Theta = zeros(N, M-1); Theta(:, 1) = PopDec(:, 1) * pi / 2;
                    for i = 2 : M-1, Theta(:, i) = (pi ./ (4 * (1 + g))) .* (1 + 2 * g .* PopDec(:, i)); end
                    PopObj = repmat(1 + g, 1, M) .* fliplr(cumprod([ones(N, 1), cos(Theta)], 2)) .* [ones(N, 1), sin(Theta(:, M-1:-1:1))];
                case 6
                    g = sum(XM.^0.1, 2);
                    Theta = zeros(N, M-1); Theta(:, 1) = PopDec(:, 1) * pi / 2;
                    for i = 2 : M-1, Theta(:, i) = (pi ./ (4 * (1 + g))) .* (1 + 2 * g .* PopDec(:, i)); end
                    PopObj = repmat(1 + g, 1, M) .* fliplr(cumprod([ones(N, 1), cos(Theta)], 2)) .* [ones(N, 1), sin(Theta(:, M-1:-1:1))];
                case 7
                    g = 1 + 9 * mean(XM, 2);
                    PopObj = zeros(N, M); PopObj(:, 1:M-1) = PopDec(:, 1:M-1);
                    h = M - sum(PopObj(:, 1:M-1) ./ (1 + repmat(g, 1, M-1)) .* (1 + sin(3 * pi * PopObj(:, 1:M-1))), 2);
                    PopObj(:, M) = (1 + g) .* h;
            end
        end

        %% Calculation Engines (WFG Standard Checked & Fixed)
        function PopObj = CalWFG(obj, PopDec, pName)
            M = obj.GlobalMaxM; [N, D] = size(PopDec);
            k = obj.WFG_k; l = obj.WFG_l;
            probIndex = str2double(pName(4:end)); z = PopDec; 
            
            switch probIndex
                case 1 % WFG1
                    % Step 1: 线性偏移 (仅作用于距离变量)
                    t1 = z; 
                    t1(:, k+1:end) = obj.wfg_s_linear(z(:, k+1:end), 0.35);
                    
                    % Step 2: 扁平化偏移 (位置变量) & 多项式偏移 (距离变量)
                    t2 = zeros(N, D);
                    t2(:, 1:k) = obj.wfg_b_flat(t1(:, 1:k), 0.8, 0.75, 0.85);
                    t2(:, k+1:end) = obj.wfg_b_poly(t1(:, k+1:end), 0.02);
                    
                    % Step 3: 加权求和降维 (注意 WFG1 的权重不是 ones，而是 2*idx)
                    t3 = zeros(N, M);
                    for i = 1:M-1
                        idx = (i-1)*k/(M-1)+1 : i*k/(M-1);
                        t3(:, i) = obj.wfg_r_sum(t2(:, idx), 2 * idx); 
                    end
                    idx_M = k+1:D;
                    t3(:, M) = obj.wfg_r_sum(t2(:, idx_M), 2 * idx_M);
                    
                    % Step 4: 形状映射 (Convex + Mixed)
                    x = t3; 
                    h = obj.shape_convex(x); 
                    h(:, M) = obj.shape_mixed(x(:, 1), 5, 1.0);
                
                case 2 % WFG2
                    t1 = z; t1(:, k+1:end) = obj.wfg_s_linear(z(:, k+1:end), 0.35);
                    t2 = zeros(N, k + l/2); t2(:, 1:k) = t1(:, 1:k); t2(:, k+1:end) = obj.wfg_r_nonsep(t1(:, k+1:end), 2);
                    t3 = zeros(N, M);
                    for i = 1:M-1, t3(:, i) = obj.wfg_r_sum(t2(:, (i-1)*k/(M-1)+1 : i*k/(M-1)), ones(1, k/(M-1))); end
                    t3(:, M) = obj.wfg_r_sum(t2(:, k+1:end), ones(1, size(t2,2)-k));
                    x = t3; h = obj.shape_convex(x); h(:, M) = obj.shape_disc(x(:, 1), 5, 1.0, 1.0);
                
                case 3 % WFG3
                    t1 = z; t1(:, k+1:end) = obj.wfg_s_linear(z(:, k+1:end), 0.35);
                    t2 = zeros(N, k + l/2); t2(:, 1:k) = t1(:, 1:k); t2(:, k+1:end) = obj.wfg_r_nonsep(t1(:, k+1:end), 2);
                    t3 = zeros(N, M);
                    for i = 1:M-1, t3(:, i) = obj.wfg_r_sum(t2(:, (i-1)*k/(M-1)+1 : i*k/(M-1)), ones(1, k/(M-1))); end
                    t3(:, M) = obj.wfg_r_sum(t2(:, k+1:end), ones(1, size(t2,2)-k));
                    x = t3; h = obj.shape_linear(x);
                
                case {4, 5, 6, 7, 8} % 【修复】WFG4-8 共享逻辑，剔除 WFG9
                    % Step 1: First Transformation
                    if probIndex == 4
                        t1 = obj.wfg_s_multi(z, 30, 10, 0.35);
                    elseif probIndex == 5
                        t1 = obj.wfg_s_decept(z, 0.35, 0.001, 0.05);
                    elseif probIndex == 6
                        t1 = z; t1(:, k+1:end) = obj.wfg_s_linear(z(:, k+1:end), 0.35);
                    elseif probIndex == 7
                        t1 = z; 
                        for i = 1:k
                            r = obj.wfg_r_sum(z(:, i+1:end), ones(1, D-i));
                            t1(:, i) = obj.wfg_b_param(z(:,i), r, 0.98/49.98, 0.02, 50);
                        end
                    elseif probIndex == 8
                        t1 = z;
                        r_pos = obj.wfg_r_sum(z(:, 1:k), ones(1, k));
                        for i = k+1:D
                             t1(:, i) = obj.wfg_b_param(z(:,i), r_pos, 0.98/49.98, 0.02, 50);
                        end
                    end
                    
                    % Step 2: Second Transformation
                    if probIndex == 6
                        reduced_dist = obj.wfg_r_nonsep(t1(:, k+1:end), l);
                        t2 = [t1(:, 1:k), reduced_dist];
                    elseif probIndex == 7 || probIndex == 8
                        t2 = t1; t2(:, k+1:end) = obj.wfg_s_linear(t1(:, k+1:end), 0.35);
                    else % WFG4, 5
                        t2 = t1; 
                    end
                    
                    % Step 3: Reduction (仅针对 WFG4-8)
                    t3 = zeros(N, M);
                    for i = 1:M-1
                        idx_start = (i-1)*k/(M-1)+1;
                        idx_end = i*k/(M-1);
                        t3(:, i) = obj.wfg_r_sum(t2(:, idx_start:idx_end), ones(1, length(idx_start:idx_end))); 
                    end
                    t3(:, M) = obj.wfg_r_sum(t2(:, k+1:end), ones(1, size(t2,2)-k));
                    x = t3; h = obj.shape_concave(x);

                case 9 % 【修复】完全独立的 WFG9 逻辑
                    t1 = zeros(N, D);
                    t2 = zeros(N, D);
                    
                    % Step 1: b_param
                    for i = 1 : D - 1
                        r = obj.wfg_r_sum(z(:, i+1:end), ones(1, D - i));
                        t1(:, i) = obj.wfg_b_param(z(:, i), r, 0.98/49.98, 0.02, 50);
                    end
                    t1(:, D) = z(:, D); 
                    
                    % Step 2: s_decept 与 s_multi
                    t2(:, 1:k) = obj.wfg_s_decept(t1(:, 1:k), 0.35, 0.001, 0.05);
                    t2(:, k+1:end) = obj.wfg_s_multi(t1(:, k+1:end), 30, 10, 0.35);
                    
                    % Step 3: r_nonsep
                    t3 = zeros(N, M);
                    for i = 1 : M - 1
                        idx_start = (i - 1) * k / (M - 1) + 1;
                        idx_end   = i * k / (M - 1);
                        t3(:, i) = obj.wfg_r_nonsep(t2(:, idx_start:idx_end), k / (M - 1));
                    end
                    t3(:, M) = obj.wfg_r_nonsep(t2(:, k+1:end), l);
                    
                    x = t3; 
                    h = obj.shape_concave(x);
            end
            
            PopObj = repmat(x(:, M), 1, M) + repmat(obj.WFG_S, N, 1) .* h;
        end
    

        %% Helpers
        function y = wfg_b_poly(obj, y, alpha), y = y.^alpha; end
        function y = wfg_b_flat(obj, y, A, B, C), y = A + min(0, floor(y - B)) .* A .* (B - y) ./ B - min(0, floor(C - y)) .* (1 - A) .* (y - C) ./ (1 - C); end
        function y = wfg_s_linear(obj, y, A), y = abs(y - A) ./ abs(floor(A - y) + A); end
        function y = wfg_s_decept(obj, y, A, B, C), y = 1 + (abs(y - A) - B) .* (floor(y - A + B) .* (1 - C + (A - B) / B) ./ (A - B) + floor(A + B - y) .* (1 - C + (1 - A - B) / B) ./ (1 - A - B) + 1 / B); end
        function y = wfg_s_multi(obj, y, A, B, C), term = abs(y - C) ./ (2 * (floor(C - y) + C)); y = (1 + cos((4 * A + 2) * pi * (0.5 - term)) + 4 * B * (term.^2)) ./ (B + 2); end
        function y = wfg_r_sum(obj, y, w), y = sum(y .* repmat(w, size(y, 1), 1), 2) ./ sum(w); end
        function y_out = wfg_r_nonsep(obj, y, A)
            [N, D] = size(y); n_blocks = ceil(D / A); y_out = zeros(N, n_blocks);
            for i = 1 : n_blocks
                sub_y = y(:, (i - 1) * A + 1 : min(i * A, D));
                y_out(:, i) = sum(sub_y, 2) + sum(abs(sub_y - repmat(1 + 2 * (0 : size(sub_y, 2) - 1), N, 1) .* repmat(sum(sub_y, 2) / size(sub_y, 2), 1, size(sub_y, 2))), 2);
            end
            y_out = y_out ./ ceil(A / 2) ./ (1 + 2 * A - 2 * ceil(A / 2));
        end
        % WFG7/8 Helper
        function y = wfg_b_param(obj, y, u, A, B, C)
             v = A - (1 - 2 * u) .* abs(floor(0.5 - u) + A);
             y = y.^(B + (C - B) * v);
        end
        
        function h = shape_linear(obj, x), M = size(x, 2); h = fliplr(cumprod([ones(size(x,1), 1), x(:, 1:M-1)], 2)) .* [ones(size(x,1), 1), 1 - x(:, M-1:-1:1)]; end
        function h = shape_convex(obj, x), M = size(x, 2); h = fliplr(cumprod([ones(size(x,1), 1), 1 - cos(x(:, 1:M-1) * pi / 2)], 2)) .* [ones(size(x,1), 1), 1 - sin(x(:, M-1:-1:1) * pi / 2)]; end
        function h = shape_concave(obj, x), M = size(x, 2); h = fliplr(cumprod([ones(size(x,1), 1), sin(x(:, 1:M-1) * pi / 2)], 2)) .* [ones(size(x,1), 1), cos(x(:, M-1:-1:1) * pi / 2)]; end
        function h = shape_mixed(obj, x1, A, alpha), h = (1 - x1 - cos(2 * A * pi * x1 + pi / 2) / (2 * A * pi)).^alpha; end
        function h = shape_disc(obj, x1, A, B, alpha), h = 1 - x1.^alpha .* cos(A * x1.^B * pi).^2; end
        
        function map = parseConfigString(obj, str, maxM)
            map = cell(1, maxM); map(:) = {'DTLZ2'}; segments = strsplit(strtrim(str), ';');
            for i = 1:length(segments)
                seg = strtrim(segments{i}); if isempty(seg), continue; end
                parts = strsplit(seg, ':'); rangeStr = strtrim(parts{1}); probName = upper(strtrim(parts{2}));
                indices = []; subParts = strsplit(rangeStr, ',');
                for k = 1:length(subParts)
                    sub = subParts{k};
                    if contains(sub, '-')
                        rangeBounds = str2double(strsplit(sub, '-')); indices = [indices, rangeBounds(1):rangeBounds(2)];
                    else, indices = [indices, str2double(sub)]; end
                end
                for idx = indices, if idx <= maxM, map{idx} = probName; end, end
            end
        end
        function groups = parseSequenceString(obj, str)
            if isempty(str), groups = {[]}; return; end
            str = strtrim(str); stages = strsplit(str, ';'); groups = cell(1, length(stages));
            for i = 1:length(stages)
                s = strtrim(stages{i});
                if isempty(s), groups{i} = []; else, numStrings = strsplit(s, ','); indices = zeros(1, length(numStrings)); for j = 1:length(numStrings), indices(j) = str2double(numStrings{j}); end; groups{i} = indices; end
            end
        end
        function DrawObj(obj, Population), PopObj = Population.objs; if obj.M == 2, Draw(PopObj, 'o', 'MarkerSize', 5, 'Markerfacecolor', [.5 .5 .5], 'Markeredgecolor', 'k', {'\it f\rm_1', '\it f\rm_2', []}); elseif obj.M == 3, Draw(PopObj, 'o', 'MarkerSize', 6, 'Markerfacecolor', [.5 .5 .5], 'Markeredgecolor', 'k', {'\it f\rm_1', '\it f\rm_2', '\it f\rm_3'}); else, Draw(PopObj); end, end
        function R = GetOptimum(obj, N)
            % 1. 动态计算当前所处的阶段 
            % 【时间差同步修复】
            CurrentGen = floor(max(0, obj.FE - 1) / obj.N);
            FirstStageGens = 300; 
            
            if CurrentGen < FirstStageGens
                t = 0;
            else
                t = 1 + floor((CurrentGen - FirstStageGens) / obj.taut);
            end
            
            % 2. 获取当前阶段真实的激活索引和维度 M
            idx = min(t + 1, length(obj.ObjGroups));
            realIndices = obj.ObjGroups{idx};
            currentM = length(realIndices);
            
            % 3. 确定当前的问题类型 (智能侦测)
            pName = obj.ObjMap{realIndices(1)};
            
            % 4. 计算真正的距离变量理论最大值 (max_xM)
            max_xM = 1.0; % WFG4,5,7,8 默认最大值为 1.0
            
            if startsWith(pName, 'WFG')
                probIndex = str2double(pName(4:end));
                
                % WFG2 和 WFG3 使用了 wfg_r_nonsep 且 A = 2
                if probIndex == 2 || probIndex == 3
                    max_xM = 4.0 / 3.0; % 理论极限 1.3333
                    
                % WFG6 和 WFG9 使用了 wfg_r_nonsep 且 A = obj.WFG_l
                elseif probIndex == 6 || probIndex == 9
                    l = obj.WFG_l;
                    c = ceil(l / 2);
                    max_xM = (l^2) / (c * (1 + 2 * l - 2 * c)); % 理论极限 1.81818...
                end
            end
            
            % 5. 生成单位均匀点 (凹面投影)
            % 注：WFG4~9 均为 Concave 形状，故使用球面投影
            R = UniformPoint(N, currentM);
            R = R ./ repmat(sqrt(sum(R.^2, 2)), 1, currentM);
            
            % 6. 获取当前激活目标的缩放因子 (S)
            currentS = obj.WFG_S(realIndices);
            
            % 7. 应用缩放因子: S_m * h_m
            R = repmat(currentS, size(R, 1), 1) .* R;
            
            % 8. 加上距离变量 x_M 的真实理论最大值
            R = R + max_xM;
            
            % 9. 转换为 Minus 形式 (映射到负值区间)
            R = R * (-1);
        end
    end
end