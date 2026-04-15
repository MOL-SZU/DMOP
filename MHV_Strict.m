function [score, HV_Pop, HV_PF] = MHV_Strict(Population, Problem)
% MHV_Strict


    persistent Last_PF_Signature Last_HV_PF

    score  = 0;
    HV_Pop = 0;
    HV_PF  = 0;

    %% 1. 提取种群原始目标值
    try
        PopObj = ExtractPop(Population);
        if isempty(PopObj)
            return;
        end
        PopObj = double(PopObj);
    catch ME
        fprintf('[MHV Error 1] %s\n', ME.message);
        return;
    end

    %% 2. 提取参考 PF
    try
        if ~ismethod(Problem, 'GetOptimum')
            return;
        end

        RawPF = Problem.GetOptimum(5000);
        PF = ExtractPF(RawPF);
        PF = double(PF);

        % fprintf('[DEBUG] after GetOptimum: %d x %d\n', size(PF,1), size(PF,2));

        if isempty(PF)
            return;
        end

        [N_pop, M_pop] = size(PopObj);
        [~, M_pf] = size(PF);

        % 维度对齐
        if M_pf > M_pop
            PF = PF(:, 1:M_pop);
        elseif M_pf < M_pop
            return;
        end

        % 清理非法值
        PF(any(isnan(PF), 2) | any(isinf(PF), 2), :) = [];
        if isempty(PF)
            return;
        end

        % 去重 + 第一前沿
        PF = unique(PF, 'rows');
        if exist('NDSort', 'file')
            [FrontNo, ~] = NDSort(PF, 1);
        else
            [FrontNo, ~] = NDSort_Local_Fast(PF, 1);
        end
        PF = PF(FrontNo == 1, :);

        if isempty(PF)
            return;
        end

        % fprintf('[DEBUG] after unique + first front: %d x %d\n', size(PF,1), size(PF,2));
        % 如果点太多，降采样时强制保留边界点
        % if size(PF,1) > 5000
        %     PF = DownsamplePF_KeepBoundary(PF, 5000);
        % 
        %     if exist('NDSort', 'file')
        %         [FrontNo, ~] = NDSort(PF, 1);
        %     else
        %         [FrontNo, ~] = NDSort_Local_Fast(PF, 1);
        %     end
        %     PF = PF(FrontNo == 1, :);
        % end

    catch ME
        fprintf('[MHV Error 2] %s\n', ME.message);
        return;
    end

    [N_pop, M] = size(PopObj);
    [N_pf, ~]  = size(PF);

    %% 3. 归一化
    try
        z_star = min(PF, [], 1);
        z_nad  = max(PF, [], 1);

        gap = z_nad - z_star;
        gap(gap < 1e-12) = 1e-12;

        NormPop = (PopObj - repmat(z_star, N_pop, 1)) ./ repmat(gap, N_pop, 1);
        NormPF  = (PF     - repmat(z_star, N_pf,  1)) ./ repmat(gap, N_pf,  1);

        % 删除非法点
        validPop = all(isfinite(NormPop), 2);
        NormPop = NormPop(validPop, :);

        if isempty(NormPop)
            score = 0;
            return;
        end

        validPF = all(isfinite(NormPF), 2);
        NormPF = NormPF(validPF, :);

        if isempty(NormPF)
            score = 0;
            return;
        end

        % 边缘误差修正
        NormPop = max(NormPop, 0);
        NormPF  = max(NormPF, 0);

    catch ME
        fprintf('[MHV Error 3] %s\n', ME.message);
        return;
    end

    %% 4. HV 计算
    try
        UpperBound = 1.1;
        RefPoint = ones(1, M) * UpperBound;

        % 只保留进入参考盒子的点
        FilteredPop = NormPop(all(NormPop <= UpperBound, 2), :);
        FilteredPF  = NormPF(all(NormPF <= UpperBound, 2), :);
        
        % fprintf('[DEBUG] FilteredPop: %d x %d\n', size(FilteredPop,1), size(FilteredPop,2));
        % fprintf('[DEBUG] FilteredPF: %d x %d\n', size(FilteredPF,1), size(FilteredPF,2));


        % ---------- Pop HV ----------
        if isempty(FilteredPop)
            HV_Pop = 0;
        else
            HV_Pop = CalHV_Hybrid(FilteredPop, RefPoint, true);
        end

        % ---------- PF HV（带缓存） ----------
        if isempty(FilteredPF)
            HV_PF = 0;
            Last_PF_Signature = [];
        else
            Current_Signature = GetPFSignature(Problem, M, size(FilteredPF,1));

            if ~isempty(Last_PF_Signature) && isequal(Current_Signature, Last_PF_Signature)
                HV_PF = Last_HV_PF;
            else
                HV_PF = CalHV_Hybrid(FilteredPF, RefPoint, true);
                Last_HV_PF = HV_PF;
                Last_PF_Signature = Current_Signature;
            end
        end

        % ---------- 原始比值 ----------
        if HV_PF < 1e-12
            score = 0;
        else
            score = HV_Pop / HV_PF;
        end

    catch ME
        fprintf('[MHV Error 4] %s\n', ME.message);
        return;
    end
end

% =========================================================================
% PF 签名：按当前阶段缓存 PF_HV
% =========================================================================
function Sig = GetPFSignature(Problem, M, PFCount)
    Sig = sprintf('M=%d|PF=%d', M, PFCount);
end

% =========================================================================
% 核心 HV 计算
% =========================================================================
function Volume = CalHV_Hybrid(Points, RefPoint, isPF)
    if nargin < 3
        isPF = false;
    end

    [~, M] = size(Points);

    % 只保留第一前沿
    if exist('NDSort', 'file')
        [FrontNo, ~] = NDSort(Points, 1);
    else
        [FrontNo, ~] = NDSort_Local_Fast(Points, inf);
    end
    Points = Points(FrontNo == 1, :);
    Points = unique(Points, 'rows');

    if isempty(Points)
        Volume = 0;
        return;
    end

    N = size(Points,1);

    % ---------- 2D 精确 ----------
    if M == 2
        [~, sortIdx] = sort(Points(:,1));
        SortedP = Points(sortIdx, :);

        Volume = (RefPoint(1) - SortedP(end,1)) * (RefPoint(2) - SortedP(end,2));
        for i = size(SortedP,1)-1 : -1 : 1
            Volume = Volume + (SortedP(i+1,1) - SortedP(i,1)) * (RefPoint(2) - SortedP(i,2));
        end
        Volume = max(Volume, 0);
        return;
    end

    % ---------- 3D 小规模精确 ----------
    if M == 3 && N < 500
        Volume = CalHV_Exact_3D(Points, RefPoint, M);
        return;
    end

    % ---------- 高维 Monte Carlo ----------
    if isPF
        SampleNum = 1000000;
    else
        SampleNum = 50000;
    end

    MinValue = min(Points, [], 1);
    MaxValue = RefPoint;

    Samples = unifrnd(repmat(MinValue, SampleNum, 1), repmat(MaxValue, SampleNum, 1));

    for i = 1:size(Points,1)
        if isempty(Samples)
            break;
        end
        domi = true(size(Samples,1),1);
        m = 1;
        while m <= M && any(domi)
            domi = domi & Points(i,m) <= Samples(:,m);
            m = m + 1;
        end
        Samples(domi,:) = [];
    end

    Volume = prod(MaxValue - MinValue) * (1 - size(Samples,1) / SampleNum);
end

% =========================================================================
% 3D 精确 HV
% =========================================================================
function Volume = CalHV_Exact_3D(Points, RefPoint, M)
    Volume = 0;
    pl = sortrows(Points);
    S = {1, pl};

    for k = 1 : M-1
        S_ = {};
        for i = 1 : size(S,1)
            Stemp = Slice(pluck(S,i,2), k, RefPoint);
            for j = 1 : size(Stemp,1)
                temp(1) = {cell2mat(Stemp(j,1)) * cell2mat(S(i,1))};
                temp(2) = Stemp(j,2);
                S_ = Add(temp, S_);
            end
        end
        S = S_;
    end

    for i = 1 : size(S,1)
        p = Head(pluck(S,i,2));
        Volume = Volume + cell2mat(S(i,1)) * abs(p(M) - RefPoint(M));
    end
end

function out = pluck(C, i, j)
    out = cell2mat(C(i,j));
end

function S = Slice(pl, k, RefPoint)
    p = Head(pl);
    pl = Tail(pl);
    ql = [];
    S = {};

    while ~isempty(pl)
        ql = Insert(p, k+1, ql);
        p_ = Head(pl);
        cell_(1,1) = {abs(p(k) - p_(k))};
        cell_(1,2) = {ql};
        S = Add(cell_, S);
        p = p_;
        pl = Tail(pl);
    end

    ql = Insert(p, k+1, ql);
    cell_(1,1) = {abs(p(k) - RefPoint(k))};
    cell_(1,2) = {ql};
    S = Add(cell_, S);
end

function ql = Insert(p, k, pl)
    ql = [];
    hp = Head(pl);

    while ~isempty(pl) && hp(k) < p(k)
        ql = [ql; hp];
        pl = Tail(pl);
        hp = Head(pl);
    end

    ql = [ql; p];
    m = length(p);

    while ~isempty(pl)
        q = Head(pl);
        flag1 = 0;
        flag2 = 0;
        for i = k : m
            if p(i) < q(i)
                flag1 = 1;
            elseif p(i) > q(i)
                flag2 = 1;
            end
        end
        if ~(flag1 == 1 && flag2 == 0)
            ql = [ql; Head(pl)];
        end
        pl = Tail(pl);
    end
end

function p = Head(pl)
    if isempty(pl)
        p = [];
    else
        p = pl(1,:);
    end
end

function ql = Tail(pl)
    if size(pl,1) < 2
        ql = [];
    else
        ql = pl(2:end,:);
    end
end

function S_ = Add(cell_, S)
    n = size(S,1);
    m = 0;
    for k = 1 : n
        if isequal(cell_(1,2), S(k,2))
            S(k,1) = {cell2mat(S(k,1)) + cell2mat(cell_(1,1))};
            m = 1;
            break;
        end
    end
    if m == 0
        S(n+1,:) = cell_(1,:);
    end
    S_ = S;
end

% =========================================================================
% 本地非支配排序（备用）
% =========================================================================
function [FrontNo, MaxFNo] = NDSort_Local_Fast(PopObj, nSort)
    [N, M] = size(PopObj);
    FrontNo = inf(1, N);
    MaxFNo = 0;

    while sum(FrontNo < inf) < min(nSort, N)
        MaxFNo = MaxFNo + 1;
        for i = 1 : N
            if FrontNo(i) == inf
                Dominated = false;
                for j = 1 : N
                    if i ~= j && FrontNo(j) == inf
                        m = 1;
                        while m <= M && PopObj(j,m) <= PopObj(i,m)
                            m = m + 1;
                        end
                        if m > M && any(PopObj(j,:) < PopObj(i,:))
                            Dominated = true;
                            break;
                        end
                    end
                end
                if ~Dominated
                    FrontNo(i) = MaxFNo;
                end
            end
        end
    end
end

% =========================================================================
% 提取种群目标
% =========================================================================
function Objs = ExtractPop(Data)
    if isempty(Data)
        Objs = [];
        return;
    end

    if iscell(Data)
        Objs = ExtractPop(Data{end});
        return;
    end

    if isa(Data, 'SOLUTION') || isa(Data, 'class_SOLUTION')
        Objs = Data.objs;
    elseif isstruct(Data) && isfield(Data, 'objs')
        Objs = Data.objs;
    elseif isnumeric(Data)
        Objs = double(Data);
    else
        Objs = [];
    end
end

% =========================================================================
% 提取 PF
% =========================================================================
function Objs = ExtractPF(Data)
    if isempty(Data)
        Objs = [];
        return;
    end

    if iscell(Data)
        TempObjs = [];
        for i = 1:length(Data)
            Chunk = ExtractPF(Data{i});
            if isempty(Chunk)
                continue;
            end
            if ~isempty(TempObjs) && size(Chunk,2) ~= size(TempObjs,2)
                continue;
            end
            TempObjs = [TempObjs; Chunk];
        end
        Objs = TempObjs;
        return;
    end

    if isa(Data, 'SOLUTION') || isa(Data, 'class_SOLUTION')
        Objs = Data.objs;
    elseif isstruct(Data) && isfield(Data, 'objs')
        Objs = Data.objs;
    elseif isnumeric(Data)
        Objs = double(Data);
    else
        Objs = [];
    end
end

function PF2 = DownsamplePF_KeepBoundary(PF, K)
% 降采样时优先保留边界点/极值点
%
% 输入:
%   PF : n x m
%   K  : 目标保留点数
%
% 输出:
%   PF2: <= K x m

    n = size(PF,1);
    m = size(PF,2);

    if n <= K
        PF2 = PF;
        return;
    end

    keep = false(n,1);

    % -------------------------------------------------
    % 1) 强制保留每一维的最小值点和最大值点
    %    对 HV 很关键
    % -------------------------------------------------
    for j = 1:m
        [~, idxMin] = min(PF(:,j));
        [~, idxMax] = max(PF(:,j));
        keep(idxMin) = true;
        keep(idxMax) = true;
    end

    % -------------------------------------------------
    % 2) 再保留一些"离中心最远"的点，增强角点覆盖
    % -------------------------------------------------
    center = mean(PF, 1);
    dist2center = sum((PF - center).^2, 2);
    [~, idxFar] = sort(dist2center, 'descend');

    numFar = min(5*m, n);   % 可调：保留前 5*m 个最远点
    keep(idxFar(1:numFar)) = true;

    fixedIdx = find(keep);

    % 如果保留的边界点已经超过 K，只能从中截取
    if numel(fixedIdx) >= K
        PF2 = PF(fixedIdx(1:K), :);
        PF2 = unique(PF2, 'rows');
        return;
    end

    % -------------------------------------------------
    % 3) 剩余名额再随机补齐
    % -------------------------------------------------
    remain = K - numel(fixedIdx);
    otherIdx = find(~keep);

    rng(1);  % 固定随机种子，可复现
    if numel(otherIdx) > remain
        pick = otherIdx(randperm(numel(otherIdx), remain));
    else
        pick = otherIdx;
    end

    finalIdx = [fixedIdx; pick(:)];
    PF2 = PF(finalIdx, :);

    % 去重
    PF2 = unique(PF2, 'rows');
end