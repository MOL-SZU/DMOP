function [score, HV_Pop, HV_PF] = MHV_Strict(Population, Problem)
% MHV_Strict
%
% 更稳妥版本：
%   1. 保留 PF_HV 缓存；
%   2. 但缓存签名不再只依赖 M 和 PFCount；
%   3. 签名加入 Problem 类名、当前环境目标、NormPF 的统计特征和代表点；
%   4. 避免不同问题 / 不同环境 / 不同 PF 被错误复用。
%
% 输出:
%   score  = HV_Pop / HV_PF
%   HV_Pop = 当前种群 HV
%   HV_PF  = 当前环境 PF 的 HV

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

    catch ME
        fprintf('[MHV Error 2] %s\n', ME.message);
        return;
    end

    [N_pop, M] = size(PopObj);
    [N_pf, ~]  = size(PF);

    %% 3. 归一化
    try
        [z_star, z_nad] = CalTruePFBounds_DMixed(Problem, M);
        if isempty(z_star) || isempty(z_nad)
            z_star = min(PF, [], 1);
            z_nad  = max(PF, [], 1);
        end

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

        % ---------- Pop HV ----------
        if isempty(FilteredPop)
            HV_Pop = 0;
        else
            HV_Pop = CalHV_Hybrid(FilteredPop, RefPoint, true);
        end

        % ---------- PF HV，使用更稳妥缓存 ----------
        TruePFHV = CalTrueNormalizedPFHV_DMixed(Problem, RefPoint);

        if ~isempty(TruePFHV)
            HV_PF = TruePFHV;
        elseif isempty(FilteredPF)
            HV_PF = 0;
            Last_PF_Signature = [];
        else
            Current_Signature = GetPFSignature_Strict(Problem, M, FilteredPF);

            if ~isempty(Last_PF_Signature) && isequal(Current_Signature, Last_PF_Signature)
                HV_PF = Last_HV_PF;
                % 如需调试缓存命中，可打开下面一行
                % fprintf('[MHV Cache Hit] %s | HV_PF = %.10f\n', Current_Signature, HV_PF);
            else
                HV_PF = CalHV_Hybrid(FilteredPF, RefPoint, true);

                Last_HV_PF = HV_PF;
                Last_PF_Signature = Current_Signature;

                % 如需调试缓存刷新，可打开下面一行
                % fprintf('[MHV Cache Update] %s | HV_PF = %.10f\n', Current_Signature, HV_PF);
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
% 更稳妥的 PF 签名
% =========================================================================
function [z_star, z_nad] = CalTruePFBounds_DMixed(Problem, M)
    z_star = [];
    z_nad  = [];

    try
        if ~isa(Problem, 'D_MIXED_SEQ') || isempty(Problem.ObjGroups)
            return;
        end

        CurrentGen = floor(max(0, Problem.FE - 1) / Problem.N);
        if CurrentGen < Problem.FirstStageGens
            t = 0;
        else
            t = 1 + floor((CurrentGen - Problem.FirstStageGens) / Problem.taut);
        end
        idx = min(t + 1, length(Problem.ObjGroups));
        realIndices = Problem.ObjGroups{idx};
        if length(realIndices) ~= M
            return;
        end

        pName = Problem.ObjMap{realIndices(1)};
        if any(~strcmp(Problem.ObjMap(realIndices), pName))
            return;
        end

        if strcmp(pName, 'DTLZ1')
            k = Problem.D - Problem.GlobalMaxM + 1;
            gMax = 100 * k * (1 + 1.202602641454022);
            z_star = ones(1, M) * (-0.5 * (1 + gMax));
            z_nad  = zeros(1, M);
        elseif any(strcmp(pName, {'DTLZ2','DTLZ4'}))
            z_star = ones(1, M) * (-(1 + 0.5^2 * 10));
            z_nad  = zeros(1, M);
        elseif strcmp(pName, 'DTLZ3')
            k = Problem.D - Problem.GlobalMaxM + 1;
            gMax = 100 * k * (1 + 1.202602641454022);
            z_star = ones(1, M) * (-(1 + gMax));
            z_nad  = zeros(1, M);
        elseif startsWith(pName, 'WFG')
            pIndex = str2double(pName(4:end));
            max_xM = 1;
            if any(pIndex == [6, 9])
                l = Problem.WFG_l;
                c = ceil(l / 2);
                max_xM = l^2 / (c * (1 + 2*l - 2*c));
            end
            scales = Problem.WFG_S(realIndices);
            z_star = -(max_xM + scales);
            z_nad  = ones(1, M) * (-max_xM);
        end
    catch
        z_star = [];
        z_nad  = [];
    end
end


function Volume = CalTrueNormalizedPFHV_DMixed(Problem, RefPoint)
    Volume = [];

    try
        if ~isa(Problem, 'D_MIXED_SEQ') || isempty(Problem.CurrentIndices)
            return;
        end

        pName = Problem.ObjMap{Problem.CurrentIndices(1)};
        M = length(RefPoint);
        if any(abs(RefPoint - RefPoint(1)) > 1e-12)
            return;
        end

        margin = RefPoint(1) - 1;
        if margin < 0
            return;
        end

        Volume = 0;
        for k = 0:M
            if strcmp(pName, 'DTLZ1')
                % y = 1-h, h >= 0, sum(h) = 1
                unitSection = 1 / factorial(k);
            elseif any(strcmp(pName, {'DTLZ2','DTLZ3','DTLZ4'})) || startsWith(pName, 'WFG')
                % y = 1-h, h >= 0, norm(h,2) = 1
                unitSection = pi^(k/2) / gamma(k/2 + 1) / 2^k;
            else
                Volume = [];
                return;
            end
            Volume = Volume + nchoosek(M, k) * margin^(M-k) * unitSection;
        end
    catch
        Volume = [];
    end
end


function Sig = GetPFSignature_Strict(Problem, M, FilteredPF)
% GetPFSignature_Strict
%
% 原来的签名：
%   Sig = sprintf('M=%d|PF=%d', M, PFCount);
%
% 问题：
%   不同问题、不同环境、不同 PF，只要 M 和 PFCount 一样，就可能错误复用 HV_PF。
%
% 新签名加入：
%   1. Problem 类名；
%   2. 当前目标数 M；
%   3. PF 点数；
%   4. 当前环境目标索引，如果 Problem 有 CurrentIndices；
%   5. PF 的 min / max / mean / std；
%   6. PF 排序后的若干个代表点。

    PFCount = size(FilteredPF, 1);

    % ------------------------------
    % 1. Problem 类名
    % ------------------------------
    problemClass = class(Problem);

    % ------------------------------
    % 2. 当前环境目标索引
    % ------------------------------
    currentIndexStr = 'UnknownIndices';

    try
        if isprop(Problem, 'CurrentIndices') && ~isempty(Problem.CurrentIndices)
            currentIndexStr = mat2str(Problem.CurrentIndices);
        end
    catch
        currentIndexStr = 'UnknownIndices';
    end

    % ------------------------------
    % 3. 当前 FE / N 信息
    % ------------------------------
    stageStr = 'UnknownStage';

    try
        if isprop(Problem, 'FE') && isprop(Problem, 'N')
            CurrentGen = floor(max(0, Problem.FE - 1) / Problem.N);
            stageStr = sprintf('Gen=%d|FE=%d|N=%d', CurrentGen, Problem.FE, Problem.N);
        end
    catch
        stageStr = 'UnknownStage';
    end

    % ------------------------------
    % 4. PF 统计特征
    % ------------------------------
    pfMin  = min(FilteredPF, [], 1);
    pfMax  = max(FilteredPF, [], 1);
    pfMean = mean(FilteredPF, 1);
    pfStd  = std(FilteredPF, 0, 1);

    % 四舍五入，避免极小浮点误差导致缓存完全失效
    pfMin  = round(pfMin,  10);
    pfMax  = round(pfMax,  10);
    pfMean = round(pfMean, 10);
    pfStd  = round(pfStd,  10);

    % ------------------------------
    % 5. PF 代表点
    % ------------------------------
    SortedPF = sortrows(FilteredPF);

    sampleNum = min(10, size(SortedPF, 1));

    if sampleNum <= 1
        sampleIdx = 1;
    else
        sampleIdx = round(linspace(1, size(SortedPF, 1), sampleNum));
    end

    samplePF = SortedPF(sampleIdx, :);
    samplePF = round(samplePF, 10);

    % ------------------------------
    % 6. 合成签名
    % ------------------------------
    Sig = sprintf(['Class=%s|' ...
                   'M=%d|' ...
                   'PFCount=%d|' ...
                   'Indices=%s|' ...
                   '%s|' ...
                   'Min=%s|' ...
                   'Max=%s|' ...
                   'Mean=%s|' ...
                   'Std=%s|' ...
                   'Samples=%s'], ...
                   problemClass, ...
                   M, ...
                   PFCount, ...
                   currentIndexStr, ...
                   stageStr, ...
                   mat2str(pfMin), ...
                   mat2str(pfMax), ...
                   mat2str(pfMean), ...
                   mat2str(pfStd), ...
                   mat2str(samplePF));
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
            Volume = Volume + ...
                (SortedP(i+1,1) - SortedP(i,1)) * ...
                (RefPoint(2) - SortedP(i,2));
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

    Samples = unifrnd(repmat(MinValue, SampleNum, 1), ...
                      repmat(MaxValue, SampleNum, 1));

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

    Volume = prod(MaxValue - MinValue) * ...
             (1 - size(Samples,1) / SampleNum);
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
% 本地非支配排序，备用
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


% =========================================================================
% 降采样时优先保留边界点/极值点
% 当前主流程中没有启用，但保留备用
% =========================================================================
function PF2 = DownsamplePF_KeepBoundary(PF, K)

    n = size(PF,1);
    m = size(PF,2);

    if n <= K
        PF2 = PF;
        return;
    end

    keep = false(n,1);

    % 1. 保留每一维的最小值点和最大值点
    for j = 1:m
        [~, idxMin] = min(PF(:,j));
        [~, idxMax] = max(PF(:,j));

        keep(idxMin) = true;
        keep(idxMax) = true;
    end

    % 2. 保留离中心最远的点
    center = mean(PF, 1);
    dist2center = sum((PF - center).^2, 2);

    [~, idxFar] = sort(dist2center, 'descend');

    numFar = min(5*m, n);
    keep(idxFar(1:numFar)) = true;

    fixedIdx = find(keep);

    if numel(fixedIdx) >= K
        PF2 = PF(fixedIdx(1:K), :);
        PF2 = unique(PF2, 'rows');
        return;
    end

    % 3. 剩余名额随机补齐
    remain = K - numel(fixedIdx);
    otherIdx = find(~keep);

    rng(1);

    if numel(otherIdx) > remain
        pick = otherIdx(randperm(numel(otherIdx), remain));
    else
        pick = otherIdx;
    end

    finalIdx = [fixedIdx; pick(:)];

    PF2 = PF(finalIdx, :);
    PF2 = unique(PF2, 'rows');
end
