classdef KTDMOEA < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation> <dynamic>
% Knowledge Transfer Dynamic Multi-Objective Evolutionary Algorithm
% theta --- 2 --- Parameter for PS expansion
%
% Modified version:
%   1. Use stage-signature-based prior environment detection.
%   2. Detection logic is consistent with the conservative DTAEA version.
%   3. No longer rely only on Problem.M or Offspring objective dimension.
%   4. If stage changes but objective number M is unchanged, only re-evaluate
%      the population without additional reconstruction.
%   5. Keep original KTDMOEA PSExpansion / PSContraction logic.

%------------------------------- Reference --------------------------------
% G. Ruan, L. L. Minku, S. Menzel, B. Sendhoff, and X. Yao, Knowledge
% Transfer for Dynamic Multi-Objective Optimization With a Changing Number
% of Objectives, IEEE Transactions on Emerging Topics in Computational
% Intelligence, 2024.
%--------------------------------------------------------------------------

    properties
        theta = 2; % Parameter theta as described in Eq. (2)
    end

    methods
        function main(Algorithm, Problem)

            %% Parameter setting
            Algorithm.theta = Algorithm.ParameterSet(2);

            %% Initialize Population
            Population = Problem.Initialization();

            %% Stage-signature prior change detection, consistent with DTAEA
            if isprop(Problem, 'FirstStageGens') && ~isempty(Problem.FirstStageGens)
                firstGens = Problem.FirstStageGens;
            elseif isprop(Problem, 'taut') && ~isempty(Problem.taut)
                firstGens = Problem.taut;
            else
                error('KTDMOEA requires Problem.FirstStageGens or Problem.taut to determine the first change generation.');
            end

            CurrentStage     = getCurrentStageSignature_KTDMOEA(Problem, firstGens);
            CurrentStageInfo = getStageInfo_KTDMOEA(Problem, Population, CurrentStage);

            lastM = CurrentStageInfo.M;

            %% Optimization Loop
            while Algorithm.NotTerminated(Population)

                % =========================================================
                % Step 0: Prior environment detection by stage signature
                % =========================================================
                NextStage = getCurrentStageSignature_KTDMOEA(Problem, firstGens);
                EnvironmentChanged = ~isequal(NextStage, CurrentStage);

                if EnvironmentChanged

                    OldStageInfo = CurrentStageInfo;
                    NewStageInfo = getStageInfo_KTDMOEA(Problem, Population, NextStage);

                    fprintf('[KTDMOEA Change Detected] FE=%d | OldStage=%s | NewStage=%s | OldM=%d | NewM=%d\n', ...
                        Problem.FE, ...
                        stageSigToString_KTDMOEA(CurrentStage), ...
                        stageSigToString_KTDMOEA(NextStage), ...
                        OldStageInfo.M, ...
                        NewStageInfo.M);

                    if NewStageInfo.M > OldStageInfo.M

                        % Objective number increases:
                        % keep KTDMOEA original PS expansion response
                        Population = Algorithm.PSExpansion(Population, Problem, OldStageInfo.M);

                    elseif NewStageInfo.M < OldStageInfo.M

                        % Objective number decreases:
                        % keep KTDMOEA original PS contraction response
                        Population = Algorithm.PSContraction(Population, Problem);

                    else

                        % Same-M stage change:
                        % conservative handling, only re-evaluate population
                        Population = Problem.Evaluation(Population.decs);

                    end

                    CurrentStage     = NextStage;
                    CurrentStageInfo = NewStageInfo;
                    lastM = NewStageInfo.M; %#ok<NASGU>

                    % Skip reproduction in the change-response iteration
                    continue;
                end

                % =========================================================
                % Step 1: Evolutionary Optimization Process
                % =========================================================

                % Calculate fitness for mating selection
                [FrontNo, ~] = NDSort(Population.objs, Population.cons, inf);
                CrowdDis     = CrowdingDistance(Population.objs, FrontNo);

                % Mating Selection
                MatingPool = TournamentSelection(2, Problem.N, FrontNo, CrowdDis);
                Offspring  = OperatorGA(Problem, Population(MatingPool));

                % Population update using KTDMOEA archive update logic
                Population = Algorithm.UpdateCA([Population, Offspring], Problem.N);
            end
        end


        %% Algorithm 1: PS Expansion Strategy
        function NewPop = PSExpansion(Algorithm, OldPop, Problem, oldM)

            % Extract Old Pareto Set
            [FrontNo, ~] = NDSort(OldPop.objs, 1);
            PSt = OldPop(FrontNo == 1);

            if isempty(PSt)
                PSt = OldPop;
            end

            % Step 1: Search Expansion Directions
            Dirs = Algorithm.SearchExpansionDirections(PSt, Problem, OldPop.decs);

            N = Problem.N;
            N_dir = size(Dirs, 1);
            TransferredDecs = [];

            if N_dir > 0

                % Step 2: Calculate N_base
                N_base = floor((N - oldM) / (N_dir * Algorithm.theta));

                if N_base < 1
                    N_base = 1;
                end

                % Evenly select base solutions from PSt using crowding distance
                [~, rank] = sort(CrowdingDistance(PSt.objs), 'descend');
                nSelect = min(length(PSt), N_base);
                BaseSols = PSt(rank(1:nSelect));

                % Step 3: Generate solutions along expansion directions
                Decs = BaseSols.decs;
                NewDecs = [];

                for j = 1 : N_dir

                    D_vec = Dirs(j, :);

                    for i = 1 : size(Decs, 1)

                        x_i = Decs(i, :);

                        for k = 1 : Algorithm.theta

                            % Calculate C
                            diff_upper = (Problem.upper - x_i) ./ D_vec;
                            diff_lower = (Problem.lower - x_i) ./ D_vec;

                            candidates = [diff_upper, diff_lower];
                            candidates(candidates < 1e-6) = inf;

                            C = min(candidates);

                            if isinf(C)
                                C = 1.0;
                            end

                            x_new = x_i + C * rand() * D_vec;
                            NewDecs = [NewDecs; x_new];
                        end
                    end
                end

                TransferredDecs = NewDecs;
            end

            % Step 4: Fill the rest
            nRem = N - size(TransferredDecs, 1);

            if nRem > 0

                SelectPool = OldPop;

                if length(SelectPool) > nRem

                    [~, rank] = sort(CrowdingDistance(SelectPool.objs), 'descend');
                    FillDecs = SelectPool(rank(1:nRem)).decs;

                else

                    nNeed = nRem - length(SelectPool);
                    FillDecs = [SelectPool.decs; ...
                        Problem.lower + rand(nNeed, Problem.D) .* (Problem.upper - Problem.lower)];

                end

                TransferredDecs = [TransferredDecs; FillDecs];
            end

            if size(TransferredDecs, 1) > N
                TransferredDecs = TransferredDecs(1:N, :);
            end

            NewPop = Problem.Evaluation(TransferredDecs);
        end


        %% Algorithm 2: Search Expansion Direction
        function Dirs = SearchExpansionDirections(Algorithm, PSt, Problem, PSt_Decs)

            Dirs = [];

            % Line 1: Find Extreme Points
            [~, maxIdx] = max(PSt.objs, [], 1);
            Pe = PSt(unique(maxIdx));

            if isempty(Pe)
                return;
            end

            % Line 2: Select random extreme point
            xe = Pe(randi(length(Pe)));
            xe_dec = xe.decs;

            % Line 3: Generate detective population around xe
            N = Problem.N;
            RepXe = repmat(xe, N, 1);
            P_var_Decs = Algorithm.PolynomialMutation(RepXe.decs, Problem.lower, Problem.upper);

            % Line 4: Evaluate P_var in new environment
            P_var  = Problem.Evaluation(P_var_Decs);
            PSt_New = Problem.Evaluation(PSt_Decs);

            % Lines 5-8: Filter dominated solutions
            P_non = [];

            for i = 1 : length(P_var)

                dominated = false;

                for j = 1 : length(PSt_New)

                    if all(PSt_New(j).objs <= P_var(i).objs) && ...
                            any(PSt_New(j).objs < P_var(i).objs)

                        dominated = true;
                        break;
                    end
                end

                if ~dominated
                    P_non = [P_non, P_var(i)];
                end
            end

            if isempty(P_non)
                return;
            end

            % Line 9: Density estimation using weight vectors
            currentM = size(P_non(1).objs, 2);
            W = Algorithm.UniformPoint(N, currentM);

            AllObjs = [PSt_New.objs; P_non.objs];

            Zmin = min(AllObjs, [], 1);
            Zmax = max(AllObjs, [], 1);
            Zmax(Zmax == Zmin) = Zmax(Zmax == Zmin) + 1e-6;

            % Identify occupied regions by PSt
            OccupiedRegions = false(1, size(W, 1));
            PSt_Norm = (PSt_New.objs - Zmin) ./ (Zmax - Zmin);

            for i = 1 : length(PSt_New)

                norm_sol = sqrt(sum(PSt_Norm(i, :).^2));

                if norm_sol == 0
                    continue;
                end

                cosine = (PSt_Norm(i, :) * W') ./ norm_sol;
                [~, regionIdx] = max(cosine);

                OccupiedRegions(regionIdx) = true;
            end

            % Filter P_non based on occupied regions
            P_non_Norm = (P_non.objs - Zmin) ./ (Zmax - Zmin);
            KeepIdx = false(1, length(P_non));

            for i = 1 : length(P_non)

                norm_sol = sqrt(sum(P_non_Norm(i, :).^2));

                if norm_sol == 0
                    KeepIdx(i) = true;
                    continue;
                end

                cosine = (P_non_Norm(i, :) * W') ./ norm_sol;
                [~, regionIdx] = max(cosine);

                if ~OccupiedRegions(regionIdx)
                    KeepIdx(i) = true;
                end
            end

            P_non = P_non(KeepIdx);

            if isempty(P_non)
                return;
            end

            % Line 13: Form expansion directions
            Dirs = zeros(length(P_non), Problem.D);

            for i = 1 : length(P_non)

                vec = P_non(i).decs - xe_dec;
                nrm = norm(vec);

                if nrm > 1e-10
                    Dirs(i, :) = vec / nrm;
                else
                    Dirs(i, :) = vec;
                end
            end

            Dirs = unique(Dirs, 'rows');
        end


        %% Algorithm 3: PS Contraction Strategy
        function NewPop = PSContraction(Algorithm, OldPop, Problem)

            % Line 1: Re-evaluate old population
            UpdatedOldPop = Problem.Evaluation(OldPop.decs);

            [FrontNo, ~] = NDSort(UpdatedOldPop.objs, 1);
            P_non = UpdatedOldPop(FrontNo == 1);

            if isempty(P_non)
                P_non = UpdatedOldPop;
            end

            TransferredDecs = P_non.decs;

            % Line 3: Find extreme points
            [~, maxIdx] = max(P_non.objs, [], 1);
            Pe_indices = unique(maxIdx);

            % Lines 4-7: Spread enhancement
            SpreadDecs = [];

            for i = 1 : length(Pe_indices)

                idx_e = Pe_indices(i);
                xe = P_non(idx_e);

                dists = pdist2(xe.objs, P_non.objs);
                dists(dists == 0) = inf;

                [~, idx_close] = min(dists);
                x_close = P_non(idx_close);

                vec = xe.decs - x_close.decs;
                nrm = norm(vec);

                if nrm > 1e-10

                    D = vec / nrm;

                    diff_upper = (Problem.upper - xe.decs) ./ D;
                    diff_lower = (Problem.lower - xe.decs) ./ D;

                    cand = [diff_upper, diff_lower];
                    cand(cand < 1e-6) = inf;

                    C = min(cand);

                    if isinf(C)
                        C = 0;
                    end

                    x_new = xe.decs + C * D;
                    SpreadDecs = [SpreadDecs; x_new];
                end
            end

            TransferredDecs = [TransferredDecs; SpreadDecs];

            % Line 8: Fill
            nNeed = Problem.N - size(TransferredDecs, 1);

            if nNeed > 0

                FillDecs = zeros(nNeed, Problem.D);

                for i = 1 : nNeed

                    p = randperm(length(P_non), 2);

                    xa = P_non(p(1)).decs;
                    xb = P_non(p(2)).decs;

                    x_temp = xa + rand() * (xa - xb);
                    x_temp = max(x_temp, Problem.lower);
                    x_temp = min(x_temp, Problem.upper);

                    FillDecs(i, :) = x_temp;
                end

                TransferredDecs = [TransferredDecs; FillDecs];
            end

            if size(TransferredDecs, 1) > Problem.N
                TransferredDecs = TransferredDecs(1:Problem.N, :);
            end

            NewPop = Problem.Evaluation(TransferredDecs);
        end


        %% Helper: CA Update Mechanism
        function Population = UpdateCA(Algorithm, MixedPop, N) %#ok<INUSD>

            % 1. Remove duplicates in objective space
            [~, uniqueIdx] = unique(MixedPop.objs, 'rows');
            MixedPop = MixedPop(uniqueIdx);

            % 2. Non-dominated sorting
            [FrontNo, ~] = NDSort(MixedPop.objs, 1);

            % 3. Select strictly non-dominated solutions
            Next = FrontNo == 1;
            Population = MixedPop(Next);

            % 4. Truncate if size exceeds N
            if length(Population) > N
                [~, Rank] = sort(CrowdingDistance(Population.objs), 'descend');
                Population = Population(Rank(1:N));
            end
        end


        %% Helper: Polynomial Mutation
        function Offspring = PolynomialMutation(Algorithm, Decs, lower, upper) %#ok<INUSD>

            [N, D] = size(Decs);

            if size(lower, 1) ~= N
                lower = repmat(lower, N, 1);
                upper = repmat(upper, N, 1);
            end

            proM = 1 / D;
            disM = 20;

            Site = rand(N, D) < proM;
            mu   = rand(N, D);

            temp = Site & mu <= 0.5;

            Decs(temp) = Decs(temp) + ...
                (upper(temp) - lower(temp)) .* ...
                ((2 .* mu(temp) + ...
                (1 - 2 .* mu(temp)) .* ...
                (1 - (Decs(temp) - lower(temp)) ./ ...
                (upper(temp) - lower(temp))).^(disM + 1)).^(1 / (disM + 1)) - 1);

            temp = Site & mu > 0.5;

            Decs(temp) = Decs(temp) + ...
                (upper(temp) - lower(temp)) .* ...
                (1 - (2 .* (1 - mu(temp)) + ...
                2 .* (mu(temp) - 0.5) .* ...
                (1 - (upper(temp) - Decs(temp)) ./ ...
                (upper(temp) - lower(temp))).^(disM + 1)).^(1 / (disM + 1)));

            Offspring = max(min(Decs, upper), lower);
        end


        %% Helper: UniformPoint
        function [W, N] = UniformPoint(Algorithm, N, M) %#ok<INUSD>

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
% Stage-signature prior change detection helpers for KTDMOEA
% =========================================================================
function signature = getCurrentStageSignature_KTDMOEA(Problem, FirstStageGens)

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


function StageInfo = getStageInfo_KTDMOEA(Problem, Population, StageSignature)

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


function str = stageSigToString_KTDMOEA(sig)

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