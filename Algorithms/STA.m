classdef STA < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation> <dynamic>
% Similarity Transfer Approach (STA)
%
% Modified version:
%   1. Use stage-signature-based prior environment detection.
%   2. Detection logic is consistent with the conservative DTAEA version.
%   3. No longer rely only on Problem.M or Offspring objective dimension.
%   4. If stage changes but objective number M is unchanged, only re-evaluate
%      the population without additional reconstruction.
%   5. Keep original STA SimilarityTransfer logic for M changes.
%
% -------------------------------------------------------------------------
% Reference:
% G. Ruan, Z. Hou, and X. Yao, "Coping With a Severely Changing Number of
% Objectives in Dynamic Multi-Objective Optimization," IEEE Transactions on
% Evolutionary Computation, 2025.
% -------------------------------------------------------------------------

    properties
        % Archive to store historical populations for different objective numbers
        % Structure: struct('M', {}, 'Pop', {})
        Archive = struct('M', {}, 'Pop', {});

        % STA specific parameters
        gamma = 0.2;  % Proportion of random solutions for diversity enhancement
        theta = 2;    % Parameter for KTDMOEA expansion
    end

    methods
        function main(Algorithm, Problem)

            %% Initialize Population
            Population = Problem.Initialization();

            %% Stage-signature prior change detection, consistent with DTAEA
            if isprop(Problem, 'FirstStageGens') && ~isempty(Problem.FirstStageGens)
                firstGens = Problem.FirstStageGens;
            elseif isprop(Problem, 'taut') && ~isempty(Problem.taut)
                firstGens = Problem.taut;
            else
                error('STA requires Problem.FirstStageGens or Problem.taut to determine the first change generation.');
            end

            CurrentStage     = getCurrentStageSignature_STA(Problem, firstGens);
            CurrentStageInfo = getStageInfo_STA(Problem, Population, CurrentStage);

            lastM = CurrentStageInfo.M;

            %% Save initial population to archive
            Algorithm.UpdateArchive(lastM, Population);

            %% Optimization Loop
            while Algorithm.NotTerminated(Population)

                % =========================================================
                % Step 0: Prior environment detection by stage signature
                % =========================================================
                NextStage = getCurrentStageSignature_STA(Problem, firstGens);
                EnvironmentChanged = ~isequal(NextStage, CurrentStage);

                if EnvironmentChanged

                    OldStageInfo = CurrentStageInfo;
                    NewStageInfo = getStageInfo_STA(Problem, Population, NextStage);

                    fprintf('[STA Change Detected] FE=%d | OldStage=%s | NewStage=%s | OldM=%d | NewM=%d\n', ...
                        Problem.FE, ...
                        stageSigToString_STA(CurrentStage), ...
                        stageSigToString_STA(NextStage), ...
                        OldStageInfo.M, ...
                        NewStageInfo.M);

                    % Save knowledge from the old environment before switching
                    Algorithm.UpdateArchive(OldStageInfo.M, Population);

                    if NewStageInfo.M ~= OldStageInfo.M

                        % Objective number changes:
                        % keep original STA SimilarityTransfer response
                        Population = Algorithm.SimilarityTransfer(Problem, NewStageInfo.M);

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
                % Step 1: Standard Evolutionary Optimization
                % =========================================================
                [FrontNo, ~] = NDSort(Population.objs, Population.cons, inf);
                CrowdDis     = CrowdingDistance(Population.objs, FrontNo);

                OffspringSize = Problem.N;
                RemainingFE = RemainingEvaluationsInCurrentStage_STA(Problem, CurrentStageInfo, firstGens);

                if isfinite(RemainingFE) && RemainingFE < Problem.N
                    OffspringSize = 2 * floor(RemainingFE / 2);

                    if OffspringSize < 2
                        if RemainingFE > 0
                            Problem.Initialization(RemainingFE);
                        end
                        continue;
                    end
                end

                MatingPool = TournamentSelection(2, OffspringSize, FrontNo, CrowdDis);
                Offspring  = OperatorGA(Problem, Population(MatingPool));

                Population = Algorithm.EnvironmentalSelection([Population, Offspring], Problem.N);
            end
        end


        %% Method: Update Archive
        function UpdateArchive(Algorithm, M, Pop)

            [FrontNo, ~] = NDSort(Pop.objs, 1);
            BestSols = Pop(FrontNo == 1);

            if isempty(BestSols)
                BestSols = Pop;
            end

            newIdx = length(Algorithm.Archive) + 1;
            Algorithm.Archive(newIdx).M = M;
            Algorithm.Archive(newIdx).Pop = BestSols;
            return;
        end


        %% Method: Similarity Transfer Strategy
        function NewPop = SimilarityTransfer(Algorithm, Problem, currentM)

            if nargin < 3 || isempty(currentM)
                currentM = Problem.M;
            end

            if isempty(Algorithm.Archive)
                NewPop = Problem.Initialization();
                return;
            end

            HistoryMs = [Algorithm.Archive.M];
            Distances = abs(HistoryMs - currentM);
            minDist = min(Distances);
            candidates = find(Distances == minDist);

            % Select the first candidate
            bestIdx = candidates(1);

            SourceEnv = Algorithm.Archive(bestIdx);
            SourcePop = SourceEnv.Pop;

            M_tr = SourceEnv.M;
            M_c  = currentM;

            TransferredDecs = [];

            if M_tr < M_c

                % Case A: Increasing NObj -> Use KTDMOEA Expansion + STA Randomization
                [TransferredDecs, N_exp] = Algorithm.KTDMOEA_Expansion(SourcePop, Problem, M_tr);

                % STA Strategy: Randomization Enhancing Diversity
                threshold = (M_c - M_tr) / 2;

                if N_exp < threshold

                    N = Problem.N;
                    nRandom = round(Algorithm.gamma * N);

                    if nRandom > 0

                        RandomDecs = repmat(Problem.lower, nRandom, 1) + ...
                            rand(nRandom, Problem.D) .* ...
                            (repmat(Problem.upper, nRandom, 1) - repmat(Problem.lower, nRandom, 1));

                        if size(TransferredDecs, 1) >= N
                            TransferredDecs(end - nRandom + 1:end, :) = RandomDecs;
                        else
                            TransferredDecs = [TransferredDecs; RandomDecs];
                        end
                    end
                end

            elseif M_tr > M_c

                % Case B: Decreasing NObj -> Use LEC Contraction
                TransferredDecs = Algorithm.LEC_Contraction(SourcePop, Problem, M_tr);

            else

                % Case C: Same NObj -> Direct Copy
                TransferredDecs = SourcePop.decs;
                TransferredDecs = Algorithm.AdjustDecsDimension(TransferredDecs, Problem.D, Problem.lower, Problem.upper);

            end

            N = Problem.N;
            [nGen, dGen] = size(TransferredDecs);

            if nGen < N

                nNeed = N - nGen;
                RandomDecs = repmat(Problem.lower, nNeed, 1) + ...
                    rand(nNeed, dGen) .* ...
                    (repmat(Problem.upper, nNeed, 1) - repmat(Problem.lower, nNeed, 1));

                TransferredDecs = [TransferredDecs; RandomDecs];

            elseif nGen > N

                perm = randperm(nGen);
                TransferredDecs = TransferredDecs(perm(1:N), :);

            end

            NewPop = Problem.Evaluation(TransferredDecs);
        end


        %% Helper 1: KTDMOEA Expansion Logic
        function [NewDecs, N_exp] = KTDMOEA_Expansion(Algorithm, SourcePop, Problem, M_tr)

            % 1. Identify expansion directions following KTDMOEA:
            % randomly select one extreme point and search around it.
            [Dirs, ~] = Algorithm.SearchExpansionDirections(SourcePop, Problem);
            N_exp = size(Dirs, 1);

            if N_exp == 0
                NewDecs = Algorithm.AdjustDecsDimension(SourcePop.decs, Problem.D, Problem.lower, Problem.upper);
                return;
            end

            % 2. Generate Solutions along Directions
            N = Problem.N;
            N_dir = size(Dirs, 1);

            N_base = floor((N - M_tr) / (N_dir * Algorithm.theta));

            if N_base < 1
                N_base = 1;
            end

            [~, rank] = sort(CrowdingDistance(SourcePop.objs), 'descend');
            BaseSols = SourcePop(rank(1:min(length(SourcePop), N_base)));

            NewDecs = [];
            BaseDecs = Algorithm.AdjustDecsDimension(BaseSols.decs, Problem.D, Problem.lower, Problem.upper);

            for i = 1 : size(BaseDecs, 1)

                x = BaseDecs(i, :);

                for j = 1 : N_dir

                    D_vec = Dirs(j, :);

                    for k = 1 : Algorithm.theta

                        x_new = Algorithm.GenerateSolution(x, D_vec, Problem.lower, Problem.upper);
                        NewDecs = [NewDecs; x_new];

                    end
                end

                if size(NewDecs, 1) >= N
                    break;
                end
            end

            if size(NewDecs, 1) < N

                nRem = N - size(NewDecs, 1);
                FillDecs = BaseDecs(1:min(nRem, size(BaseDecs, 1)), :);

                if size(FillDecs, 1) < nRem

                    nNeed = nRem - size(FillDecs, 1);
                    RandFill = repmat(Problem.lower, nNeed, 1) + ...
                        rand(nNeed, Problem.D) .* ...
                        (repmat(Problem.upper, nNeed, 1) - repmat(Problem.lower, nNeed, 1));

                    FillDecs = [FillDecs; RandFill];
                end

                NewDecs = [NewDecs; FillDecs];
            end
        end


        %% Helper 2: LEC Contraction Logic
        function NewDecs = LEC_Contraction(Algorithm, SourcePop, Problem, M_tr)

            BaseDecs = Algorithm.AdjustDecsDimension(SourcePop.decs, Problem.D, Problem.lower, Problem.upper);

            C_con = [];

            try

                [Coeff, ~, ~] = pca(BaseDecs);
                nEig = min(size(Coeff, 2), M_tr - 1);

                if nEig > 0
                    C_con = Coeff(:, 1:nEig)';
                end

            catch

                C_con = [];

            end

            if isempty(C_con)
                C_con = 2 * rand(1, Problem.D) - 1;
            end

            % Selection of most promising contraction directions
            D_con = [];

            SourceInTarget = Problem.Evaluation(BaseDecs);
            idx = randi(length(SourceInTarget));
            x_base = SourceInTarget(idx);

            for i = 1 : size(C_con, 1)

                D_vec = C_con(i, :);

                nGen = floor(Problem.N / max(1, (M_tr - 1)));

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

            if isempty(D_con)
                D_con = C_con;
            end

            % Execution
            lambda_param = floor(Problem.N / 20);

            if lambda_param < 1
                lambda_param = 1;
            end

            perm = randperm(size(BaseDecs, 1));
            P_base = BaseDecs(perm(1:min(size(BaseDecs, 1), lambda_param)), :);

            NewDecs = [];
            N_con = size(D_con, 1);
            Ng = floor(Problem.N / (N_con * size(P_base, 1)));

            if Ng < 1
                Ng = 1;
            end

            for i = 1 : size(P_base, 1)

                x_i = P_base(i, :);

                for j = 1 : N_con

                    D_vec = D_con(j, :);

                    for k = 1 : Ng

                        x_new = Algorithm.GenerateSolution(x_i, D_vec, Problem.lower, Problem.upper);
                        NewDecs = [NewDecs; x_new];

                    end
                end
            end

            if size(NewDecs, 1) < Problem.N

                nRem = Problem.N - size(NewDecs, 1);
                indices = randi(size(BaseDecs, 1), nRem, 1);
                Fill = BaseDecs(indices, :);

                NewDecs = [NewDecs; Fill];
            end
        end


        %% Utility: KTDMOEA-style Search of Expansion Directions
        function [Dirs, Valid] = SearchExpansionDirections(Algorithm, SourcePop, Problem)

            Dirs = [];
            Valid = false;

            if isempty(SourcePop)
                return;
            end

            SourceDecs = Algorithm.AdjustDecsDimension(SourcePop.decs, Problem.D, Problem.lower, Problem.upper);

            [~, maxIdx] = max(SourcePop.objs, [], 1);
            ExtremeDecs = SourceDecs(unique(maxIdx), :);

            if isempty(ExtremeDecs)
                return;
            end

            xe = ExtremeDecs(randi(size(ExtremeDecs, 1)), :);

            N = Problem.N;
            D = Problem.D;
            lower = Problem.lower;
            upper = Problem.upper;

            if size(lower, 1) > 1
                lower = lower';
            end

            if size(upper, 1) > 1
                upper = upper';
            end

            P_var_Decs = repmat(xe, N, 1);
            LowerRep = repmat(lower, N, 1);
            UpperRep = repmat(upper, N, 1);

            Site = rand(N, D) < 1 / D;
            mu = rand(N, D);

            temp = Site & mu <= 0.5;
            P_var_Decs(temp) = P_var_Decs(temp) + ...
                (UpperRep(temp) - LowerRep(temp)) .* ...
                ((2 .* mu(temp) + ...
                (1 - 2 .* mu(temp)) .* ...
                (1 - (P_var_Decs(temp) - LowerRep(temp)) ./ ...
                (UpperRep(temp) - LowerRep(temp))).^21).^(1 / 21) - 1);

            temp = Site & mu > 0.5;
            P_var_Decs(temp) = P_var_Decs(temp) + ...
                (UpperRep(temp) - LowerRep(temp)) .* ...
                (1 - (2 .* (1 - mu(temp)) + ...
                2 .* (mu(temp) - 0.5) .* ...
                (1 - (UpperRep(temp) - P_var_Decs(temp)) ./ ...
                (UpperRep(temp) - LowerRep(temp))).^21).^(1 / 21));

            P_var_Decs = max(min(P_var_Decs, UpperRep), LowerRep);

            P_var = Problem.Evaluation(P_var_Decs);
            SourceInTarget = Problem.Evaluation(SourceDecs);

            IsDominated = false(length(P_var), 1);
            TargetObjs = SourceInTarget.objs;
            ScoutObjs = P_var.objs;

            for i = 1 : length(P_var)
                dominators = all(TargetObjs <= ScoutObjs(i, :), 2) & ...
                    any(TargetObjs < ScoutObjs(i, :), 2);
                IsDominated(i) = any(dominators);
            end

            KeepIdx = find(~IsDominated);

            if isempty(KeepIdx)
                return;
            end

            currentM = size(TargetObjs, 2);
            [W, ~] = Algorithm.UniformPoint(N, currentM);

            AllObjs = [TargetObjs; ScoutObjs(KeepIdx, :)];
            Zmin = min(AllObjs, [], 1);
            Zmax = max(AllObjs, [], 1);
            Zmax(Zmax == Zmin) = Zmax(Zmax == Zmin) + 1e-6;

            SourceNorm = (TargetObjs - Zmin) ./ (Zmax - Zmin);
            ScoutNorm = (ScoutObjs(KeepIdx, :) - Zmin) ./ (Zmax - Zmin);

            OccupiedRegions = false(1, size(W, 1));

            for i = 1 : size(SourceNorm, 1)
                regionIdx = Algorithm.AssociateByPerpendicularDistance(SourceNorm(i, :), W);
                OccupiedRegions(regionIdx) = true;
            end

            FinalKeep = false(size(KeepIdx));

            for i = 1 : length(KeepIdx)
                regionIdx = Algorithm.AssociateByPerpendicularDistance(ScoutNorm(i, :), W);

                if ~OccupiedRegions(regionIdx)
                    FinalKeep(i) = true;
                end
            end

            ValidIdx = KeepIdx(FinalKeep);

            if isempty(ValidIdx)
                return;
            end

            RawDirs = P_var_Decs(ValidIdx, :) - repmat(xe, length(ValidIdx), 1);
            len = sqrt(sum(RawDirs.^2, 2));
            Dirs = RawDirs ./ (len + 1e-10);
            Dirs = unique(Dirs, 'rows');
            Valid = ~isempty(Dirs);
        end


        %% Utility: Associate by perpendicular distance to weight vectors
        function regionIdx = AssociateByPerpendicularDistance(~, sol, W)

            if all(abs(sol) < 1e-12)
                regionIdx = 1;
                return;
            end

            WNorm2 = sum(W.^2, 2);
            WNorm2(WNorm2 < 1e-10) = 1e-10;
            scalar = (W * sol') ./ WNorm2;
            Proj = W .* repmat(scalar, 1, size(W, 2));
            Dist = sqrt(sum((Proj - repmat(sol, size(W, 1), 1)).^2, 2));
            [~, regionIdx] = min(Dist);
        end


        %% Utility: Find Expansion Directions
        function [Dirs, Valid] = GetExpansionDirs(Algorithm, SourcePop, Problem)

            Dirs = [];
            Valid = false;

            % 1. Prepare Data
            D = Problem.D;
            lower = Problem.lower;
            upper = Problem.upper;

            SourceDecs = Algorithm.AdjustDecsDimension(SourcePop.decs, D, lower, upper);
            SourceInTarget = Problem.Evaluation(SourceDecs);

            % 2. Find Extreme Points
            ObjsOrigin = SourcePop.objs;
            ExtremePointsDecs = [];

            for i = 1 : size(ObjsOrigin, 2)

                [~, rank] = sort(ObjsOrigin(:, i), 'descend');
                ExtremePointsDecs = [ExtremePointsDecs; SourceDecs(rank(1), :)];

            end

            ExtremePointsDecs = unique(ExtremePointsDecs, 'rows');

            % 3. Generate Scout Solutions via Mutation
            P_var_Decs = [];
            Xe_List = [];

            if size(lower, 1) > 1
                lower = lower';
            end

            if size(upper, 1) > 1
                upper = upper';
            end

            for i = 1 : size(ExtremePointsDecs, 1)

                xe = ExtremePointsDecs(i, :);

                N_mutants = 100;
                Mus = repmat(xe, N_mutants, 1);

                LowerRep = repmat(lower, N_mutants, 1);
                UpperRep = repmat(upper, N_mutants, 1);

                Site = rand(size(Mus)) < 1 / D;
                mu = rand(size(Mus));

                temp = Site & mu <= 0.5;

                Mus(temp) = Mus(temp) + ...
                    (UpperRep(temp) - LowerRep(temp)) .* ...
                    ((2 .* mu(temp) + ...
                    (1 - 2 .* mu(temp)) .* ...
                    (1 - (Mus(temp) - LowerRep(temp)) ./ ...
                    (UpperRep(temp) - LowerRep(temp))).^21).^(1 / 21) - 1);

                temp = Site & mu > 0.5;

                Mus(temp) = Mus(temp) + ...
                    (UpperRep(temp) - LowerRep(temp)) .* ...
                    (1 - (2 .* (1 - mu(temp)) + ...
                    2 .* (mu(temp) - 0.5) .* ...
                    (1 - (UpperRep(temp) - Mus(temp)) ./ ...
                    (UpperRep(temp) - LowerRep(temp))).^21).^(1 / 21));

                Mus = max(min(Mus, UpperRep), LowerRep);

                P_var_Decs = [P_var_Decs; Mus];
                Xe_List = [Xe_List; repmat(xe, N_mutants, 1)];
            end

            if isempty(P_var_Decs)
                return;
            end

            % 4. Evaluate Scout Solutions in Target Environment
            P_var_Pop = Problem.Evaluation(P_var_Decs);

            % 5. Filter Dominated Solutions
            IsDominated = false(length(P_var_Pop), 1);

            TargetObjs = SourceInTarget.objs;
            ScoutObjs  = P_var_Pop.objs;

            for i = 1 : length(P_var_Pop)

                currentObj = ScoutObjs(i, :);
                dominators = all(TargetObjs <= currentObj, 2) & ...
                    any(TargetObjs < currentObj, 2);

                if any(dominators)
                    IsDominated(i) = true;
                end
            end

            KeepIdx = find(~IsDominated);

            if ~isempty(KeepIdx)

                % Generate weight vectors based on current objective number
                currentM = size(TargetObjs, 2);
                [W, ~] = Algorithm.UniformPoint(Problem.N, currentM);

                % Normalize all known objectives
                AllObjs = [TargetObjs; ScoutObjs(KeepIdx, :)];

                Zmin = min(AllObjs, [], 1);
                Zmax = max(AllObjs, [], 1);
                Zmax(Zmax == Zmin) = Zmax(Zmax == Zmin) + 1e-6;

                % Identify occupied regions by source solutions
                OccupiedRegions = false(1, size(W, 1));
                SourceNorm = (TargetObjs - Zmin) ./ (Zmax - Zmin);

                for i = 1 : size(SourceNorm, 1)

                    norm_sol = norm(SourceNorm(i, :));

                    if norm_sol == 0
                        continue;
                    end

                    cosine = (SourceNorm(i, :) * W') / norm_sol;
                    [~, regionIdx] = max(cosine);

                    OccupiedRegions(regionIdx) = true;
                end

                % Filter P_var based on occupied regions
                FinalKeep = false(size(KeepIdx));
                ScoutNorm = (ScoutObjs(KeepIdx, :) - Zmin) ./ (Zmax - Zmin);

                for i = 1 : length(KeepIdx)

                    norm_sol = norm(ScoutNorm(i, :));

                    if norm_sol == 0
                        FinalKeep(i) = true;
                        continue;
                    end

                    cosine = (ScoutNorm(i, :) * W') / norm_sol;
                    [~, regionIdx] = max(cosine);

                    if ~OccupiedRegions(regionIdx)
                        FinalKeep(i) = true;
                    end
                end

                ValidIdx = KeepIdx(FinalKeep);

            else

                ValidIdx = [];

            end

            % 7. Form Directions
            if ~isempty(ValidIdx)

                ValidPVar = P_var_Decs(ValidIdx, :);
                ValidXe   = Xe_List(ValidIdx, :);

                RawDirs = ValidPVar - ValidXe;
                len = sqrt(sum(RawDirs.^2, 2));
                Dirs = RawDirs ./ (len + 1e-10);

                Dirs = unique(Dirs, 'rows');
                Valid = true;
            end
        end


        %% Utility: Generate Solution with Step Size
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

            para(maskPos) = (upper(maskPos) - x(maskPos)) ./ D_vec(maskPos);
            para(maskNeg) = (lower(maskNeg) - x(maskNeg)) ./ D_vec(maskNeg);

            ss = min(para);

            if ss < 0
                ss = 0;
            end

            x_new = x + ss * rand() * D_vec;
            x_new = max(min(x_new, upper), lower);
        end


        %% Utility: Adjust Dimension
        function NewDecs = AdjustDecsDimension(~, OldDecs, TargetD, lower, upper)

            [N, OldD] = size(OldDecs);

            NewDecs = zeros(N, TargetD);
            minD = min(OldD, TargetD);

            NewDecs(:, 1:minD) = OldDecs(:, 1:minD);

            if isscalar(lower)
                lower = repmat(lower, 1, TargetD);
            end

            if isscalar(upper)
                upper = repmat(upper, 1, TargetD);
            end

            if TargetD > OldD

                low_exp = repmat(lower(OldD + 1:end), N, 1);
                upp_exp = repmat(upper(OldD + 1:end), N, 1);

                NewDecs(:, OldD + 1:end) = ...
                    low_exp + rand(N, TargetD - OldD) .* (upp_exp - low_exp);
            end
        end


        %% Helper: Environmental Selection
        function Population = EnvironmentalSelection(~, Population, N)

            [FrontNo, MaxFNo] = NDSort(Population.objs, Population.cons, N);

            Next = FrontNo < MaxFNo;
            Last = find(FrontNo == MaxFNo);

            [~, Rank] = sort(CrowdingDistance(Population(Last).objs), 'descend');

            Next(Last(Rank(1:N - sum(Next)))) = true;

            Population = Population(Next);
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
% Stage-signature prior change detection helpers for STA
% =========================================================================
function signature = getCurrentStageSignature_STA(Problem, FirstStageGens)

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


function StageInfo = getStageInfo_STA(Problem, Population, StageSignature)

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


function str = stageSigToString_STA(sig)

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


function RemainingFE = RemainingEvaluationsInCurrentStage_STA(Problem, CurrentStageInfo, FirstStageGens)

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
    nextChangeFE  = nextChangeGen * Problem.N;
    RemainingFE   = nextChangeFE - Problem.FE;

    if RemainingFE < 0
        RemainingFE = 0;
    end
end
