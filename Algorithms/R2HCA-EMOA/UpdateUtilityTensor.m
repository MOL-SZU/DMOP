function tensor = UpdateUtilityTensor(Problem,PopObj,worst,r,tensor,W,row)
% Utility tensor update

%--------------------------------------------------------------------------
% Author: Ke Shang
% Version: 0.1
% The paper has been submitted to TEVC.
% The code uses PlatEMO published in "Ye Tian, Ran Cheng, Xingyi Zhang, 
% and Yaochu Jin, PlatEMO: A MATLAB Platform for Evolutionary 
% Multi-Objective Optimization [Educational Forum], IEEE Computational 
% Intelligence Magazine, 2017, 12(4): 73-87".
%-------------------------------------------------------------------------- 

    s = PopObj(worst,:);

    OriginalPop = PopObj;
    r = r*max(PopObj,[],1);
    PopObj(worst,:) = [];
    m = min([size(OriginalPop,2),size(W,2),length(s),length(r)]);
    if m <= 0
        return;
    end
    s = s(:,1:m);
    r = r(:,1:m);
    OriginalPop = OriginalPop(:,1:m);
    PopObj = PopObj(:,1:m);
    
    utilityPoints = zeros(row,Problem.N+1);
    lastPoints = zeros(row,Problem.N+1);
    for j = 1:row
        currWeights = repmat(max(W(j,1:m),1e-12), Problem.N+1,1);
        utilityPoints(j,:) = max((s-OriginalPop)./currWeights,[],2)';
        lastPoints(j,:) = min(abs(r-OriginalPop)./currWeights,[],2)';
    end
    
    if worst == 1
        tensor(worst,:,:) = utilityPoints;
    elseif worst == Problem.N+1
        tensor(worst-1,:,:) = utilityPoints;    
    else
        tensor(worst-1,:,1:worst-1) = utilityPoints(:,1:worst-1);
        tensor(worst,:,worst+1:end) = utilityPoints(:,worst+1:end);
    end
    tensor(end,:,:) = lastPoints;
    
    utilityPoints = zeros(Problem.N,row);
    for j = 1:row
        currWeights = repmat(max(W(j,1:m),1e-12), Problem.N, 1);
        utilityPoints(:,j) = max((PopObj-s)./currWeights,[],2);
    end
    
    tensor(1:end-1,:,worst) = utilityPoints;
end
