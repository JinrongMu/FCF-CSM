function [varargout] = DPC_v2(data, K, varargin)
% DPC_v2 - Density Peak Clustering (DPC) V2.
% 
% Version V2: Select the cluster center according to the number of
% clusters, the number of cluster centers is less than the number of clusters.
% 
% Syntax
% =================
% clusterLabels = DPC_v2(data)
% [clusterLabels, centerIdxs] = DPC_v2(data)
% [clusterLabels, centerIdxs, rho, delta] = DPC_v2(data,'percent',4,'kernel','gaussian','debug_mode',1)
%
% Input Arguments
% =================
% data          Data to be clustered, NxM matrix.
% K             Number of clusters, K=0 for automatic clustering.
%
% percent       Percentage used to select dc distance, default is 4.
% kernel        The type of kernel to calculate the density rho, ['gaussian','cutoff'].
% debug_mode    Control debug information, 0: Silent, 1: Call information, 2: Call details.	
%
% Output Arguments
% =================
% out           L or {L, idxs} or {L, idxs, rho, delta}.
%
% clusterLabels	Clustering results, i.e., clustering labels.
% centerIdxs    The index of clustering centers.
% rho           Density matrix, N x 1 vector.
% delta         Minimum distance matrix, N x 1 vector.
% 
% References
% =================
% [1] A. Rodriguez, and A. Laio, 
% "Clustering by fast search and find of density peaks," 
% Science, vol. 344, no. 6191, pp. 1492-1496, 2014.
% [2] https://blog.csdn.net/qq_37055672/article/details/130000567
%
% Parameter Initialization
% =========================================================
arg = inputParser; fun_name = 'DPC_v2';                        
addParameter(arg,'percent',4);
addParameter(arg,'kernel','gaussian'); 
addParameter(arg,'debug_mode',0); 
parse(arg,varargin{:});

if arg.Results.debug_mode == 1
    fprintf('\nCall functions:\t%s\n', fun_name)  
elseif arg.Results.debug_mode == 2     
    fprintf('\nCall functions:\t%s\n', fun_name)
    fprintf('----------------------------------------');
    fprintf('\nDefault Parameters:\n'); disp(arg.Results); 
end
    
% Method Implementation 
% =========================================================

% 计算距离矩阵
distances = pdist2(data, data); 

% 计算邻域截断距离
percent = arg.Results.percent;
dc = computeNeighborhoodCutoff(distances, percent);

% 计算局部密度
kernel = arg.Results.kernel; 
rho = computeLocalDensity(distances, dc, kernel);

% 计算最小距离和最近邻
[delta, nneigh] = computeMinDistance(distances, rho);
delta(isinf(delta)) = max(distances(:));

% 选择聚类中心
centerIdxs = selectClusterCenters(rho, delta, K);

% 分配聚类标签
clusterLabels = assignClusterLabels(nneigh, centerIdxs);

% Output Settings
% =========================================================
% clusterLabels = process_labels(clusterLabels);
if nargout == 2
    varargout = {clusterLabels, centerIdxs};
elseif nargout == 4
    varargout = {clusterLabels, centerIdxs, rho, delta};
else
    varargout = {clusterLabels};
end

end

function dc = computeNeighborhoodCutoff(distances, percent)
    % 计算邻域截断距离
    n = size(distances, 1); idx = triu(true(n,n),1);
    dc = prctile(distances(idx), percent);
end

function rho = computeLocalDensity(distances, dc, kernel)
    % 计算局部密度
    switch kernel
        case 'gaussian'         % Gaussian kernel
            rho = sum(exp(-(distances./dc).^2), 2) - 1;
        case 'cutoff'           % Cutoff kernel
            rho = sum(distances < dc, 2);
    end
end

function [minDistance, nearestNeighbor] = computeMinDistance(distances, rho)
    % 计算最小距离和最近邻点
    [~, sortedIndices] = sort(rho, 'descend');
    numPoints = size(distances, 1);
    minDistance = Inf(numPoints, 1);
    nearestNeighbor = zeros(numPoints, 1);
    
    maxRhoIdx = sortedIndices(1);
    minDistance(maxRhoIdx) = max(distances(:));
    dists = distances(maxRhoIdx, :); dists(maxRhoIdx) = Inf;
    [~, nearestIdx] = min(dists);
    nearestNeighbor(maxRhoIdx) = nearestIdx;
    
    for i = 2:numPoints
        for j = 1:i-1
            if distances(sortedIndices(i), sortedIndices(j)) < minDistance(sortedIndices(i))
                minDistance(sortedIndices(i)) = distances(sortedIndices(i), sortedIndices(j));
                nearestNeighbor(sortedIndices(i)) = sortedIndices(j);
            end
        end
    end
end

function clusterCenters = selectClusterCenters(rho, delta, K)
    % 选择聚类中心
    if K > 0
        [~, sortedIndices] = sort(rho .* delta, 'descend');
        clusterCenters = sortedIndices(1:K); % K为聚类数量
    else
        rho_thresh = prctile(rho, 70); delta_thresh = prctile(delta, 70);
        clusterCenters = find((rho>rho_thresh) & (delta>delta_thresh));
        % clusterIndices = 1:length(rho);
    end
end

% function clusterLabels = assignClusterLabels(nearestNeighbor, clusterCenters)
%     % 分配聚类标签
%     numPoints = length(nearestNeighbor);
%     clusterLabels = zeros(numPoints, 1);
%     
%     for i = 1:length(clusterCenters)
%         clusterLabels(clusterCenters(i)) = i;
%     end
%     
%     while nnz(clusterLabels) < numPoints
%         unLabelIdx = find(clusterLabels~=0);
%         disp(unLabelIdx); disp(size(nearestNeighbor));
%         for i = 1:length(unLabelIdx)
%             idx = unLabelIdx(i); disp([i idx nearestNeighbor(idx)])
%             neighLabel = clusterLabels(nearestNeighbor(idx));
%             if neighLabel > 0
%                 clusterLabels(idx) = neighLabel;
%             end
%         end
%     end
% end

function clusterLabels = assignClusterLabels(nearestNeighbor, clusterCenters)
    % 分配聚类标签
    numPoints = length(nearestNeighbor);
    clusterLabels = zeros(numPoints, 1);
    
    for i = 1:numPoints
        if ismember(i, clusterCenters)
            clusterLabels(i) = find(clusterCenters == i);
        else
            clusterLabels(i) = clusterLabels(nearestNeighbor(i));
        end
    end
end
