function [varargout] = DPC_v4(data, varargin)
% DPC_v4 - Density Peak Clustering (DPC) V4. (Error)
% 
% Version V4: Find data points with high local density and large gradient changes as cluster centers.
% 
% Syntax
% =================
% clusterLabels = DPC_v4(data)
% [clusterLabels, centerIdxs] = DPC_v4(data)
% [clusterLabels, centerIdxs, rho, delta] = DPC_v4(data,'percent',2,'kernel','cutoff','num_neigh',5,'debug_mode',1)
% 
% Input Arguments
% =================
% data          Data to be clustered, NxM matrix.
% percent       Percentage used to select dc distance, default is 4.
% kernel        The type of kernel to calculate the density rho, ['gaussian','cutoff'].
% num_neigh     Minimum Neighborhood Points
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
arg = inputParser; fun_name = 'DPC_v4';                         
addParameter(arg,'percent',2);
addParameter(arg,'kernel','cutoff');  
addParameter(arg,'num_neigh',5);  
addParameter(arg,'debug_mode',0); 
parse(arg,varargin{:});

% Method Implementation 
% =========================================================

% Calculate the distance matrix
distMatrix = pdist2(data, data);

% Setting parameters
% epsilon = 0.3;  % density threshold dc
% minPts = 5;     % The number of minimum neighbor points
epsilon = prctile(distMatrix(:), arg.Results.percent);
minPts = arg.Results.num_neigh;   

% Calculate density
kernel = lower(arg.Results.kernel); 
switch kernel    
    case 'gaussian'         % Gaussian kernel
        density = sum(exp(-(distMatrix./epsilon).^2), 2) - 1;      
    case 'cutoff'           % Cutoff kernel
        density = sum(distMatrix < epsilon, 2);
end

% Calculate the gradient of the local density
gradient = zeros(size(data, 1), 1);
for i = 1:size(data, 1)
    neighborIdx = knnsearch(data, data(i, :), 'K', minPts+1);
    neighborDist = distMatrix(i, neighborIdx(2:end));
    gradient(i) = sum(density(neighborIdx(2:end)) > density(i) & (neighborDist' < epsilon));
end

% Find the density peak as the cluster center
clusterCenterIdx = find(gradient==max(gradient));

% Assign data points to cluster centers
clusterIdx = zeros(size(data, 1), 1);
for i = 1:length(clusterCenterIdx)
    clusterIdx(clusterCenterIdx(i)) = i;
end
clusterId = 1;
for i = 1:size(data, 1)
    if i == clusterCenterIdx
        continue;
    end
    neighborIdx = knnsearch(data, data(i, :), 'K', minPts+1);
    neighborDist = distMatrix(i, neighborIdx(2:end));
    [~, nearestCenterIdx] = min(neighborDist);
    if density(i) > density(neighborIdx(nearestCenterIdx+1)) || density(i) == density(neighborIdx(nearestCenterIdx+1)) && i < neighborIdx(nearestCenterIdx+1)
        clusterId = clusterId + 1;
    end
    clusterIdx(i) = clusterId;
end

% Output Settings
% =========================================================
if nargout == 2
    varargout = {clusterIdx, clusterCenterIdx};
elseif nargout == 4
    varargout = {clusterIdx, clusterCenterIdx, density, gradient};
else
    varargout = {clusterIdx};  
end

% Debug Information
% =========================================================
if arg.Results.debug_mode == 1
    fprintf('\nCall functions:\t%s\n', fun_name)  
elseif arg.Results.debug_mode == 2     
    fprintf('\nCall functions:\t%s\n', fun_name)
    fprintf('----------------------------------------');
    fprintf('\nDefault Parameters:\n'); disp(arg.Results);    

    clusterLabels = clusterIdx;
    selectIdx = clusterCenterIdx; 
    rho = density; delta = gradient;
    selectData = [selectIdx, rho(selectIdx), delta(selectIdx), clusterLabels(selectIdx)];
    varNames = {'index','rho','delta','label'};
    debug_info = array2table(selectData,'VariableNames',varNames);
    fprintf('\nClustering Center Information:\n'); disp(debug_info)
    
    figure;
    plot(rho,delta,'o','MarkerSize',4,'MarkerFaceColor','k','MarkerEdgeColor','k');
    hold on;     
    plot(rho(selectIdx),delta(selectIdx),'o','MarkerSize',10,...
        'MarkerFaceColor','red','MarkerEdgeColor','green');
    title ('Decision Graph'); xlabel ('\rho'); ylabel ('\delta');  
    
end

end
