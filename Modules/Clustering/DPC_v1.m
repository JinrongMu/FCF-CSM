function [varargout] = DPC_v1(data, varargin)
% DPC_v1 - Density Peak Clustering (DPC) V1 .
% 
% Version V1: Automatically selects the cluster center and is equal to the
% number of clusters, does not support the specified number of clusters.
% 
% Syntax
% =================
% clusterLabels = DPC_v1(data)
% [clusterLabels, centerIdxs] = DPC_v1(data)
% [clusterLabels, centerIdxs, rho, delta] = DPC_v1(data,'percent',4,'kernel','gaussian','debug_mode',1)
% 
% Input Arguments
% =================
% data          Data to be clustered, NxM matrix.
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
arg = inputParser; fun_name = 'DPC_v1';                         
addParameter(arg,'percent',4);
addParameter(arg,'kernel','gaussian');  
addParameter(arg,'debug_mode',0); 
parse(arg,varargin{:});

% Method Implementation 
% =========================================================
[n, ~] = size(data);

% Calculate pairwise distances
D = squareform(pdist(data));

% Calculate the neighborhood cut-off distance
percent = arg.Results.percent;
idx = triu(true(n,n),1); 
dc = prctile(D(idx), percent);
% triu_D = triu(D,1); 
% sort_D = sort(triu_D(triu_D~=0)); % sort_D = sort(D(idx)); 
% position = round(n*(n-1)*percent/100); % position = round(nnz(idx)*percent/100); 
% dc = sort_D(position); % or dc = prctile(sort_D, percent);

% Calculate density and sort in descending order
kernel = lower(arg.Results.kernel); 
switch kernel    
    case 'gaussian'         % Gaussian kernel
        rho = sum(exp(-(D./dc).^2), 2) - 1;      
    case 'cutoff'           % Cutoff kernel
        rho = sum(D < dc, 2);
end
[~, rhoIdx] = sort(rho, 'descend');

% Calculate delta for each point
delta = inf(n, 1); radius = dc; 
for i = 1:n
    for j = i+1:n
        if D(i,j) < radius && rho(i) < rho(j)
            delta(i) = min(delta(i), D(i,j));
            delta(j) = min(delta(j), D(i,j));
        end
    end
end
% delta = zeros(n, 1); delta(rhoIdx(1)) = max(D(:)); 
% for j = 2:n
%     jIdx = rho > rho(rhoIdx(j));
%     delta(rhoIdx(j)) = min(D(rhoIdx(j),jIdx));
% end

% Assign cluster indices
clusterLabels = zeros(n, 1); 
currentLabel = 0;
for i = 1:n
    if delta(rhoIdx(i)) == inf
        continue;
    end
    if clusterLabels(rhoIdx(i)) == 0
        currentLabel = currentLabel + 1;
        clusterLabels(rhoIdx(i)) = currentLabel;
    end
    neighborIdx = find(D(rhoIdx(i), :) < radius);
    for j = 1:length(neighborIdx)
        if clusterLabels(neighborIdx(j)) == 0
            clusterLabels(neighborIdx(j)) = clusterLabels(rhoIdx(i));
        elseif delta(neighborIdx(j)) > delta(rhoIdx(i))
            clusterLabels(neighborIdx(j)) = clusterLabels(rhoIdx(i));
        end
    end
end

% Assign cluster centers
labels = unique(clusterLabels);
centerIdxs = zeros(length(labels), 1); 
for i = 1:length(labels)
    labelCnt = labels(i); labelRho = zeros(n, 1);
    labelIdx = find(clusterLabels==labelCnt);
    labelRho(labelIdx) = rho(labelIdx);
    [~, centerIdx] = max(labelRho); 
    centerIdxs(i) = centerIdx;   
end

% Output Settings
% =========================================================
delta(isinf(delta)) = max(D(:));
% clusterLabels = process_labels(clusterLabels);
if nargout == 2
    varargout = {clusterLabels, centerIdxs};
elseif nargout == 4
    varargout = {clusterLabels, centerIdxs, rho, delta};
else
    varargout = {clusterLabels};  
end

% Debug Information
% =========================================================
if arg.Results.debug_mode == 1
    fprintf('\nCall functions:\t%s\n', fun_name)  
elseif arg.Results.debug_mode == 2     
    fprintf('\nCall functions:\t%s\n', fun_name)
    fprintf('----------------------------------------');
    fprintf('\nDefault Parameters:\n'); disp(arg.Results);    

    selectIdx = centerIdxs; gamma = rho .* delta;
    selectData = [selectIdx, rho(selectIdx), delta(selectIdx), ...
        gamma(selectIdx), clusterLabels(selectIdx)];
    varNames = {'index','rho','delta','gamma','label'};
    debug_info = array2table(selectData,'VariableNames',varNames);
    fprintf('\nClustering Center Information:\n'); disp(debug_info)
    
    figure;
    subplot(1,2,1),
    plot(rho,delta,'o','MarkerSize',4,'MarkerFaceColor','k','MarkerEdgeColor','k');
    hold on;     
    plot(rho(selectIdx),delta(selectIdx),'o','MarkerSize',10,...
        'MarkerFaceColor','red','MarkerEdgeColor','green');
    title ('Decision Graph'); xlabel ('\rho'); ylabel ('\delta');  
    
    subplot(1,2,2),
    plot(1:length(rho),gamma,'o','MarkerSize',4,'MarkerFaceColor','k','MarkerEdgeColor','k');
    hold on; 
    plot(selectIdx,gamma(selectIdx),'o','MarkerSize',10,...
        'MarkerFaceColor','red','MarkerEdgeColor','green');
    title ('Decision Graph'); xlabel ('id'); ylabel ('\gamma');
    
end

end
