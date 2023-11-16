function varargout = cluster_DPC(data, K, varargin)
% cluster_DPC - Ensemble of density peak clustering algorithms.
%
% Syntax
% =================
% L = cluster_DPC(data);                % K = 0
% [L, idxs] = cluster_DPC(data, K);
% [L, idxs, rho, delta] = cluster_DPC(data, K, 'percent', 4, ...
%     'kernel', 'gaussian', 'version', 'v1', 'debug_mode',1)
%
% Input Arguments
% =================
% data          Data to be clustered, NxM matrix.
% K             Number of clusters, K=0 for automatic clustering.
% 
% percent       Percentage used to select dc distance, default is 4.
% kernel        The type of kernel to calculate the density rho, ['gaussian','cutoff'].
% version       Clustering algorithm version.
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
%
% Parameter Initialization    
% =========================================================
% addpath('./Modules/Clustering/')
% addpath('./Functions/')      % Required for the function 'process_labels'.

if ~exist('K','var'),  K = 0; end

arg = inputParser; fun_name = 'cluster_DPC';
addParameter(arg,'percent',4);
addParameter(arg,'kernel','gaussian');
addParameter(arg,'num_neigh',5);  
addParameter(arg,'version','v1');
addParameter(arg,'debug_mode',0); 
parse(arg,varargin{:});

debug_mode = arg.Results.debug_mode;
if debug_mode == 1
    fprintf('\nCall functions:\t%s\n', fun_name);
elseif debug_mode == 2     
    fprintf('\nCall functions:\t%s\n', fun_name);
    fprintf('----------------------------------------');
    fprintf('\nDefault Parameters:\n'); disp(arg.Results);
end

% Method Implementation 
% =========================================================
percent = arg.Results.percent;
kernel = arg.Results.kernel;

if strcmpi(arg.Results.version, 'v1')
    [clusterLabels, centerIdxs, rho, delta] = DPC_v1(data, ...
        'percent', percent, 'kernel', kernel, 'debug_mode', debug_mode-1);
elseif strcmpi(arg.Results.version, 'v2')
    [clusterLabels, centerIdxs, rho, delta] = DPC_v2(data, K, ...
        'percent', percent, 'kernel', kernel, 'debug_mode', debug_mode-1);
elseif strcmpi(arg.Results.version, 'v3')
    [clusterLabels, centerIdxs, rho, delta] = DPC_v3(data, K, ...
        'percent', percent, 'kernel', kernel, 'debug_mode', debug_mode-1);   
else % 'v4'
    num_neigh = arg.Results.num_neigh;   
    [clusterLabels, centerIdxs, rho, delta] = DPC_v4(data, 'percent', percent, ...
        'kernel', kernel, 'num_neigh', num_neigh,'debug_mode', debug_mode-1);
end

clusterLabels = process_labels(clusterLabels);

% Output Settings
% =========================================================
if nargout == 4
    varargout = {clusterLabels, centerIdxs, rho, delta};
elseif nargout == 3
    K = length(centerIdxs);
    varargout = {clusterLabels, centerIdxs, K};
elseif nargout == 2
    varargout = {clusterLabels, centerIdxs};
else
    varargout = {clusterLabels};  
end

% Debug Information
% =========================================================
if debug_mode == 2
    selectIdx = centerIdxs; gamma = rho .* delta;
    selectLabel = double(clusterLabels(selectIdx)); 
    selectData = [selectIdx, rho(selectIdx), delta(selectIdx), ...
        gamma(selectIdx), selectLabel];
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
    hold off
    
    subplot(1,2,2),
    plot(1:length(rho),gamma,'o','MarkerSize',4,'MarkerFaceColor','k','MarkerEdgeColor','k');
    hold on; 
    plot(selectIdx,gamma(selectIdx),'o','MarkerSize',10,...
        'MarkerFaceColor','red','MarkerEdgeColor','green');
    title ('Decision Graph'); xlabel ('id'); ylabel ('\gamma');
    hold off
    
end

end