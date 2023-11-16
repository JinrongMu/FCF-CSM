function [varargout] = DPC_v3(data, K, varargin)
% DPC_v3 - Density Peak Clustering (DPC) V3.
% 
% Version V3: Select the cluster center according to the number of clusters, 
% the number of cluster centers is less than the number of clusters, 
% support manual selection of cluster centers.
% 
% Syntax
% =================
% clusterLabels = DPC_v3(data)
% [clusterLabels, centerIdxs] = DPC_v3(data)
% [clusterLabels, centerIdxs, rho, delta] = DPC_v3(data,'percent',4,'kernel','gaussian','debug_mode',1)
%
% Input Arguments
% =================
% data          Data to be clustered, NxM matrix.
% K             Number of clusters, manually selected when it is 0.
%
% percent       Percentage used to select dc distance, default is 4.
% kernel        The type of kernel to calculate the density rho, ['gaussian','cutoff'].
% debug_mode    Control debug information, 0: Silent, 1: Call information, 2: Call details.	
%
% Output Arguments:
% ================
% cluster_labels	The labels of the clusters. Lable equals to 0 means it's in the halo region
% center_idxs       Index of cluster center.

% Parameter Initialization
% =========================================================
arg = inputParser; fun_name = 'DPC_v3';                       
addParameter(arg,'percent',2);
addParameter(arg,'kernel','gaussian'); 
addParameter(arg,'debug_mode',0); 
addParameter(arg,'is_debug',0); 
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

%% Density Peak Clustering Procedure.
fprintf('\nStart Density Peak Clustering\n');
fprintf('----------------------------------\n');

dist = squareform(pdist(data));

% Estimate dc
disp('Estimating dc...');
percent = arg.Results.percent;
N = size(dist,1); position = round(N*(N-1)*percent/100);
tri_u = triu(dist,1); sda = sort(tri_u(tri_u~=0));
dc = sda(position); 

% Compute rho(density)
fprintf('Computing Rho with gaussian kernel of radius: %12.6f\n', dc);
kernel = lower(arg.Results.kernel);
switch kernel    
    case 'gaussian'
        rho = sum(exp(-(dist./dc).^2),2)-1;      
    case 'cutoff' 
        rho = sum((dist-dc)<0, 2);
end
[~, ordrho] = sort(rho, 'descend');

% Compute delta
disp('Computing delta...');
delta = zeros(size(rho)); nneigh = zeros(size(rho));
delta(ordrho(1)) = -1; nneigh(ordrho(1)) = 0;
for i = 2:size(dist,1)
    range = ordrho(1:i-1);
    [delta(ordrho(i)), tmp_idx] = min(dist(ordrho(i),range));
    nneigh(ordrho(i)) = range(tmp_idx); 
end
delta(ordrho(1)) = max(delta(:));

%% Decision graph, choose min rho and min delta
if K == 0
    figure;
    plot(rho(:),delta(:),'o','MarkerSize',4,'MarkerFaceColor','k','MarkerEdgeColor','k');
    hold on; title ('Decision Graph'); xlabel ('\rho'); ylabel ('\delta');
    rect = getrect(gcf); % get user decision
    rho_min = rect(1); delta_min = rect(2); bad_idx = find(rho < rho_min);
    rho_range = rect(3); delta_range = rect(4);
    target_rho = rho>rho_min & rho<rho_min+rho_range;
    target_delta = delta>delta_min & delta<delta_min+delta_range;
    select_index = find(target_rho & target_delta);
    select_rho = rho(select_index); select_delta = delta(select_index);
else
    [~, orddelta] = sort(delta, 'descend');
    select_index = orddelta(1:K);
    select_rho = rho(select_index); select_delta = delta(select_index);
    rho_min = min(select_rho); delta_min = min(select_delta);
    bad_idx = find(rho < rho_min);    
end

if arg.Results.debug_mode == 2
    figure;
    plot(rho(:),delta(:),'o','MarkerSize',4,'MarkerFaceColor','k','MarkerEdgeColor','k');
    hold on; title ('Decision Graph'); xlabel ('\rho'); ylabel ('\delta');
    plot(rho(select_index),delta(select_index),'o','MarkerSize',10,...
        'MarkerFaceColor','red','MarkerEdgeColor','green');
    hold off;
end

% Find cluster centers
disp('Finding cluster centers...');
center_idxs = (find(delta>delta_min));
% delete centers whose rho is samller than rho_min
for i = 1:length(center_idxs)
    if ~isempty(find(center_idxs(i)==bad_idx, 1))
        center_idxs(i) = -1;
    end
end
center_idxs(center_idxs==-1) = [];
disp([num2str(length(center_idxs)), ' cluster centers found...']);


%% Assignment

% raw assignment
disp('Assigning data-points into clusters...');
cluster_lables = -1*ones(size(dist,1),1);
for i = 1:length(center_idxs)
    cluster_lables(center_idxs(i)) = i;
end
for i=1:length(cluster_lables)
    if (cluster_lables(ordrho(i))==-1)
        cluster_lables(ordrho(i)) = cluster_lables(nneigh(ordrho(i)));
    end
end
raw_cluster_lables = cluster_lables;

% find and cut off halo region
disp('Cut off halo regions...');
for i = 1:length(center_idxs)
    tmp_idx = find(raw_cluster_lables==i);
    tmp_dist = dist(tmp_idx,:);
    tmp_dist(:,tmp_idx) = max(dist(:));
    tmp_rho = rho(tmp_idx);
    tmp_lables = raw_cluster_lables(tmp_idx);
    tmp_border = find(sum(tmp_dist<dc,2)>0);
    if ~isempty(tmp_border)
        rho_b = max(tmp_rho(tmp_border));
        halo_idx = rho(tmp_idx) < rho_b;
        tmp_lables(halo_idx) = 0;
        % lable equals to 0 means it's in the halo region
        cluster_lables(tmp_idx) = tmp_lables;
    end
end
cluster_lables = cluster_lables + 1; center_idxs = select_index;

% Output Settings
% =========================================================
% cluster_lables = process_labels(cluster_lables);
if nargout == 4
    varargout = {cluster_lables, center_idxs, rho, delta};
elseif nargout == 2
    varargout = {cluster_lables, center_idxs};
else
    varargout = {cluster_lables};
end

% Debug Information
% =========================================================
if arg.Results.debug_mode == 2
    center_data = [select_index,select_rho,select_delta,cluster_lables(select_index)];
    T = array2table(center_data, 'VariableNames',{'index','rho','delta','label'});
    disp(T);
end

end