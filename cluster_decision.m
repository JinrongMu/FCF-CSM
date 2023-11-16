function varargout = cluster_decision(img_uint8, lab_gt, lab_cluster, varargin)
% Clustering Label Decision
%
% Syntax
% =================
% [gtLabels, tarLabels] = cluster_decision(img_uint8, lab_gt, lab_cluster, ...
%     'method', 'CDIV2', 'debug_mode', 2);
%
% Input Arguments
% =================
% img_uint8 	RGB Image, uint8, [0, 255].
% lab_gt    	Ground truth, uint8, {0, 255}.
% lab_cluster   Clustering label image, [1, N].
% method        Clustering decision method, {'CDIV1','CDIV2','HDIV1','HDIV2','HueV1','HueV2','LABV1'}
% debug_mode    Print debug information, 0: Silent, 1: Call information, 2: Call details.

    % Parameter Initialization
    % =========================================================    
    arg = inputParser; fun_name = 'cluster_decision'; 
    addParameter(arg,'method','CDIV2');
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
    img_float = double(img_uint8) / 255.;
    img_gt = double(lab_gt) / 255.;
    [lab_color, cmap] = imregional(img_float, lab_cluster); 
    [gtLabels, gtProbs, gtThresh] = decisionProbThresh(img_gt, ...
        lab_cluster, 0.5 , 'scale', 0.9, 'debug_mode', debug_mode-1);
    
    method = arg.Results.method;
    if strcmpi(method, 'CDIV1')
        probMap = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
        [tarLabels, tarProbs, tarThresh] = decisionProbThresh(probMap, ...
            lab_cluster, 0.8 , 'scale', 0.9, 'debug_mode', debug_mode-1);
        
    elseif strcmpi(method, 'CDIV2')
        optMap = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
        deltaE_tol = 10;
        [tarLabels, meanValues, optLabel, deltaE] = decisionOptLabel(optMap, ...
            lab_cluster,cmap,'deltaE_tol',deltaE_tol,'debug_mode', debug_mode-1);
        
    elseif strcmpi(method, 'HueV1')
        probMap = get_features(img_uint8, 'Hue', 'HueNorm');
        [tarLabels, tarProbs, tarThresh] = decisionProbThresh(probMap, ...
            lab_cluster, 0.8 , 'scale', 0.9, 'debug_mode', debug_mode-1);
        
    elseif strcmpi(method, 'HueV2')
        img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
        img_hue = get_features(img_uint8, 'Hue', 'HueNorm');
        probMap = 0.75 * img_cdi + 0.25 * img_hue;
        [tarLabels, tarProbs, tarThresh] = decisionProbThresh(probMap, ...
            lab_cluster, 0.8 , 'scale', 0.8, 'debug_mode', debug_mode-1);
        
    elseif strcmpi(method, 'HDIV1')
        probMap = get_features(img_uint8, 'Hue', 'HueProb');
        [tarLabels, tarProbs, tarThresh] = decisionProbThresh(probMap, ...
            lab_cluster, 0.8 , 'scale', 0.9, 'debug_mode', debug_mode-1);
        
    elseif strcmpi(method, 'HDIV2')
        img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
        img_hdi = get_features(img_uint8, 'Hue', 'HueProb');
        probMap = 0.75 * img_cdi + 0.25 * img_hdi;
        [tarLabels, tarProbs, tarThresh] = decisionProbThresh(probMap, ...
            lab_cluster, 0.8 , 'scale', 0.9, 'debug_mode', debug_mode-1);
        
    elseif strcmpi(method, 'LABV1')
        img_lab = rgb2lab(img_float); 
        optMap = img_lab(:,:,2); deltaE_tol = 10;
        [tarLabels, meanValues, optLabel, deltaE] = decisionOptLabel(optMap, ...
            lab_cluster,cmap,'deltaE_tol',deltaE_tol,'debug_mode', debug_mode-1);
    end
    
    
    % Output Settings
    % =========================================================
    if nargout == 2
        varargout = {gtLabels, tarLabels};
    else
        varargout = {gtLabels};
    end
    
    % Debug Information
    % =========================================================
    if debug_mode == 2
        
        fprintf('TargetLabels: %s.\n', mat2str(tarLabels));
        
        figure
        subplot(1,2,1), imshow(img_uint8), title('Image')
        subplot(1,2,2), imshow(lab_color), title('Cluster Lable')
        
        groups = {'LABV1', 'CDIV2'};  % Optimal Label Method
        if any(strcmpi(method, groups))
            figure
            ax1 = subplot(1,3,1); xy_labels = {'Label', 'Lesion Ratios'};
            draw_decision_graph(gtProbs, gtThresh, cmap, ax1, ...
                'direction','none','xy_labels',xy_labels,'titles','GT');     
            ax2 = subplot(1,3,2); xy_labels = {'Label', 'Decision Indexs'};
            optThresh = 0; str_title = sprintf("Optical Label: %d", optLabel);
            draw_decision_graph(meanValues, optThresh, cmap, ax2, ...
                'direction','none','xy_labels',xy_labels,'titles',str_title); 
            ax3 = subplot(1,3,3); xy_labels = {'Label', 'deltaE'};
            optValues = deltaE + 5;
            optThresh = optValues(optLabel) + deltaE_tol; 
            str_title = sprintf("Target Labels: %s", mat2str(tarLabels'));
            draw_decision_graph(optValues, optThresh, cmap, ax3, ...
                'direction','ascend','xy_labels',xy_labels,'titles', str_title);
            
            suptitle('Decision Graph');

        else   % Probability Threshold Method
            figure
            ax1 = subplot(1,2,1); xy_labels = {'Label', 'Lesion Ratios'};
            draw_decision_graph(gtProbs, gtThresh, cmap, ax1, ...
                'direction','none','xy_labels',xy_labels,'titles','GT');
            ax2 = subplot(1,2,2); xy_labels = {'Label', 'Decision Indexs'};
            draw_decision_graph(tarProbs, tarThresh, cmap, ax2, ...
                'direction','none','xy_labels',xy_labels,'titles',method);

        end
        
    end

end

function varargout = decisionProbThresh(probMap, labelImage, thresh, varargin)
% Probability Threshold Method
%
% Syntax
% =================
% [targetLabels, meanProbs, thresh] = decisionProbThresh(probMap, ...
%     labelImage, thresh,'scale',0.9,'debug_mode',0);

    % Parameter Initialization
    % =========================================================    
    arg = inputParser; fun_name = 'decisionProbThresh'; 
    addParameter(arg,'scale',0.9); 
    addParameter(arg,'debug_mode',0);
    parse(arg,varargin{:});
    
    % Method Implementation
    % =========================================================
    [~, meanProbs] = imregional(probMap, labelImage); 
    targetLabels = find(meanProbs > thresh);
    if isempty(targetLabels)
        thresh = max(meanProbs(:)) * arg.Results.scale;
        targetLabels = find(meanProbs > thresh);
    end
    
    % Output Settings
    % =========================================================
    if nargout == 3
        varargout = {targetLabels, meanProbs, thresh};
    elseif nargout == 2
        varargout = {targetLabels, meanProbs};
    else
        varargout = {targetLabels};
    end
    
    % Debug Information
    % =========================================================
    if arg.Results.debug_mode == 1
        fprintf('\nCall functions:\t%s\n', fun_name)      
    elseif arg.Results.debug_mode == 2
        fprintf('\nCall functions:\t%s\n', fun_name)
        fprintf('----------------------------------------\n');
        fprintf('targetLabels: %s.\n', mat2str(targetLabels));
        fprintf('probValues: %s.\n', mat2str(meanProbs));
        fprintf('thresh: %.2f.\n', thresh);   
    end
   
end

function varargout = decisionOptLabel(optMap, labelImage, cmap, varargin)
% Optimal Category Method deltaE_tol = 10
%
% Syntax
% =================
% [targetLabels, meanValues, optLabel, deltaE] = decisionOptLabel(optMap, ...
%     labelImage, cmap, 'deltaE_tol', 10, 'debug_mode', 2);

    % Parameter Initialization
    % =========================================================    
    arg = inputParser; fun_name = 'decisionProbThresh'; 
    addParameter(arg,'deltaE_tol',10); 
    addParameter(arg,'debug_mode',0);
    parse(arg,varargin{:});
    
    % Method Implementation
    % =========================================================
    [~, meanValues] = imregional(optMap, labelImage); 
    [~, optLabel] = max(meanValues);
    cmap_lab = rgb2lab(cmap);
    deltaE = pdist2(cmap_lab, cmap_lab(optLabel,:));
    targetLabels = find(deltaE < arg.Results.deltaE_tol);

    
    % Output Settings
    % =========================================================
    if nargout == 4
        varargout = {targetLabels, meanValues, optLabel, deltaE};
    elseif nargout == 3
        varargout = {targetLabels, meanValues, optLabel};
    elseif nargout == 2
        varargout = {targetLabels, meanValues};
    else
        varargout = {targetLabels};
    end
    
    % Debug Information
    % =========================================================
    if arg.Results.debug_mode == 1
        fprintf('\nCall functions:\t%s\n', fun_name)      
    elseif arg.Results.debug_mode == 2
        fprintf('\nCall functions:\t%s\n', fun_name)
        fprintf('----------------------------------------\n');
        fprintf('targetLabels: %s.\n', mat2str(targetLabels));
        fprintf('OptimalLable: %d.\n', optLabel);   
    end
   
end

function ax = draw_decision_graph(values, thresh, cmap, ax, varargin)
% Drawing Decision Graph
%
% Syntax
% =================
% draw_decision_graph(values, thresh, cmap, ax, 'direction', 'none', ...
%     'xy_labels',{'Label', 'Values'}, 'titles', 'Decision Graph'); 
%
    % Parameter Initialization
    % =========================================================
    if isempty(ax), ax = gca; end
    % ax = gcf; set(ax, 'Name', 'Decision Graph');
    if isempty(cmap), cmap = get_colors(length(values)); end
    
    arg = inputParser; fun_name = 'decisionProbThresh'; 
    addParameter(arg,'direction','none'); % {'ascend', 'descend', 'none'}
    addParameter(arg,'xy_labels',{'Label', 'Values'});
    addParameter(arg,'titles','Decision Graph');
    addParameter(arg,'debug_mode',0);
    parse(arg,varargin{:});
    
    if arg.Results.debug_mode > 0
        fprintf('\nCall functions:\t%s\n', fun_name)   
    end
    
    % Method Implementation
    % =========================================================
    direction = arg.Results.direction; 
    if strcmpi(direction, 'ascend') || strcmpi(direction, 'descend') 
        [~, idxs] = sort(values, direction);
        values = values(idxs); cmap = cmap(idxs, :);
        is_sorted = true;
    else
        is_sorted = false;
    end

    b = bar(ax, values); b.FaceColor = 'flat'; 
    for k=1:length(values), b.CData(k,:) = cmap(k,:); end
    
    if is_sorted
        xticks(1:length(values));
        xticklabels(idxs)
    end
    
    xy_labels = arg.Results.xy_labels; 
    [xstr, ystr] = deal(xy_labels{:});
    xlabel(xstr); ylabel(ystr); title(arg.Results.titles)
    
    if thresh ~= 0
        hold on;
        x = [0, length(values)+1]; y = [thresh, thresh];
        line(ax, x, y, 'Color', [0 0 0], 'LineWidth', 0.5, 'LineStyle', '--')
        hold off
    end
    
    % Set graphics style: IEEE Journal Style
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1.5);
    set(gcf, 'Color', 'w');
end

