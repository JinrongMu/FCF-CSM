% addpath('./Functions/')  % Required for the function 'imregional' and 'process_labels'.
function varargout = pixel_blocks_merge(img_uint8, initLabel, varargin)
% Initial merging of pixel blocks.

    % Parameter Initialization
    % =========================================================
    arg = inputParser; fun_name = 'pixel_blocks_merge';
    addParameter(arg,'merge_method','CDIV1');
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
    merge_method = arg.Results.merge_method;
    if strcmpi(merge_method, 'ColorV1')
        img_lab = rgb2lab(double(img_uint8)/255.); deltaE_tol = 2.3;
        [mergeLabel, mergeN] = mergeColorDifference(img_lab,initLabel,deltaE_tol);
        
    elseif strcmpi(merge_method, 'ColorV2')
        img_rgb = double(img_uint8)/255.; deltaE_tol = 0.05;
        [mergeLabel, mergeN] = mergeColorDifference(img_rgb,initLabel,deltaE_tol);
        
    elseif strcmpi(merge_method, 'ColorV3')
        img_lab = rgb2lab(double(img_uint8)/255.); deltaE_tol = 1;
        img_prob = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
        lab_cdi = cat(3, img_lab, img_prob);
        [mergeLabel, mergeN] = mergeColorGuide(lab_cdi, initLabel, deltaE_tol);
        
    elseif strcmpi(merge_method, 'ColorV4')
        img_lab = rgb2lab(double(img_uint8)/255.); deltaE_tol = 1;
        mu = 18.2676; sigma_ = 5.6107; max_prob = 0.0710;
        x = img_lab(:,:,2);
        y = exp(-(x-mu).^2/(2*sigma_^2))/(sqrt(2*pi)*sigma_);
        lab_prob = cat(3, img_lab, y/max_prob);
        [mergeLabel, mergeN] = mergeColorGuide(lab_prob, initLabel, deltaE_tol);
        
    elseif strcmpi(merge_method, 'ColorV5')
        img_lab = rgb2lab(double(img_uint8)/255.); 
        [mergeLabel, mergeN] = mergeProbGuide(img_lab(:,:,2), initLabel, ...
            'thresh_policy', 'Static', 'merge_policy', 'Unilateral', ...
            'lowerT', 0, 'upperT', 60);
        
    elseif strcmpi(merge_method, 'ColorV6')
        img_lab = rgb2lab(double(img_uint8)/255.); deltaE_tol = 2.3;
        [mergeLabel, mergeN] = mergeColorDifference(img_lab(:,:,2:3),initLabel,deltaE_tol);
        
    elseif strcmpi(merge_method, 'CDIV1')
        img_prob = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
        [mergeLabel, mergeN] = mergeProbGuide(img_prob, initLabel, ...
            'thresh_policy', 'Dynamic', 'merge_policy', 'Unilateral', ...
            'lowerT', 0.4, 'upperT', 0.9);
        
    elseif strcmpi(merge_method,'CDIV2')
        img_prob = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
        [mergeLabel, mergeN] = mergeProbGuide(img_prob, initLabel, ...
            'thresh_policy', 'Dynamic', 'merge_policy', 'Bilateral', ...
            'lowerT', 0.5, 'upperT', 0.9);
        
    elseif strcmpi(merge_method,'HDIV1')
        img_prob = get_features(img_uint8, 'Hue', 'HueProb');
        [mergeLabel, mergeN] = mergeProbGuide(img_prob, initLabel, ...
            'thresh_policy', 'Dynamic', 'merge_policy', 'Unilateral', ...
            'lowerT', 0.6, 'upperT', 0.95);
        
    elseif strcmpi(merge_method,'HDIV2')
        img_prob = get_features(img_uint8, 'Hue', 'HueProb');
        [mergeLabel, mergeN] = mergeProbGuide(img_prob, initLabel, ...
            'thresh_policy', 'Dynamic', 'merge_policy', 'Bilateral', ...
            'lowerT', 0.6, 'upperT', 0.95);
        
    else % merge_method: 'None' or 'default'
        [mergeLabel, mergeN] = deal(initLabel, length(unique(initLabel)));
        
    end
    
    % Output Settings
    % =========================================================
    if nargout == 2
        varargout = {mergeLabel, mergeN};
    else
        varargout = {mergeLabel};
    end
    
    % Debug Information
    % =========================================================
    if arg.Results.debug_mode == 2        
        figure
        
        [frows,fcols] = deal(1, 2);
        subplot(frows,fcols,1); 
        BW = boundarymask(initLabel); 
        imshow(imoverlay(img_uint8,BW,'cyan'),'InitialMagnification','fit')
        xstr = sprintf('N=%d', length(unique(initLabel)));
        xlabel(xstr), title('Superpixel Generation');
        
        subplot(frows,fcols,2); 
        BW = boundarymask(mergeLabel);
        imshow(imoverlay(img_uint8,BW,'cyan'),'InitialMagnification','fit')
        xlabel(sprintf('N=%d', mergeN)), title('Superpixel Merging'); 
        
        sgtitle('Coarse Clustering'); % supertitle(merge_method); % 
        
        figure
        subplot(1,2,1)
        imshow(mergeLabel, [])
        colormap('jet'), colorbar();
        xlabel(sprintf('N=%d', mergeN)), title('Merged Labels'); 
        subplot(1,2,2)
        colorLabel = imregional(double(img_uint8)/255., mergeLabel);
        imshow(colorLabel);
        xlabel(sprintf('N=%d', mergeN)), title('Color Labels');
        
        sgtitle('Coarse Clustering'); % supertitle(merge_method); % 
 
    end

end

function [mergeLabel, mergeN] = mergeColorGuide(lab_cdi, initLabel, deltaE_tol)
% Merging of labels of color images based on color difference.
% lab color difference (deltaE) + gamut dynamic threshold (gamutT)

    if nargin < 3
        deltaE_tol = 2;  % For Lab color difference
    end

    initN = length(unique(initLabel));
    [~, superFeatures] = imregional(lab_cdi, initLabel);
    
    superLabels = zeros(initN, 1); currN = 1;   
    while nnz(superLabels == 0)
        deltaE = pdist2(superFeatures(:,1:3), superFeatures(currN,1:3));
        maxGamut = max(superFeatures(:,4), superFeatures(currN,4));
        gamutT = (2 * (maxGamut - 1) .^2 + 1) * deltaE_tol; % [deltaE, 3*deltaE]
        % gamutT = (3 * (maxGamut - 1) .^2 + 1) * deltaE_tol; % [deltaE, 4*deltaE]
        % gamutT = - 2 * maxGamut + 4; % [deltaE, 4*deltaE]
        loc = [1:initN]';
        indexs = (loc >= currN) & (deltaE < gamutT);
%         curLabel = superLabels(currN);
%         if curLabel == 0
%             curLabel = max(superLabels) + 1;
%             % superLabels(currN) = curLabel;
%         end
        if superLabels(currN) == 0
            curLabel = max(superLabels) + 1;
        else
            curLabel = superLabels(currN);
        end
        superLabels(indexs) = curLabel;
        currN = currN + 1;
    end
    
    mergeLabel = superLabels(initLabel);
    mergeN = max(superLabels);
 
end

function [mergeLabel, mergeN] = mergeColorDifference(img_color, initLabel, deltaE_tol)
% Merging of labels of color images based on color difference.

    if nargin < 3
        deltaE_tol = 4;  % For Lab color difference
    end

    initN = length(unique(initLabel));
    currN = initN;  mergeN = 0;
    
    while mergeN==0 || mergeN < currN
        if mergeN == 0
            currLabel = initLabel; 
        else
            currLabel = mergeLabel; currN = mergeN;
        end
        
        [~, mean_color] = imregional(img_color, currLabel);
        dists = squareform(pdist(mean_color));
        area_no = 1:currN;
        
        for i = 1:currN-1
            for j = i+1:currN
                if area_no(j)==j && dists(i,j) < deltaE_tol
                    area_no(j) = i;
                end
            end
        end
 
        mergeLabel = area_no(currLabel);  
        mergeLabel = process_labels(mergeLabel);
        mergeN = length(unique(mergeLabel));
    end
 
end

function [mergeLabel, mergeN] = mergeProbGuide(prob_map, initLabel, varargin)
% Merging of labels of color images based on probabilistic map.
    
    arg = inputParser; 
    addParameter(arg,'thresh_policy','Dynamic');    % {'Static ', 'Dynamic'}
    addParameter(arg,'merge_policy','Unilateral');  % {'Unilateral', 'Bilateral'}
    addParameter(arg,'lowerT',0.5);
    addParameter(arg,'upperT',0.9);
    parse(arg,varargin{:});

    prob_image = imregional(prob_map, initLabel);
    mergeLabel = initLabel;
    
    % Static thresholds
    lowerT = arg.Results.lowerT;
    upperT = arg.Results.upperT;
    
    % Dynamic thresholds
    if strcmpi(arg.Results.thresh_policy, 'Dynamic')
        [min_prob, max_prob] = bounds(prob_image(:));
        lowerT = (1 - lowerT) * min_prob + lowerT * max_prob;
        upperT = (1 - upperT) * min_prob + upperT * max_prob;       
    end
    
    % Unilateral merging: background suppression
    mergeLabel(prob_image <= lowerT) = 0;
    
    % Bilateral merging: foreground and background downscaling
    if strcmpi(arg.Results.merge_policy, 'Bilateral')
        mergeLabel(prob_image > upperT) = max(initLabel(:))+1;
    end
    
    mergeLabel = process_labels(mergeLabel);
    mergeN = length(unique(mergeLabel));
   
end
