function varargout = PreprocessImages(image, label, varargin)
% Crop ROI + Padding + Limet Size
%
% Syntax
% =================
% [img_proc, lab_proc, bounding_box] =  PreprocessImages(image, label, ...
%     'crop_guide', 'GT', 'padSize', 100, 'padType', 'Loose', ...
%     'size_policy', 'PixelNum', 'size_limit', 500, 'size_tol', 20, ...
%     'size_align', 'None', 'debug_mode', debug_mode);
%

    % Parameter Initialization
    % =========================================================
    arg = inputParser; fun_name = 'PreprocessImages';
    addParameter(arg,'crop_guide','GT');        % {'None', 'GT', 'CDI', 'Hue', 'HDI'}
    addParameter(arg,'padSize',100);
    addParameter(arg,'padType','Loose');        % {'Loose', 'Tight'}
    addParameter(arg,'size_policy','PixelNum');	% {'PixelNum', 'SizeWidth', 'None'}
    addParameter(arg,'size_limit',500);
    addParameter(arg,'size_tol',20);
    addParameter(arg,'size_align','None');      % {'Resize', 'None'}
    addParameter(arg,'debug_mode',0);
    parse(arg,varargin{:});

    debug_mode = arg.Results.debug_mode;
    if debug_mode == 1
        fprintf('\nCall functions:\t%s\n', fun_name)
    elseif debug_mode == 2
        fprintf('\nCall functions:\t%s\n', fun_name)
        fprintf('----------------------------------------');
        fprintf('\nDefault Parameters:\n'); disp(arg.Results);
    end

    % ---------------------------------------------------------------------
    % Crop ROI
    % ---------------------------------------------------------------------
    crop_guide = arg.Results.crop_guide;
    padSize = arg.Results.padSize; 
    padType = arg.Results.padType;
    [img_roi, lab_roi, bounding_box, prob_map] = crop_pws_image(image, ...
        label, 'crop_guide', crop_guide, 'padSize', padSize, ...
        'padType', padType, 'debug_mode', debug_mode-1);
    
    % ---------------------------------------------------------------------
    % Restrict Image Size
    % ---------------------------------------------------------------------
    size_policy = lower(arg.Results.size_policy);
    size_limit = arg.Results.size_limit;
    size_tol = arg.Results.size_tol;
    size_align = arg.Results.size_align;
    [img_resize,lab_resize,size_scale] = restrict_image_size(img_roi,lab_roi,...
        'size_policy', size_policy, 'size_limit', size_limit,...
        'size_tol', size_tol, 'size_align', size_align, 'debug_mode', debug_mode-1);

    % Output Settings
    % =========================================================
    if nargout == 4
        varargout = {img_resize, lab_resize, bounding_box, size_scale};
    elseif nargout == 3
        varargout = {img_resize, lab_resize, bounding_box};
    elseif nargout == 2
        varargout = {img_resize, lab_resize};
    else
        varargout = {img_resize};
    end
    
    % Debug Information
    % =========================================================
    if debug_mode == 2     
        figure;
        subplot(2,2,1), imshow(image);
        rectangle('Position', bounding_box, 'EdgeColor', 'b', 'LineWidth', 1);
        xlabel(mat2str(bounding_box)); title('Image + ROI');
        subplot(2,2,2), imshow(prob_map, []), title(crop_guide);
        xlabel(mat2str(size(prob_map)));
        subplot(2,2,3), imshow(img_resize), title('Resized Image');
        xlabel(mat2str(size(img_resize)));
        subplot(2,2,4), imshow(lab_resize, []), title('Resized Label');
        xlabel(mat2str(size(lab_resize)));
    end

end

function varargout = crop_pws_image(image, label, varargin)

% Syntax
% =================
% [img_roi, lab_roi, bounding_box] = crop_pws_image(image, label, ...
%     'crop_guide', crop_guide, 'padSize', padSize, ...
%     'padType', padType, 'debug_mode', debug_mode);


    % Parameter Initialization
    % =========================================================    
    arg = inputParser; fun_name = 'crop_pws_image';
    addParameter(arg,'crop_guide','GT');    % {'None', 'GT', 'CDI', 'Hue', 'HDI'}
    addParameter(arg,'padSize',100);
    addParameter(arg,'padType','Loose');    % {'Loose', 'Tight'}
    addParameter(arg,'debug_mode',0);
    parse(arg,varargin{:});
    
    debug_mode = arg.Results.debug_mode;
    if debug_mode == 1
        fprintf('\nCall functions:\t%s\n', fun_name)
    elseif debug_mode == 2
        fprintf('\nCall functions:\t%s\n', fun_name)
        fprintf('----------------------------------------');
        fprintf('\nDefault Parameters:\n'); disp(arg.Results);
    end
    
    % Method Implementation
    % =========================================================
    
    % 1. get crop mask
    crop_guide = arg.Results.crop_guide;
    if strcmpi(crop_guide, 'GT')
        prob_map = label; mask = prob_map; 
    elseif strcmpi(crop_guide, 'CDI')
        prob_map = get_features(image, 'CDI', 'CIE1931-1Sigma');
        mask = prob2mask(prob_map, 0.6);
    elseif strcmpi(crop_guide, 'HDI')
        prob_map = get_features(image, 'Hue', 'HueProb');
        mask = prob2mask(prob_map, 0.75);
    elseif strcmpi(crop_guide, 'Hue')
        prob_map = get_features(image, 'Hue', 'HueNorm');
        mask = prob2mask(prob_map, 0.8);
    else  % crop_guide: 'None'
        mask = []; prob_map = label; 
    end
    
    % 2. crop and padding
    if isempty(mask)
        [x,y,w,h] = deal(1, 1, size(image,2), size(image,1));
        img_roi = image; lab_roi = label;
    else
        padSize = arg.Results.padSize; padType = arg.Results.padType;
        [x,y,w,h] = calc_bounding_box(mask, 'padSize', padSize, ...
            'padType', padType, 'debug_mode', debug_mode-1);
        img_roi = image(y:y+h-1,x:x+w-1,:);
        lab_roi = label(y:y+h-1,x:x+w-1);
    end
    
    % Output Settings
    % =========================================================
    bounding_box = [x,y,w,h];
    if nargout == 4 
        varargout = {img_roi, lab_roi, bounding_box, prob_map};
    elseif nargout == 3
        varargout = {img_roi, lab_roi, bounding_box};
    elseif nargout == 2
        varargout = {img_roi, lab_roi};
    else
        varargout = {img_roi};
    end
    
    % Debug Information
    % =========================================================
    if (debug_mode == 2) && (~isempty(mask))
        figure;
        subplot(2,2,1), imshow(image);
        rectangle('Position', [x,y,w,h], 'EdgeColor', 'b', 'LineWidth', 1);
        xlabel(mat2str([x,y,w,h])); title('Image + ROI');
        subplot(2,2,2), imshow(prob_map, []), title(crop_guide);
        xlabel(mat2str(size(prob_map)));
        subplot(2,2,3), imshow(img_roi), title('ROI Image');
        xlabel(mat2str(size(img_roi)));
        subplot(2,2,4), imshow(lab_roi, []), title('ROI Label');
        xlabel(mat2str(size(lab_roi)));
    end

end

function mask = prob2mask(prob_map, prob_thresh)

    binary_mask = prob_map > prob_thresh;

    % Remove small objects from the foreground.
    P = round(numel(prob_map) * 0.0005); conn = 4;
    binary_mask = bwareaopen(binary_mask, P, conn);     
    mask = uint8(binary_mask) * 255;
end

function varargout = restrict_image_size(image, label, varargin)

    % Parameter Initialization
    % =========================================================  
    arg = inputParser; fun_name = 'restrict_image_size';
    addParameter(arg,'size_policy','PixelNum'); % {'PixelNum', 'SizeWidth', 'None'}
    addParameter(arg,'size_limit',500);
    addParameter(arg,'size_tol',20);
    addParameter(arg,'size_align','Resize');    % {'Resize', 'None'}
    addParameter(arg,'debug_mode',0);
    parse(arg,varargin{:});

    debug_mode = arg.Results.debug_mode;
    if debug_mode == 1
        fprintf('\nCall functions:\t%s\n', fun_name)
    elseif debug_mode == 2
        fprintf('\nCall functions:\t%s\n', fun_name)
        fprintf('----------------------------------------');
        fprintf('\nDefault Parameters:\n'); disp(arg.Results);
    end

    % Method Implementation
    % =========================================================
    [rows, cols, ~] = size(image);
    size_policy = lower(arg.Results.size_policy);
    size_limit = arg.Results.size_limit;
    size_tol = arg.Results.size_tol;
    size_align = arg.Results.size_align;
    
    size_scale = 0;
    switch size_policy
        case 'sizewidth'
            max_size = max(rows, cols);
            if (max_size - size_limit)^2 > (size_tol^2)
                size_scale = round(size_limit / max_size, 6);
            end
        case 'pixelnum'
            num_pixel = rows * cols;
            max_pixel = (size_limit + size_tol)^2;
            min_pixel = (size_limit - size_tol)^2;
            if (num_pixel < min_pixel) || (num_pixel > max_pixel)
                pixel_limit = size_limit^2;
                size_scale = round(sqrt(pixel_limit / num_pixel), 6);
            end
        otherwise
            fprintf("Undefined Parameter 'size_policy': %s\n",size_policy);
    end
    
    if size_scale > 1 && strcmpi(size_align, 'None')
        size_scale = 0;
    end
    if size_scale == 0 || size_scale == 1
        img_resize = image; lab_resize = label; 
    else
        img_resize = imresize(image, size_scale, 'bilinear');
        lab_resize = imresize(label, size_scale, 'nearest');
    end
    
    % Output Settings
    % =========================================================
    if nargout == 3
        varargout = {img_resize, lab_resize, size_scale};
    else
        varargout = {img_resize, lab_resize};       
    end

end