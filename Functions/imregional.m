function varargout = imregional(img_float, img_label, varargin)
% imregional - Image partition statistics.
%
% Syntax
% =================
% region_map = imregional(img_float, img_label);
% [region_map, region_feature] = imregional(img_float, img_label);
% [region_map, region_feature, pixel_count] = imregional(img_float, img_label);
% [region_map, region_feature, pixel_count, img_label] = imregional(img_float, img_label);
%
% Input Arguments
% =================
% img_float         Single or multi-channel processing of images..
% img_label         Partition labels for images.
% debug_mode        Print debug information, 0: Silent, 1: Call information, 2: Call details.
%
% Output Arguments
% =================
% region_map        Area-averaged images, such as colorized labels.
% region_feature	Area-averaged arrays, such as cmap array.
% pixel_count       Region pixel counts.
% img_label         Labeled images after processing.
% debug_mode        Print debug information, 0: Silent, 1: Call information, 2: Call details.
% 
% Parameter Initialization    
% =========================================================
arg = inputParser; fun_name = 'imregional';
addParameter(arg,'debug_mode',0);              
parse(arg,varargin{:});

debug_mode = arg.Results.debug_mode;
if (debug_mode == 1) || (debug_mode == 2)
    fprintf('\nCall functions:\t%s\n', fun_name)
end

% Method Implementation    
% =========================================================

% Process labels
[img_label, numLabels] = process_labels(img_label, 'debug_mode', debug_mode);

  
% Calculate the mean array after regionalization (cmap array).
idxs = label2idx(img_label);
region_feature = zeros(numLabels,size(img_float,3)); 
pixel_count = zeros(numLabels,1);
for i = 1:numLabels
    pixel_count(i) = length(idxs{i});  
    for j = 1:size(img_float,3)
        img_slice = img_float(:,:,j);
        region_feature(i,j) = mean(img_slice(idxs{i}));
    end
end

% Calculate the mean image after regionalization (colorLabel).
region_map = zeros(size(img_float), 'double');
for j = 1:size(img_float,3)
    slice_value = region_feature(:,j);
    region_map(:,:,j) = slice_value(img_label);
end

% Output Settings
% =========================================================
if nargout == 2
    varargout = {region_map, region_feature};
elseif nargout == 3
    varargout = {region_map, region_feature, pixel_count};
elseif nargout == 4
    varargout = {region_map, region_feature, pixel_count, img_label};
else
    varargout = {region_map};
end
   
end
