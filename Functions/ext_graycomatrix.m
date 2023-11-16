function varargout = ext_graycomatrix(I, varargin)
% ext_graycomatrix - Extend the graycomatrix of matlab.
%
% Syntax
% =================
% img_glcm = ext_graycomatrix(I,'OffsetAngle',45);
% [img_glcm,glcms] = ext_graycomatrix(I,'OffsetAngle',45);
% [img_glcm,glcms,SI] = ext_graycomatrix(I,'OffsetAngle',45,'NumLevels',8,'isNorm',true);
%
% Input Arguments
% =================
% I             Single channel 2D image, (uint8 or double).
% OffsetAngle   OffsetAngle,{0, 45, 90, 135, 180}
% NumLevels     Number of gray levels, specified as an integer, [1, 8].
% isNorm        Whether to normalize the output image.
% debug_mode    Print debug information, 0: Silent, 1: Call information, 2: Call details.
%
% Output Arguments
% =================
% img_glcm      Probabilistic output image for graycomatrix.
% glcms         Gray-level co-occurrence matrix (or matrices)
% SI            Scaled image used in calculation of GLCM
%
% References
% =================
% [1] Haralick, R.M., K. Shanmugan, and I. Dinstein, 
% "Textural Features for Image Classification", IEEE Transactions on 
%  Systems, Man, and Cybernetics, Vol. SMC-3, 1973, pp. 610-621.

% Parameter Initialization
% =========================================================
arg = inputParser; fun_name = 'ext_graycomatrix'; 
addParameter(arg,'OffsetAngle',180); 
addParameter(arg,'NumLevels',8);
addParameter(arg,'isNorm',true);
addParameter(arg,'debug_mode',0);
parse(arg,varargin{:});

assert(size(I,3)==1, "The input to '%s' must be a 2D image.", fun_name);

if arg.Results.debug_mode == 1
    fprintf('\nCall Function:\t%s\n', fun_name)
elseif arg.Results.debug_mode == 2
    fprintf('\nCall Function:\t%s\n', fun_name)
    fprintf('----------------------------------------');
    fprintf('\nDefault Parameters:\n'); disp(arg.Results);       
end

% Method Implementation    
% =========================================================
offset = [0 1;-1 1;-1 0;-1 -1]; 
NumLevels = arg.Results.NumLevels;
[glcms, SI] = graycomatrix(I, 'Offset', offset, 'Symmetric', true, ...
    'NumLevels', NumLevels, 'GrayLimits', []);

[H, W] = size(I); img_glcm = zeros([H,W], 'double'); 
SI_cols = transpose(im2col(SI, [2 2]));

switch arg.Results.OffsetAngle
    case 0
        i=1; glcm_i = glcms(:,:,i); id_sub = SI_cols(:,[1 3]);
        id = sub2ind([NumLevels NumLevels], id_sub(:,1), id_sub(:,2));
        img_glcm(1:end-1,1:end-1) = reshape(glcm_i(id), [H-1 W-1]);
    case 45
        i=2; glcm_i = glcms(:,:,i); id_sub = SI_cols(:,[2 3]);
        id = sub2ind([NumLevels NumLevels], id_sub(:,1), id_sub(:,2));
        img_glcm(2:end,1:end-1) = reshape(glcm_i(id), [H-1 W-1]);
    case 90
        i=3; glcm_i = glcms(:,:,i); id_sub = SI_cols(:,[2 1]);
        id = sub2ind([NumLevels NumLevels], id_sub(:,1), id_sub(:,2));
        img_glcm(2:end,1:end-1) = reshape(glcm_i(id), [H-1 W-1]);
    case 135
        i=4; glcm_i = glcms(:,:,i); id_sub = SI_cols(:,[4 1]);
        id = sub2ind([NumLevels NumLevels], id_sub(:,1), id_sub(:,2));
        img_glcm(2:end,2:end) = reshape(glcm_i(id), [H-1 W-1]);
    otherwise
        img_glcm = SI;
end

if arg.Results.isNorm
    img_glcm = img_glcm ./ max(img_glcm(:));
end

% Output Settings
% =========================================================
if nargout == 3
    varargout = {img_glcm, glcms, SI};
elseif nargout == 2
    varargout = {img_glcm, glcms};
else
    varargout = {img_glcm};
end

end