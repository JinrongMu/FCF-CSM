function H = imhomogeneity(I, varargin)
% imhomogeneity - Computing the local homogeneity of an image.
%
% Syntax
% =================
% H = imhomogeneity(I);
% H = imhomogeneity(I, 'grad_space', 'gray', 'grad_op', 'default', 'homo_space', 'hsv');
% H = imhomogeneity(I, 'grad_space', 'gray', 'grad_op', 'default', ...
%     'homo_space', 'hsv', 'ksize', 5, 'norm_level', 'image-wise', 'debug_mode', 1);
% H = imhomogeneity(I, 'grad_space', 'homo', 'grad_op', 'sobel', ...
%     'homo_space', 'hsv', 'ksize', 5, 'norm_level', 'image-wise', 'debug_mode',1);
%
% Input Arguments
% =================
% I             RGB Image, float [0,1] or uint8 [0, 255].
% grad_space    Gradient calculation space, {'gray', 'homo'}.
% grad_op       Gradient operator,{'sobel', 'prewitt','roberts','laplacian', 'default'}.
% homo_space    Homogeneous computing space, {'hsv', 'lab', 'rgb'}
% ksize         Local window size, odd number, default is 5.
% norm_level    Normalization strategy, {'channel-wise', 'image-wise'}
% debug_mode    Print debug information, 0: Silent, 1: Call information, 2: Call details.
%
% Output Arguments
% =================
% H             Local homogeneity
%
% Reference papers 
% =================
% [1] H. D. Cheng, and Y. Sun, 
%     "A hierarchical approach to color image segmentation using homogeneity," 
%     IEEE Trans. Image Process., vol. 9, no. 12, pp. 2071-2082, 2000.
% [2] T. Sa?, and M. ?unka?, 
%     "Color image segmentation based on multiobjective artificial bee colony optimization,"  
%     Appl. Soft Comput., vol. 34, pp. 389-401, 2015.
% [3] Z. Zhou, X. Zhao, and S. Zhu,
%     "K-harmonic means clustering algorithm using feature weighting for 
%     color image segmentation," Multimedia Tools and Applications, 
%     vol. 77, no. 12, pp. 15139-15160, 2018.
% [4] A. G. Oskouei, M. Hashemzadeh, B. Asheghi, and M. A. Balafar, 
%     "CGFFCM: CLuster-weight and group-local feature-weight learning in 
%     fuzzy C-means clustering algorithm for color image segmentation," 
%     Appl. Soft Comput., vol. 113, pp. 108005, 2021.
%
% Algotithm
% =================
% H = 1 - E.* V
% E             Discontinuity, Normalize magnitude of the gradient.   
% V             Normalize standard deviation.

% Parameter Initialization    
% =========================================================
arg = inputParser; fun_name = 'imhomogeneity'; 
addParameter(arg,'grad_space','gray');
addParameter(arg,'grad_op','sobel');
addParameter(arg,'homo_space','hsv');
addParameter(arg,'ksize',5);
addParameter(arg,'norm_level','channel-wise');
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

% RGB -> HSV
if isa(I, 'uint8'), I = double(I) / 255.; end

if strcmpi(arg.Results.homo_space, 'rgb')
    img_homo = I;
elseif strcmpi(arg.Results.homo_space, 'lab')
    img_homo = rgb2lab(I);
else
    img_homo = rgb2hsv(I);
end

% Calculate discontinuity
grad_op = arg.Results.grad_op;
if strcmpi(arg.Results.grad_space, 'gray')
    img_gray = rgb2gray(I);
    E = imgradient(img_gray, 'grad_op', grad_op, 'debug_mode', debug_mode-1);
else
    E = imgradient(img_homo, 'grad_op', grad_op, 'debug_mode', debug_mode-1);
end

% Calculate the normalized standard deviation
pool_size = [arg.Results.ksize arg.Results.ksize];
V = impool(img_homo, 'std', 'pool_size', pool_size, 'debug_mode', debug_mode-1);

% Normalized gradient and standard deviation.
if strcmpi(arg.Results.norm_level, 'image-wise')
    E = ext_normalize(E, 'max');
    V = ext_normalize(V, 'max');
else
    if size(E,3) == 1, E = repmat(E, 1, 1 ,3); end
    for i = 1:3
        E(:,:,i) = ext_normalize(E(:,:,i), 'max');
        V(:,:,i) = ext_normalize(V(:,:,i), 'max');
    end
end
    
% Calculate homogeneity
H = 1 - E .* V;

end