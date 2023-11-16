function img_gauss = imsmooth(img_uint8, varargin)
% imsmooth - Smooth an image using a Gaussian filter.
%
% Syntax
% =================
% img_gauss = imsmooth(img_uint8);
% img_gauss = imsmooth(img_uint8,'sigma',2,'ksize',5,'domain','spatial','debug_mode',1);
%
% Input Arguments
% =================
% img_uint8 	Smooth image, uint8. 
% sigma         Standard deviation of the Gaussian distribution, default 2.
% ksize         Size of the Gaussian filter, default size is 2*ceil(2*sigma)+1.
% domain        omain in which to perform filtering, ['auto','spatial','frequency'].
% debug_mode    Print debug information, 0: Silent, 1: Call information, 2: Call details.
%
% Output Arguments
% =================
% img_gauss 	smoothed image.

% Parameter Initialization    
% =========================================================
arg = inputParser; fun_name = 'imsmooth';
addParameter(arg,'sigma',2);
addParameter(arg,'ksize',5); 
addParameter(arg,'domain','spatial');
addParameter(arg,'debug_mode',0);
parse(arg,varargin{:});

if arg.Results.debug_mode == 1
    fprintf('\nCall functions:\t%s\n', fun_name)
elseif arg.Results.debug_mode == 2
    fprintf('\nCall functions:\t%s\n', fun_name)
    fprintf('----------------------------------------');
    fprintf('\nDefault Parameters:\n'); disp(arg.Results);
end

sigma = arg.Results.sigma;
ksize = arg.Results.ksize;
domain = arg.Results.domain;

% Method Implementation    
% =========================================================
if size(img_uint8, 3) == 1
    img_gauss = imgaussfilt(img_uint8,sigma,...
        'FilterSize',ksize,'FilterDomain',domain);
else
    img_gauss = zeros(size(img_uint8), 'like', img_uint8);
    for i = 1:size(img_uint8, 3)
        Iblur = imgaussfilt(img_uint8(:,:,i),sigma,...
            'FilterSize',ksize,'FilterDomain',domain);
        img_gauss(:,:,i) = Iblur;
    end
end

end