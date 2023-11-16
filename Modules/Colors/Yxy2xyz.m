function varargout = Yxy2xyz(arr, varargin)
% Yxy2xyz - Convert chromaticity coordinates Yxy to CIE XYZ (double).
%
% Syntax
% ================
% arr_xyz = Yxy2xyz(arr)
% arr_xyz = Yxy2xyz(arr, 'cieType', 'CIE1931')
% arr_xyz = Yxy2xyz(arr, 'cieType', 'CIE1931', 'debug_mode', 1)
% [X, Y, Z] = Yxy2xyz(arr, 'cieType', 'CIE1931', 'debug_mode', 1)
%
% Input Arguments
% =================
% arr           Yxy Array, Y in [0,100], x and y in [0,1]. 
% cieType       Chromaticity Diagram Standards, ['CIE1931','CIE1960','CIE1976'].
% debug_mode    Print debug information, 0: Silent, 1: Call information, 2: Call details.
%
% Output Arguments
% =================
% out        	XYZ Array or [X, Y, Z].

% Parameter Initialization
% =========================================================
arg = inputParser; fun_name = 'Yxy2xyz';
addParameter(arg,'cieType','CIE1931');
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

FlipDims = (size(arr,3) == 1);
if FlipDims, arr = permute(arr,[1,3,2]); end

% Yuv -> Yxy
cieType = lower(arg.Results.cieType);
if strcmp(cieType, 'cie1960')
    c_u = 3; c_v = 2; coeffs = [2 -8 4];
elseif strcmp(cieType, 'cie1976')
    c_u = 9; c_v = 4; coeffs = [6 16 12];
else
    assert(strcmp(cieType, 'cie1931'));
end

if strcmp(cieType, 'cie1960') || strcmp(cieType, 'cie1976')
    arr_sum = coeffs(1)*arr(:,:,2) + coeffs(2)*arr(:,:,3) + coeffs(3);
    arr(:,:,2) = c_u * arr(:,:,2) ./ arr_sum;
    arr(:,:,3) = c_v * arr(:,:,3) ./ arr_sum;
    arr(isnan(arr)) = 0;
end

% Yxy -> XYZ
out = zeros(size(arr));
out(:,:,2) = arr(:,:,1) / 100;
out(:,:,1) = out(:,:,2) .* arr(:,:,2) ./ arr(:,:,3);
out(:,:,3) = out(:,:,2) .* (1-arr(:,:,2)-arr(:,:,3)) ./ arr(:,:,3);
out(isnan(out)) = 0;

% Output Settings
% =========================================================
if nargout == 3
    varargout = {out(:,:,1), out(:,:,2), out(:,:,3)};
else
    if FlipDims, out = permute(out,[1,3,2]); end
    varargout = {out};
end

end