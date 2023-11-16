function varargout = rgb2Yxy(arr, varargin)
% rgb2Yxy - Calculation of luminance and chromaticity coordinates under different CIE standards.
%
% Syntax
% =================
% Yxy = rgb2Yxy(arr_rgb);
% Yxy = rgb2Yxy(arr_rgb, 'cieType', 'CIE1976');
% Yxy = rgb2Yxy(arr_rgb, 'cieType', 'CIE1976', 'debug_mode', 1);
% [arr_Yxy, arr_xyz] = rgb2Yxy(arr_rgb, 'cieType', 'CIE1976', 'debug_mode', 1);
%
% Input Arguments
% =================
% arr           RGB Array, uint8 [0, 255] or float [0, 1].
% cieType       Chromaticity Diagram Standards, ['CIE1931','CIE1960','CIE1976'].
% debug_mode    Print debug information, 0: Silent, 1: Call information, 2: Call details.
%
% Output Arguments
% =================
% out           Yxy Array or XYZ Array

% Parameter Initialization
% =========================================================
arg = inputParser; fun_name = 'rgb2Yxy';
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

% uint8 -> float
if isa(arr, 'uint8')
    arr = double(arr) / 255.;
end

% RGB -> XYZ
arr_xyz = rgb2xyz(arr);

% XYZ -> Yxy
arr_Yxy = xyz2Yxy(arr_xyz, 'cieType', arg.Results.cieType, ...
    'debug_mode', arg.Results.debug_mode-1);

% Output Settings
% =========================================================
if nargout == 2
    varargout = {arr_Yxy, arr_xyz};
else
    varargout = {arr_Yxy};
end

end