function varargout = Yxy2rgb(arr, varargin)
% Yxy2rgb - Convert Yxy to RGB (double, uint8)
%
% Syntax
% =================
% arr_rgb = Yxy2rgb(arr)
% arr_rgb = Yxy2rgb(arr, 'cieType', 'CIE1931')
% arr_rgb = Yxy2rgb(arr, 'cieType', 'CIE1931', 'OutputType', 'uint8', 'debug_mode', 1)
% [arr_rgb, arr_xyz] = Yxy2rgb(arr, 'cieType', 'CIE1931', 'OutputType', 'uint8', 'debug_mode', 1)
%
% Input Arguments
% =================
% arr           Yxy Array, Y in [0,100], x and y in [0,1]. 
% cieType       Chromaticity Diagram Standards, ['CIE1931','CIE1960','CIE1976'].
% OutputType    Data type of returned RGB values, ['double', 'uint8', 'uint16'].
% debug_mode    Print debug information, 0: Silent, 1: Call information, 2: Call details.
%
% Output Arguments
% =================
% out        	RGB Array or XYZ Array.

% Parameter Initialization    
% =========================================================
arg = inputParser; fun_name = 'Yxy2rgb';
addParameter(arg,'cieType','CIE1931');
addParameter(arg,'OutputType', 'double');
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

% Yxy -> XYZ
arr_xyz = Yxy2xyz(arr, 'cieType', arg.Results.cieType, ...
    'debug_mode', arg.Results.debug_mode-1);

% XYZ -> RGB
arr_rgb = xyz2rgb(arr_xyz, 'OutputType', arg.Results.OutputType);   
if strcmpi(arg.Results.OutputType, 'double')
    arr_rgb(arr_rgb<0) = 0;
    arr_rgb(arr_rgb>1) = 1;
end

% Output Settings
% =========================================================
if nargout == 2
    varargout = {arr_rgb, arr_xyz};
else
    varargout = {arr_rgb};
end

end