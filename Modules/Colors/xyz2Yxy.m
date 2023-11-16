function varargout = xyz2Yxy(arr, varargin)
% xyz2Yxy - Calculation of luminance and chromaticity coordinates under different CIE standards.
%
% Syntax
% =================
% Yxy = xyz2Yxy(arr_xyz);
% Yxy = xyz2Yxy(arr_xyz, 'cieType', 'CIE1976');
% Yxy = xyz2Yxy(arr_xyz, 'cieType', 'CIE1976', 'debug_mode', 1);
% [Y, x, y] = xyz2Yxy(arr_xyz, 'cieType', 'CIE1976', 'debug_mode', 1);
%
% Input Arguments
% =================
% arr           XYZ Array, float, [0, 1.08].
% cieType       Chromaticity Diagram Standards, ['CIE1931','CIE1960','CIE1976'].
% debug_mode    Control debug information, 0: Silent, 1: Call information, 2: Call details.
%
% Output Arguments
% =================
% out        	Yxy Array or [Y, x, y]

% Parameter Initialization
% =========================================================
arg = inputParser; fun_name = 'xyz2Yxy'; 
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

FlipDims = (size(arr, 3) == 1);
if FlipDims, arr = permute(arr,[1,3,2]); end

% Calculate Yxy
out = zeros(size(arr));
cieType = lower(arg.Results.cieType);
switch cieType
    case 'cie1931'
        c_u = 1; c_v = 1; coeffs = [1.,1.,1.];
    case 'cie1960'
        c_u = 4; c_v = 6; coeffs = [1.,15.,3.];
    case 'cie1976'
        c_u = 4; c_v = 9; coeffs = [1.,15.,3.];
    otherwise
        error("Undefined CIE standard 'cieType': %s", arg.Results.cieType);
end

arr_sum = coeffs(1)*arr(:,:,1) + coeffs(2)*arr(:,:,2) + coeffs(3)*arr(:,:,3);
out(:,:,1) = 100 * arr(:,:,2);
out(:,:,2) = c_u * arr(:,:,1) ./ arr_sum;
out(:,:,3) = c_v * arr(:,:,2) ./ arr_sum;
% Special Cases: black xyz (0,0,0) -> Yxy (0,0,0)
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