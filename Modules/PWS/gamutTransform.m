function varargout = gamutTransform(arr, varargin)
% gamutTransform - Elliptical Gamut Transformation
% 
% Syntax
% =================
% arr_gamut = gamutTransform(arr,'ell_id','CIE1931-1Sigma','alpha',0.5,'debug_mode',1);
%
% Input Arguments
% =================
% arr           RGB Array, uint8 [0, 255] or float [0, 1].
% ell_id        Gamut Ellipse Identifier, such as 'CIE1931-1Sigma'.
% alpha         The base of the exponential mapping, the default is 0.5.
% debug_mode    Print debug information, 0: Silent, 1: Call information, 2: Call details.
%
% Output Arguments
% =================
% out           p or {p, D} or {p, D, arr_xyz, ell}.
%
% p             Gamut mapped array. 
% D             Mahalanobis distance map.
% arr_xyz       XYZ Array.
% ell           Elliptic parameters in the form of structures.

% Parameter Initialization
% ========================================================= 
arg = inputParser; fun_name = 'gamutTransform';
addParameter(arg,'ell_id','CIE1931-1Sigma');
addParameter(arg,'alpha',0.5);   
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

% RGB -> Yxy
cieType = arg.Results.ell_id(1:7);	% ['CIE1931' 'CIE1960' 'CIE1976']
[arr_Yxy, arr_xyz] = rgb2Yxy(arr, 'cieType', cieType, 'debug_mode', debug_mode-1);

% Process parameter
FlipDims = (size(arr_Yxy,3) == 1);
if FlipDims, arr_Yxy = permute(arr_Yxy,[1,3,2]); end

ell = GamutEllipse(arg.Results.ell_id,'debug_mode', debug_mode-1);
x = arr_Yxy(:,:,2) - ell.center(1);
y = arr_Yxy(:,:,3) - ell.center(2);
a = ell.axes(1) / 2; b = ell.axes(2) / 2;
theta = ell.angle * pi / 180.;  %  theta = (180-ell.angle) * pi / 180.; % 

% Calculate the Mahalanobis distance
A = (b * cos(theta))^2 + (a * sin(theta))^2;
B = (b^2 - a^2) * sin(2 * theta);
C = (b * sin(theta))^2 + (a * cos(theta))^2;
R = (a * b)^2;
D = (A * (x.^2) + B * (x .* y) + C * (y.^2)) / R;
D = sqrt(D);

% Exponential Distance Mapping
alpha = arg.Results.alpha;
p = alpha .^ D;

% Output Settings
% =========================================================
if nargout == 4
    varargout = {p, D, arr_xyz, ell};
elseif nargout == 3
    varargout = {p, D, arr_xyz};
elseif nargout == 2
    varargout = {p, D};
else
    varargout = {p};
end

end

