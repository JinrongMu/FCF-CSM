% Calculate Hue, Saturation, and Luminance Based on a Chromaticity Diagram.
% H - Hue after rotation based on reference white.
% S - Ratio of chromaticity coordinates to reference white to sRGB edge.
% Y - Luminance as defined in the chromaticity diagram
function varargout = rgb2hsy(arr, varargin)
%
% Syntax
% =================
% [Hue, Sat, Y] = rgb2hsy(arr,'cieType','CIE1931','StartAngle',-60,...
%     'ClockWise',true,'hueType','Hue','edgeType','sRGB','debug_mode',0);
% HSY = rgb2hsy(arr,'cieType','CIE1931','StartAngle',-60,...
%     'ClockWise',true,'hueType','Hue','edgeType','sRGB','debug_mode',0);
%
% Input Arguments
% =================
% arr           RGB Array, uint8 [0, 255] or float [0, 1].
% cieType       Chromaticity Diagram Standards, ['CIE1931','CIE1960','CIE1976'].
% StartAngle    Zero degree hue angle, the default is -60.
% ClockWise     Hue Angle Rotation Direction, the default is true.
% hueType       The type of hue value, {'Hue', 'HueNorm', 'HueProb'}
% edgeType      Edge of saturation curve, sL(spectral line) or sRGB.
% debug_mode    Print debug information, 0: Silent, 1: Call information, 2: Call details.
%
% Output Arguments
% =================
% out           HSY Array or [Hue, Sat, Y]
%
% Hue           Custom Hue
% Sat           Custom Saturation
% Y             Custom Luminance

% Parameter Initialization
% =========================================================  
arg = inputParser; fun_name = 'rgb2hsy';
addParameter(arg,'cieType','CIE1931');  
addParameter(arg,'StartAngle',-60); 
addParameter(arg,'ClockWise',true);
addParameter(arg,'hueType','Hue'); 
addParameter(arg,'edgeType','sRGB'); 
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

cieType = arg.Results.cieType;
StartAngle = arg.Results.StartAngle;
ClockWise  = arg.Results.ClockWise;
hueType   = arg.Results.hueType;
edgeType   = arg.Results.edgeType;

% Method Implementation    
% =========================================================

% RGB -> Yxy
arr_Yxy = rgb2Yxy(arr, 'cieType', cieType, 'debug_mode', debug_mode-1);

FlipDims = (size(arr_Yxy,3) == 1);
if FlipDims, arr_Yxy = permute(arr_Yxy,[1,3,2]); end

% Calculate Luminance (Y) and Chroma (xy)
Y = arr_Yxy(:,:,1); arr_xy = arr_Yxy(:,:,2:3); 

% Calculate Hue
chromDiag = ChromaDiagram(cieType);
wp_xy = chromDiag.getRefWhitePoint();

Hue = calc_hue_angle(arr_xy, 'Centers', wp_xy, 'StartAngle', StartAngle, ...
    'ClockWise', ClockWise, 'debug_mode', debug_mode-1);

if strcmpi(hueType, 'HueNorm') || strcmpi(hueType, 'HueProb')
    Hue = Hue / 360;
    if strcmpi(hueType, 'HueProb')   
        mu = 0.8109; sigma = 0.0428;
        HDI = exp(-(Hue-mu).^2/(2*sigma^2))/(sqrt(2*pi)*sigma);
        Hue = HDI / 10;
    end
end

% Calculate Saturation
x1 = arr_xy(:,:,1) - wp_xy(1);
y1 = arr_xy(:,:,2) - wp_xy(2);
d1 = sqrt(x1 .^ 2 + y1 .^ 2); 

if strcmpi(edgeType, 'sRGB')
    [x2, y2] = calc_srgb_edges(arr_xy, 'cieType', cieType, ...
        'inputType', 'ChromaCoord', 'debug_mode', debug_mode-1);
%     stdHue = calc_hue_angle(arr_xy, 'Centers', ref_wp, 'StartAngle', 0, ...
%         'debug_mode', debug_mode-1);
%     [x2, y2] = calc_srgb_edges(stdHue, 'cieType', cieType, ...
%         'inputType', 'HueAngle', 'debug_mode', debug_mode-1);
end
d2 = sqrt((x2 - wp_xy(1)).^2 + (y2 - wp_xy(2)).^2); 
Sat = d1 ./ d2; Sat(isnan(Sat)) = 0;
% Process pure black (0,0,0) -> 0(S)
Sat(Y==0) = 0;

% Output Settings
% =========================================================
if nargout == 3
    varargout = {Hue, Sat, Y};
else
    out = cat(3, Hue, Sat, Y);
    if FlipDims, out = permute(out,[1,3,2]); end  
    varargout = {out};
end
   
end

%% Helper functions for calculating hue and saturation.
function [x, y] = calc_srgb_edges(arr, varargin)
%calc_srgb_edges - Calculate the intersection with the sRGB color gamut.
%
% Syntax
% =================
% [x, y] = calc_srgb_edges(arr,'cieType','CIE1931',...
%     'inputType','ChromaCoord','debug_mode',0);
%
% Input Arguments
% =================
% cieType       Chromaticity Diagram Standards, ['CIE1931','CIE1960','CIE1976'].
% inputType     Input type, {'HueAngle', 'ChromaCoord'}.
% debug_mode    Print debug information, 0: Silent, 1: Call information, 2: Call details.
%
% Output Arguments
% =================
% x             Horizontal coordinates. 
% y             Vertical coordinates.

% Parameter Initialization
% =========================================================
arg = inputParser; fun_name = 'calc_srgb_edges';
addParameter(arg,'cieType','CIE1931'); 
addParameter(arg,'inputType','ChromaCoord');  
addParameter(arg,'debug_mode', 0);
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
chromDiag = ChromaDiagram(arg.Results.cieType);
srgb_xy = chromDiag.getSRGBGamutVertices();
wp_xy = chromDiag.getRefWhitePoint();
srgb_hue = calc_hue_angle(srgb_xy, 'Centers', wp_xy);

coords.red = srgb_xy(1,:);
coords.green = srgb_xy(2,:);
coords.blue = srgb_xy(3,:);
coords.white_point = wp_xy;

angles.red = srgb_hue(1);
angles.green = srgb_hue(2);
angles.blue = srgb_hue(3);

% Radial saturation
if strcmpi(arg.Results.inputType, 'ChromaCoord')
    [k1, b1] = calc_point_slope(arr, coords.white_point);
    A = calc_hue_angle(arr, 'Centers', coords.white_point);
elseif strcmpi(arg.Results.inputType, 'HueAngle')
    if max(arr)>360
        warning("Angles over 360 degrees are truncated.")
        arr = mod(arr, 360);
    end
    k1 = tan(arr ./ 180 .* pi);
    b1 = -1 * coords.white_point(1) .* k1 + coords.white_point(2);
    A = arr;
else
    error("Function '%s' has undefined parameter 'inputType': '%s'.", ...
        fun_name, arg.Results.inputType);
end

% sRGB gamut edge
k2 = zeros(size(k1)); b2 = zeros(size(b1)); 
mask = A >= angles.red & A < angles.green;
[k_rg,b_rg] = calc_point_slope(coords.red,coords.green);
k2(mask) = k_rg; b2(mask) = b_rg;

mask = A >= angles.green & A < angles.blue;
[k_gb,b_gb] = calc_point_slope(coords.green,coords.blue);
k2(mask) = k_gb; b2(mask) = b_gb;

mask = A >= angles.blue | A < angles.red;
[k_br,b_br] = calc_point_slope(coords.red,coords.blue);
k2(mask) = k_br; b2(mask) = b_br;

if numel(A)==1 && A==90
    x = coords.white_point(1); y = k_rg * x + b_rg;
elseif numel(A)==1 && A==270
    x = coords.white_point(1); y = k_br * x + b_br;
else
    [x, y] = calc_cross_point(k1, b1, k2, b2);
    if nnz(A==90)>0
        x(A==90) = coords.white_point(1);
        y(A==90) = k_rg * coords.white_point(1) + b_rg;
    end
    
    if nnz(A==270)>0
        x(A==270) = coords.white_point(1);
        y(A==270) = k_br * coords.white_point(1) + b_br;
    end
end

end

function [x, y] = calc_cross_point(k1, b1, k2, b2)

    x = -1 * (b2 - b1) ./ (k2 - k1);
    y = (b1 .* k2 - k1 .* b2) ./ (k2 - k1);

end

function [k, b] = calc_point_slope(p1, p2)

    if size(p1,3)==1, p1 = permute(p1,[1,3,2]); end
    if size(p2,3)==1, p2 = permute(p2,[1,3,2]); end

    k = (p2(:,:,2) - p1(:,:,2)) ./ (p2(:,:,1) - p1(:,:,1));
    % k(isinf(k))=0;
    b = -1 * k .* p1(:,:,1) + p1(:,:,2);

end

function angles = calc_hue_angle(arr, varargin)
%calc_hue_angle -  Compute Hue Angle Helper Function.
%
% Syntax
% =================
% angles = calc_hue_angle(arr_xy, 'Centers', [0,0], ...
%     'StartAngle', 0, 'ClockWise', false, 'debug_mode', 0);
%
% Input Arguments
% =================
% arr           Chromaticity Coordinates xy Array, float.
% Centers       Center point of chromaticity coordinates.
% StartAngle    Zero degree hue angle.
% ClockWise     Clockwise or not, the default is false.
% debug_mode    Print debug information, 0: Silent, 1: Call information, 2: Call details.
%
% Output Arguments
% =================
% angles        Custom Hue Angle. 

% Parameter Initialization
% ========================================================= 
arg = inputParser; fun_name = 'calc_hue_angle';
addParameter(arg,'Centers',[0,0]); 
addParameter(arg,'StartAngle',0); 
addParameter(arg,'ClockWise',false);
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
if FlipDims, arr = permute(arr, [1,3,2]); end

x = arr(:,:,1) - arg.Results.Centers(1);
y = arr(:,:,2) - arg.Results.Centers(2);
angles = atan2(y, x) * 180 / pi;
angles = mod(angles - arg.Results.StartAngle, 360);
if arg.Results.ClockWise
    angles = 360 - angles; 
end

% Process pure black (0,0,0) or pure white (1,1,1)
blackMask = ((arr(:,:,1)==0) & (arr(:,:,2)==0));
whiteMask = ((x==0) & (y==0));
zeroMask = blackMask | whiteMask;
angles(zeroMask) = 0;

end
