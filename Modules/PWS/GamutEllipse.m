% Chromaticity Statistical Model (CSM)
function varargout = GamutEllipse(identifier, varargin)
% GamutEllipse - Covariance Ellipse Fitting Parameters.
% 
% Syntax
% =================
% ell = GamutEllipse('CIE1931-1Sigma');
% [x y] = GamutEllipse('CIE1931-1Sigma', 'num_points', 1000);
% [center, axes, angle] = GamutEllipse('CIE1931-1Sigma');
%
% Input Arguments
% =================
% identifier	Identifier of the elliptic model.
% num_points    The number of ellipse coordinate points.
% debug_mode    Print debug information, 0: Silent, 1: Call information, 2: Call details.
%
% Output Arguments
% =================
% out           ell struct or [x y] or [center, axes, angle]
%
% ell           Elliptic parameters in the form of structures.
% center        Coordinates of the ellipse center.
% axes          The length of the axis of the ellipse.
% angle         The rotation angle of the ellipse.

% Parameter Initialization
% =========================================================
arg = inputParser; fun_name = 'GamutEllipse';
addParameter(arg,'num_points',1000);
addParameter(arg,'debug_mode',0);
parse(arg,varargin{:});

if arg.Results.debug_mode == 1
    fprintf("\nCall functions:\t%s('%s')\n", fun_name, identifier)
elseif arg.Results.debug_mode == 2
    fprintf("\nCall functions:\t%s('%s')\n", fun_name, identifier)
    fprintf('\nDefault Parameters:\n'); disp(arg.Results);
end

% Method Implementation    
% =========================================================
ell_id = lower(identifier);
switch ell_id
    case 'cie1931'
        center=[0.38161240,0.33639282];axes=[0.09500340,0.18772064];angle=106.60669708;
    case 'cie1960'
        center=[0.25175714, 0.32131788];axes=[0.04768985, 0.15592331];angle=99.86700439;
    case 'cie1976'
        center=[0.24919438, 0.48175445];axes=[0.07039746, 0.14328287];angle=106.38826752; 
    case 'cie1931-thresh1'      % 100th percentile
        center=[0.422332, 0.337247];axes=[0.124325, 0.314064];angle=99.593544;
    case 'cie1931-thresh10'
        center=[0.418769, 0.337485];axes=[0.116605, 0.274305];angle=99.514809;
    case 'cie1931-thresh100'
        center=[0.414046, 0.339062];axes=[0.104364, 0.237594];angle=100.398834;
    case 'cie1931-thresh1000'
        center=[0.410152, 0.340793];axes=[0.090336, 0.198324];angle=100.886482;
    case 'cie1931-thresh17'     % 75th percentile
        center=[0.417087, 0.337846];axes=[0.115116, 0.268349];angle=99.780502;
    case 'cie1931-thresh494'    % 50th percentile
        center=[0.409966, 0.340483];axes=[0.095619, 0.208950];angle=99.930878;
    case 'cie1931-1sigma'
        center=[0.417823, 0.334388];axes=[0.065945, 0.158484];angle=99.281525;
    case 'cie1931-2sigma'
        center=[0.417823, 0.334388];axes=[0.118701, 0.285272];angle=99.281525; 
    case 'cie1931-3sigma'
        center=[0.417823, 0.334388];axes=[0.197836, 0.475453];angle=99.281525;
    case 'cie1960-1sigma'
        center=[0.276593, 0.322232];axes=[0.030448, 0.122358];angle=99.343834;
    case 'cie1960-2sigma'
        center=[0.276593, 0.3222322];axes=[0.054806, 0.220244];angle=99.343834; 
    case 'cie1960-3sigma'
        center=[0.276593, 0.3222322];axes=[0.091343, 0.367074];angle=99.343834;
    case 'cie1976-1sigma'
        center=[0.276338, 0.483609];axes=[0.043686, 0.122027];angle=104.863647;
    case 'cie1976-2sigma'
        center=[0.276338, 0.483609];axes=[0.078635, 0.219648];angle=104.863647;  
    case 'cie1976-3sigma'
        center=[0.276338, 0.483609];axes=[0.131059, 0.366080];angle=104.863647;  
    otherwise
        error("Undefined identifier '%s' in GamutEllipse.", identifier);
end

% Output Settings
% =========================================================
if nargout == 3
    varargout = {center, axes, angle};
elseif nargout == 2
    [x, y] = calc_ellipse_coords(center, axes, angle, arg.Results.num_points);
    varargout = {x, y};
else
    ell = struct('center',center,'axes',axes,'angle',angle);
    varargout = {ell};
end

end

function varargout = calc_ellipse_coords(center, axes, angle, N)
% calc_ellipse_coords - Calculate Ellipse Coordinates.
% 
% Syntax
% =================
% [x, y] = calc_ellipse_coords(center, axes, angle);
% points = calc_ellipse_coords(center, axes, angle, 1000);
%
% Input Arguments
% =================
% center	The coordinates of the center of the ellipse.
% axes      Long axis (X-axis) and short axis (Y-axis).
% angle     Angle of rotation.
% N         The number of coordinate points, which defaults to 1000.
%
% Output Arguments
% =================
% points	struct
% x         [N¡Á1 double]
% y         [N¡Á1 double]

% Parameter Initialization
% =========================================================
if ~exist('N', 'var')
    N = 1000;
end

% Method Implementation    
% =========================================================
a = axes(1) / 2;  
b = axes(2) / 2;  % semi-axis
rot = angle * pi / 180.;
theta = linspace(0, 2*pi, N+1); theta = theta(1:N);

xy = transpose([a.*cos(theta); b.*sin(theta)]);
rot_matrix = transpose([cos(rot), -sin(rot); sin(rot), cos(rot)]);

xy = xy * rot_matrix + center;
points.x = xy(:,1); points.y = xy(:,2);

% Output Settings
% =========================================================
if nargout == 2
    varargout = {points.x, points.y};
else
    varargout = {points};
end

end
