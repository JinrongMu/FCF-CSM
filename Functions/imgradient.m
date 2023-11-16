function varargout = imgradient(I, varargin)
% imgradient - Calculate the gradient, magnitude, and direction of an image.
%
% Syntax
% =================
% M = imgradient(I);
% [M, A] = imgradient(I);
% [M, A, Fx, Fy] = imgradient(I);
% [M, A, Fx, Fy] = imgradient(I,'grad_op','sobel','debug_mode',1);
%
% Input Arguments
% =================
% I             Input Array, gray or color image.
% grad_op       Gradient operator,{'sobel', 'prewitt', 'roberts','laplacian', 'default'}
% debug_mode    Print debug information, 0: Silent, 1: Call information, 2: Call details.
%
% Output Arguments
% =================
% out        	M or {M, A} or {M, A, Fx, Fy}.
% 
% M             Gradient magnitude
% A             Gradient direction, [-pi, pi].
% Fx, Fy        Numerical gradients

% Parameter Initialization    
% =========================================================
arg = inputParser; fun_name = 'imgradient';  
addParameter(arg,'grad_op','sobel');
addParameter(arg,'debug_mode',0);
parse(arg,varargin{:});

if arg.Results.debug_mode == 1
    fprintf("\nCall Function:\t%s\n", fun_name)
elseif arg.Results.debug_mode == 2
    fprintf('\nCall Function:\t%s\n', fun_name)
    fprintf('----------------------------------------');
    fprintf('\nDefault Parameters:\n'); disp(arg.Results);   
end

% Method Implementation    
% =========================================================

% uint8 -> float
% if isa(I, 'uint8'), I = double(I) / 255.; end
if isa(I, 'uint8'), I = double(I); end

% Calculate the gradient
grad_op = lower(arg.Results.grad_op);
opts = {'sobel', 'prewitt', 'roberts','laplacian'};
if any(strcmp(grad_op, opts))
    if strcmp(grad_op, 'roberts')
        h1 = [-1 0; 0 1]; h2 = [0 -1; 1 0];
    elseif strcmp(grad_op, 'laplacian')
        h1 = fspecial('laplacian',0);
        h2 = fspecial('laplacian',1);
    else
        h1 = fspecial(grad_op); h2 = h1';
    end
    Fx = imfilter(I, h1, 'replicate');
    Fy = imfilter(I, h2, 'replicate');
else
    [Fx, Fy] = gradient(I);
end

M = sqrt(Fx .* Fx + Fy .* Fy);
A = atan(Fy ./ Fx);

% Output Settings
% =========================================================
if nargout == 4
    varargout = {M, A, Fx, Fy};
elseif nargout == 2
    varargout = {M, A};
else
    varargout = {M};
end

end