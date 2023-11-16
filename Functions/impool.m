function out = impool(I, pool_type, varargin)
% impool - Perform pooling on multi-channel images.
%
% Syntax
% =================
% out = impool(I, pool_type);
% out = impool(I, pool_type, 'pool_size', [5,5], 'debug_mode',1);
%
% Input Arguments
% =================
% I             Multi-channel images. 
% pool_type     Local operation type. {'max', 'min', 'mean', 'std', 'var', 'median', 'contrast'}
% 
% pool_size     Window size, the default is [5,5].
% debug_mode    Print debug information, 0: Silent, 1: Call information, 2: Call details.
%
% Output Arguments
% =================
% out        	Filtered image.

% Parameter Initialization    
% =========================================================
arg = inputParser; fun_name = 'impool';
addParameter(arg,'pool_size', [5,5]);
addParameter(arg,'debug_mode',0); 
parse(arg,varargin{:});

if arg.Results.debug_mode == 1
    fprintf("\nCall functions:\t%s(I, '%s')\n", fun_name, pool_type);
elseif arg.Results.debug_mode == 2
    fprintf("\nCall functions:\t%s(I, '%s')\n", fun_name, pool_type);
    fprintf('----------------------------------------');
    fprintf('\nDefault Parameters:\n'); disp(arg.Results);
end

% Method Implementation
% =========================================================

% Channel-by-channel pooling
out = zeros(size(I)); pool_size = arg.Results.pool_size;
for i = 1:size(I,3)
    out(:,:,i) = pool2d(I(:,:,i), pool_type, 'pool_size', pool_size, ...
        'debug_mode', arg.Results.debug_mode-1);
end

end

function out = pool2d(I, pool_type, varargin)
% pool2d - Perform pooling on single-channel images.
%
% Syntax
% =================
% out = pool2d(I, pool_type);
% out = pool2d(I, pool_type, 'pool_size', [5,5], 'debug_mode',1);
%
% Input Arguments
% =================
% I             Single-channel images. 
% pool_type     Local operation type, {'max', 'min', 'mean', 'std', 'var',
%                                      'median', 'contrast'}.
% 
% pool_size     Window size, the default is [5,5].
% debug_mode    Print debug information, 0: Silent, 1: Call information, 2: Call details.
%
% Output Arguments
% =================
% out        	Filtered image.

% Parameter Initialization    
% =========================================================
arg = inputParser; fun_name = 'pool2d';
addParameter(arg,'pool_size', [5,5]);
addParameter(arg,'debug_mode',0); 
parse(arg,varargin{:});

if arg.Results.debug_mode == 1
    fprintf("\nCall functions:\t%s(I, '%s')\n", fun_name, pool_type);
elseif arg.Results.debug_mode == 2
    fprintf("\nCall functions:\t%s(I, '%s')\n", fun_name, pool_type);
    fprintf('----------------------------------------');
    fprintf('\nDefault Parameters:\n'); disp(arg.Results);
end

% Method Implementation    
% =========================================================
pool_type = lower(pool_type);
pool_size = arg.Results.pool_size;
switch pool_type 
    case 'max'       
        out = colfilt(I,pool_size,'sliding',@max);
    case 'min'       
        out = colfilt(I,pool_size,'sliding',@min);
    case 'mean'       
        out = colfilt(I,pool_size,'sliding',@mean);
    case 'std'       
        out = colfilt(I,pool_size,'sliding',@std);
    case 'var'       
        out = colfilt(I,pool_size,'sliding',@var);
    case 'median'       
        out = colfilt(I,pool_size,'sliding',@median);
    case 'contrast'
        fun = @(x) max(x(:))-min(x(:));
        out = nlfilter(I,pool_size,fun);
    otherwise
        error("Undefined 'pool_type': '%s'.", pool_type);
end

end