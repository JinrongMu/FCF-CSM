function out = ext_normalize(arr, method, varargin)
% ext_normalize - Extended matlab function 'normalize'.
% 
% Syntax
% =================
% N = ext_normalize(arr, 'max');
% N = ext_normalize(arr, 'scale');
%
% Input Arguments
% =================
% arr           The array to process. 
% method        normalized method, {'max','sum','zscore','norm','scale','range','center'}.
% debug_mode    Print debug information, 0: Silent, 1: Call information, 2: Call details.
%
% Output Arguments
% =================
% out        	Normalized array.

% Parameter Initialization    
% =========================================================
arg = inputParser; fun_name = 'ext_normalize';
addParameter(arg,'debug_mode',0); 
parse(arg,varargin{:});

debug_mode = arg.Results.debug_mode;
if (debug_mode == 1) || (debug_mode == 2)
    fprintf('\nCall functions:\t%s\n', fun_name)
end

if ~exist('method', 'var')
    method = 'max';
else
    method = lower(method);
    opts = {'max','sum','zscore','norm','scale','range','center'};
    if ~any(strcmp(opts, method))
        fprintf("The available options for parameter '%s' are:\n", 'method');
        disp(opts)
        error("Normalization method '%s' is undefined in function '%s'.", ...
            method, fun_name)
    end
end

% Method Implementation    
% =========================================================
if strcmp(method, 'max')
    out = arr/max(arr(:));
elseif strcmp(method, 'sum')
    out = arr/sum(arr(:));
else
    out = normalize(arr, method);
end

end