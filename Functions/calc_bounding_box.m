function varargout = calc_bounding_box(mask, varargin)
% calc_bounding_box -  Calculate the bounding box of the mask.
%
% Syntax
% =================
% out = calc_bounding_box(mask);
% [x,y,w,h] = calc_bounding_box(mask);
% out = calc_bounding_box(mask,'padSize',50,'padType','loose','debug_mode',2);
%
% Input Arguments
% =================
% mask          Mask image, values in {0,1} or {0,255}.
% padSize       0:No padding,(0,1):Padding ratio; int: padsize.
% padType       Priority cropping method,{'Loose', 'Tight'}.
% debug_mode    Print debug information, 0: Silent, 1: Call information, 2: Call details.
%
% Output Arguments
% =================
% out           bounding box: [x, y, w, h].

% Parameter Initialization
% =========================================================
arg = inputParser; fun_name = 'calc_bounding_box';
addParameter(arg,'padSize',0);
addParameter(arg,'padType','loose');
addParameter(arg,'debug_mode',0);
parse(arg,varargin{:});

padSize = arg.Results.padSize;
padType = arg.Results.padType;

% Method Implementation    
% =========================================================

% Get the bounding box
[rowIndexs, colIndexs] = find(mask);
[xmin, xmax] = bounds(colIndexs);
[ymin, ymax] = bounds(rowIndexs);
boxWidth = xmax - xmin + 1;
boxHeigh = ymax - ymin + 1; 
oldBoundingBox = [xmin, ymin, boxWidth, boxHeigh];

% Add border to bounding box

% process padsize
if padSize <= 0
    padSize = 0;
elseif padSize < 1
    padSize = round(max(boxWidth, boxHeigh) * padSize);
end

if padSize > 0
    % Adjust the aspect ratio
    aspectRatio = boxWidth / boxHeigh;
    if strcmpi(padType, 'loose') && (aspectRatio > 4/3)
        xoffset = 2 * round(padSize * 0.25);
        yoffset = 2 * round(padSize * 0.75);
    elseif strcmpi(padType, 'loose') && (aspectRatio < 3/4)
        xoffset = 2 * round(padSize * 0.75);
        yoffset = 2 * round(padSize * 0.25);
    else
        [xoffset, yoffset] = deal(padSize);
    end
    
    % Update bounding box parameters
    [maskHeigh, maskWidth] = size(mask);
    xmin = max(xmin - xoffset, 1);
    ymin = max(ymin - yoffset, 1);
    xmax = min(xmax + yoffset, maskWidth);
    ymax = min(ymax + yoffset, maskHeigh);
    boxWidth = xmax - xmin + 1;
    boxHeigh = ymax - ymin + 1;   
end

% Output Settings
% =========================================================
newBoundingBox = [xmin, ymin, boxWidth, boxHeigh];
if nargout == 4
    varargout = {xmin, ymin, boxWidth, boxHeigh};
else
    varargout = {newBoundingBox};
end

% Debug Information
% =========================================================
if arg.Results.debug_mode > 0
    fprintf('\nCall functions:\t%s\n', fun_name)
    fprintf('----------------------------------------\n');
    if arg.Results.debug_mode == 2
        fprintf('Default Parameters:\n'); disp(arg.Results);
        fprintf('Mask Bounding Box:%s.\n', mat2str(oldBoundingBox));
    end
    fprintf('Padding Bounding Box:%s.\n', mat2str(newBoundingBox));
end

end

