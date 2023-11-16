function varargout = process_labels(inLabels, varargin)
% process_labels - Process the label so that it goes from 1 to N.
% 
% Syntax
% =================
% outLabels = process_labels(inLabels);
% [outLabels, numLabels] = process_labels(inLabels,'method','label_mapping','debug_mode',1); 
%
% Input Arguments
% =================
% inLabels      Array of labels to be processed. 
% method        The method of handling labels, {'LabelMapping','DataSinking'}.
% debug_mode    Print debug information, 0: Silent, 1: Call information, 2: Call details.
%
% Output Arguments
% =================
% outLabels     Array of processed labels
% numLabels     Number of label classes after processing.

% Parameter Initialization    
% =========================================================
arg = inputParser; fun_name = 'process_labels';
addParameter(arg,'method','LabelMapping');
addParameter(arg,'debug_mode',0);              
parse(arg,varargin{:});


% Method Implementation
% =========================================================
inLabels = int32(inLabels);
if strcmpi(arg.Results.method, 'LabelMapping')
    outLabels = process_labels_v2(inLabels);
else
    outLabels = process_labels_v1(inLabels);
end

if arg.Results.debug_mode == 1
    fprintf('\nCall functions:\t%s\n', fun_name)
elseif arg.Results.debug_mode == 2
    fprintf('\nCall functions:\t%s\n', fun_name)
    fprintf('----------------------------------------');
    fprintf('\nDefault Parameters:\n'); disp(arg.Results);
    warning('Labels are discontinuous and empty clusters appear.')
    [minLabel, maxLabel] = bounds(inLabels(:));
    fprintf('\nProcess Labels: [%d,%d]->[1,%d]\n',minLabel,maxLabel,numLabels);
end


% Output Control
if nargout == 2   
    numLabels = length(unique(outLabels));
    varargout = {outLabels, numLabels};
else
    varargout = {outLabels};
end

end


%% V2: Label Mapping
function outLabels = process_labels_v2(inLabels)

    % in_label = int32(in_label);

    % Make sure that the minimum label is 1.
    Labels = unique(inLabels); [minLabel, maxLabel] = bounds(Labels);
    if minLabel ~= 1
        % fprintf("Warning: The labels contains: %d.\n", L_min);
        inLabels = inLabels - minLabel + 1;
        Labels = unique(inLabels); [minLabel, maxLabel] = bounds(Labels);
    end
    
    % Modify discontinuous labels
    if length(Labels) == (maxLabel - minLabel + 1)
        outLabels = inLabels;
    else
        labelMap = [minLabel:maxLabel]'; [~, labelIdx] = sort(Labels);
        for i = 1:length(Labels)
            labelMap(Labels(i)) = labelIdx(i);
        end
        outLabels = labelMap(inLabels);
    end

end

function out_label = process_labels_v3(in_label) %#ok
    unique_labels = unique(in_label);
    num_labels = length(unique_labels);
    [~, sorted_indices] = sort(unique_labels);
    label_map = zeros(num_labels, 1);
    label_map(sorted_indices) = 1:num_labels;
    out_label = label_map(in_label);
end

function out_label = process_labels_v4(in_label) %#ok
    unique_labels = unique(in_label);
    num_labels = length(unique_labels);
    out_label = zeros(size(in_label));
    for i = 1:num_labels
        out_label(ismember(in_label, unique_labels(i))) = i;
    end
end

%% V1: Data Sinking
function outLabels = process_labels_v1(inLabels) 

    % in_label = int32(in_label);

    % Make sure that the minimum label is 1.
    minLabel = min(inLabels(:)); 
    if minLabel == 1
        outLabels = inLabels;    
    else
        % fprintf("Warning: The labels contains: %d.\n", label_min);
        outLabels = inLabels - (minLabel - 1);
    end

    % Make sure the labels are consecutive integers.
    Labels = unique(outLabels);
    if (length(Labels)~= max(Labels)) && (length(Labels)>1)
        currLabel = 1;
        for i=2:length(Labels)
            deltaLabels = Labels(i) - Labels(i-1);
            if deltaLabels > 1
                outLabels = outLabels - (deltaLabels-1)*int32(outLabels>currLabel);
            end
            currLabel = currLabel + 1;
        end
    end

end
