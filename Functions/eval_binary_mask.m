function varargout = eval_binary_mask(y_true, y_pred, varargin)
%eval_binary_mask -  Evaluate Binary Mask.
%
% Syntax
% =================
% 1) Multiple evaluation metrics
%     metrics = eval_binary_mask(y_true, y_pred);
%     metrics = eval_binary_mask(y_true, y_pred, 'TargetLabel', 255, 'debug_mode',2);    
%     metrics = eval_binary_mask(y_true, y_pred, 'ProbThresh', 0.5, 'debug_mode',2);
%       
% 2) Single evaluation metric
%     [value, name] = eval_binary_mask(y_true, y_pred, ...
%         'EvalMetrics', 'JA', 'debug_mode', 2);
%
% Input Arguments
% =================
% y_true        Binary image, the value is {0, 1}, {0, 255} or logical value.
% y_pred        Evaluate images, probability maps or binary maps.
% EdgeWidth     Pixel width for image erosion and dilation,[innerb, outerb].
% TargetLabel	The label of the target class, such as 1 or 255.
% ProbThresh    Threshold for probability map binarization, such as 0.5.
% EvalMetrics   Abbreviation for Evaluation Metrics.
% debug_mode    Control debug information, 0: Silent, 1: Call information, 2: Call details.	
%
% Output Arguments
% =================
% out           metrics or {values, names}
%
% metrics       Evaluation metrics struct.
% values        Evaluation metrics values.
% names         Evaluation metrics names.
%
% Abbreviations
% =================
% JA            Jaccard
% DI          	Dice
% AC            Accuracy
% SE          	Sensitivity, Recall
% PR          	Precision
% SP            Specificity
% ROI           roi ratio

% Parameter Initialization
% =========================================================
arg = inputParser; fun_name = 'eval_binary_mask'; 
addParameter(arg,'TargetLabel',0); 
addParameter(arg,'ProbThresh',0.5);
addParameter(arg,'EvalMetrics','all'); 
addParameter(arg,'debug_mode',0);
parse(arg,varargin{:});

% Parameter Check
defaultMetrics = {'JA','DI','AC','SE','PR','SP','ROI'};
EvalMetrics = arg.Results.EvalMetrics;
if strcmpi(EvalMetrics, 'all')
    metricIdxs = length(EvalMetrics) + 1;
elseif any(strcmpi(EvalMetrics, defaultMetrics))
    metricIdxs = find(strcmpi(EvalMetrics, defaultMetrics));
else
    fprintf("Available 'EvalMetrics' are: \n"); disp(defaultMetrics);
    error("Undefined parameter 'EvalMetrics': '%s'.", EvalMetrics);
end


% Method Implementation 
% =========================================================

% 1. Get the binary mask: G (GT) and P (Prediction).
TargetLabel = arg.Results.TargetLabel;
ProbThresh = arg.Results.ProbThresh;
if TargetLabel ~= 0
    G = (y_true == TargetLabel); P = (y_pred == TargetLabel);
elseif ProbThresh > 0
    G = (y_true > ProbThresh); P = (y_pred > ProbThresh);
else
    G = (y_true ~= 0); P = (y_pred ~= 0);
end

% 2. Calculate mask-based evaluation metrics: [JA DI AC SE PR SP ROI]
TN = nnz(~(G | P)); FN = nnz(G & (~P)); 
FP = nnz((~G) & P); TP = nnz(G & P); 
numCorrect = nnz(G==P); numError = nnz(G~=P); %#ok
numSamples = numel(G); 
numTruthFg = nnz(G); numTruthBg = nnz(~G); 
numPredFg = nnz(P); numPredBg = nnz(~P);       %#ok

JA = TP / (TP + FN + FP + eps);
DI = 2 * TP / (2 * TP + FN + FP + eps); 
AC = numCorrect / (numSamples);
SE = TP / (numTruthFg + eps);
PR = TP / (numPredFg + eps);
SP = TN / (numTruthBg + eps);
ROI = numTruthFg / (numSamples); 

% Output Settings
% =========================================================
valueMetrics = [JA DI AC SE PR SP ROI];
if metricIdxs > length(EvalMetrics)   
    values = valueMetrics;
    names = defaultMetrics;
else
    values = valueMetrics(metricIdxs);
    names = defaultMetrics{metricIdxs};
end

if nargout == 2
    varargout = {values, names};
else
    propertyNames = names;
    propertyValues = num2cell(values); 
    metrics = cell2struct(propertyValues, propertyNames, 2);
    varargout = {metrics};
end

% Debug Information
% =========================================================
if arg.Results.debug_mode == 1
    fprintf('\nCall functions:\t%s\n', fun_name)
elseif arg.Results.debug_mode == 2 
    fprintf('\nCall functions:\t%s\n', fun_name)
    fprintf('----------------------------------------');
    fprintf('\nDefault Parameters:\n'); disp(arg.Results);
    
    figure; nrows = 1; ncols = 2;
    
    subplot(nrows, ncols, 1), imshow(G, []), title('GT');
    subplot(nrows, ncols, 2), imshow(P, []), title('Prediction');
    xstr = sprintf("JA=%.2f%%,SE=%.2f%%,PR=%.2f%%,ROI=%.2f%%", ...
        100*JA, 100*SE, 100*PR, 100*ROI); xlabel(xstr); 
    
end

end
    