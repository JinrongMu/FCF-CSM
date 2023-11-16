function varargout = eval_binary_boundary(y_true, y_pred, varargin)
% eval_binary_boundary - Evaluate Binary Mask Boundaries.
%
% Syntax
% =================
% 1) Multiple evaluation metrics
%    metrics = eval_binary_boundary(y_true, y_pred);
%    metrics = eval_binary_boundary(y_true, y_pred, 'EdgeWidth',[10,5]);
%    metrics = eval_binary_boundary(y_true, y_pred, 'EdgeWidth',[10,5], ...
%        'TargetLabel',255,'EvalMetrics','all','debug_mode',2);
%    metrics = eval_binary_boundary(y_true, y_pred, 'EdgeWidth',[10,5], ...
%        'ProbThresh',0.5,'EvalMetrics','all','debug_mode',2);
% 2) Single evaluation metric
%    [value, name] = eval_binary_boundary(y_true, y_pred, 'EdgeWidth', [10,5],...
%        'EvalMetrics', 'BIOU', 'debug_mode', 2);
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
% References
% =================
% [1] B. Cheng, R. Girshick, P. Doll¨¢r, A. C. Berg, and A. Kirillov, 
% "Boundary IoU: Improving object-centric image segmentation evaluation." 
% CVPR 2021, pp. 15334-15342.
% 
% Abbreviations
% =================
% MIOU          Mask IoU
% TIOU          Trimap IoU
% FM            F-measure
% BIOU          Boundary IoU
% BCR           Boundary Correct Rate
% BER           Boundary Error Rate
% BLF           Boundary Local Fit
%
% Variables
% =================
% G             Binary Ground Truth Mask.
% P             Prediction Binary Mask.
% G1, P1        A set of pixels on the contour line of the binary mask.
% Gd, Pd        A set of pixels in the boundary region of the binary mask
% d             The pixel width of the boundary region.

% Parameter Initialization
% =========================================================
arg = inputParser; fun_name = 'eval_binary_boundary'; 
addParameter(arg,'EdgeWidth',[10,5]); 
addParameter(arg,'TargetLabel',0); 
addParameter(arg,'ProbThresh',0.5);
addParameter(arg,'EvalMetrics','all'); 
addParameter(arg,'debug_mode',0);
parse(arg,varargin{:});

% Parameters Check
defaultMetrics = {'MIOU','TIOU','FM','BIOU','BCR','BER','BLF'};
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

% 2. Get the edge contour: Gd (GT) and Pd (Prediction).
conn = 4; 
G1 = boundarymask(G, conn); P1 = boundarymask(P, conn);

EdgeWidth = arg.Results.EdgeWidth;
if EdgeWidth(1) > 0
    stru_elem = strel('disk', EdgeWidth(1));
    G_inner = imerode(G, stru_elem);
    P_inner = imerode(P, stru_elem);
else
    G_inner = G; P_inner = P;  
end
if EdgeWidth(2) > 0
    stru_elem = strel('disk', EdgeWidth(2));
    G_outer = imdilate(G, stru_elem);
    P_outer = imdilate(P, stru_elem);
else
    G_outer = G; P_outer = P;  
end  
Gd = G_outer & (~G_inner); Pd = P_outer & (~P_inner);
P1_correct = P1 & Gd; P1_error =  P1 - P1_correct;
Pd_correct = Pd & Gd; Pd_error =  Pd - Pd_correct; 

% 3. Calculating Evaluation Metrics.
mIOU = nnz(G&P) / (nnz(G|P) + eps);         % Mask IoU
tIOU = nnz(Gd&P) / (nnz(Gd&(G|P)) + eps);   % Trimap IoU
p = nnz(P1&Gd) / (nnz(P1) + eps); r = nnz(G1&Pd) / (nnz(G1) + eps);
FM = 2 * p * r / (p + r + eps);             % F-measure
bIOU = nnz(Gd&Pd) / (nnz(Gd|Pd) + eps);     % Boundary IoU
BCR = nnz(P1_correct) / (nnz(P1)+eps);      % Boundary Correct Rate
BER = nnz(P1_error) / (nnz(P1)+eps);        % Boundary Error Rate
BLF = nnz(Gd&Pd_correct) / (nnz(Gd|Pd_correct)+eps); % Boundary Local Fit 

% Output Settings
% =========================================================
valueMetrics = [mIOU, tIOU, FM, bIOU, BCR, BER, BLF];
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
    
    figure; nrows = 2; ncols = 3;
    
    subplot(nrows, ncols, 1), imshow(G, []), title('GT');
    subplot(nrows, ncols, ncols+1), imshow(P, []), title('Prediction');
    xstr = sprintf("Mask IOU = %.2f%%", 100*mIOU); xlabel(xstr);
    
    subplot(nrows, ncols, 2), imshow(Gd, []), title('Gd');
    subplot(nrows, ncols, ncols+2), imshow(Pd, []), title('Pd');
    xstr = sprintf("Boundary IOU = %.2f%%", 100*bIOU); xlabel(xstr);
    
    subplot(nrows, ncols, 3), imshow(Pd_correct, []), title('Correct Edge');
    xstr = sprintf("BLF = %.2f%%", 100*BLF); xlabel(xstr);
    subplot(nrows, ncols, ncols+3), imshow(Pd_error, []), title('Error Edge');
    xstr = sprintf("BER = %.2f%%", 100*BER); xlabel(xstr);
  
end

end  