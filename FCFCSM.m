% function LabelImage = FCFCSM(img_uint8)
clear
close all
clc

addpath('./Data/') 
addpath('./Functions/') 
addpath('./Modules/Superpixels/') 
addpath('./Modules/Clustering/') 
addpath('./Modules/Colors/') 
addpath('./Modules/PWS/') 


debug_mode = 2;
% -------------------------------------------------------------------------
% Step 1. Load Image and Preprocessing
% -------------------------------------------------------------------------
image = imread('pws_00001.JPG'); % 5792-1-1
label = imread('pws_00001_Lesion.png');

crop_guide = 'None';  % {'None', 'GT', 'CDI', 'Hue', 'HDI'}
[img_proc, lab_proc, bbox, size_scale] =  PreprocessImages(image, label, ...
    'crop_guide', crop_guide, 'padSize', 80, 'padType', 'Loose', ...
    'size_policy', 'PixelNum', 'size_limit', 400, 'size_tol', 20, ...
    'size_align', 'None', 'debug_mode', debug_mode);


% -------------------------------------------------------------------------
% Step 2. Superpixel Partition and Merging
% -------------------------------------------------------------------------
super_method = 'gGMMSP'; super_size = 50; 
merge_method = 'ColorV3'; % {'ColorV3', 'CDIV1', 'None'}

img_gauss = imsmooth(img_proc); img_float = double(img_gauss) / 255.;
N0 = round(numel(img_gauss) / 3 / super_size);

[L1, N1] = ext_superpixels(img_gauss, N0, 'method', super_method, ...
    'debug_mode', 0);

[L2, N2] = pixel_blocks_merge(img_gauss, L1, 'merge_method', ...
    merge_method, 'debug_mode', debug_mode);

% -------------------------------------------------------------------------
% Step 3. Extract Region Features
% -------------------------------------------------------------------------
featureID = 'F5V5';%'F7V1', 'F7V3', 'F5V4';'F5V5'; 'F8V2'
[F, V, G] = extract_features(img_gauss, featureID, 'debug_mode', debug_mode);
[~, X] = imregional(F, L2);

% -------------------------------------------------------------------------
% Step 4. CGFFCM-Based Clustering
% -------------------------------------------------------------------------
n_cluster = 5;
if n_cluster == 0
    seed_policy = 'DPC-Auto';
else
    seed_policy = 'DPC';
end

if (n_cluster == 0) && strcmpi(seed_policy, 'DPC-Auto')
    [sp_dp, c_idx, K] = cluster_DPC(X,n_cluster,'version','v1','debug_mode',0);  
elseif (n_cluster > 0) && strcmpi(seed_policy, 'DPC')
    [sp_dp, c_idx, K] = cluster_DPC(X,n_cluster,'version','v2','debug_mode',0);
else % 'Random'
    sp_dp = []; c_idx = randperm(length(X), n_cluster); K = n_cluster; 
end
sp_fcm = RGFFCM(X, K, V, G, 'c_idx', c_idx, 'debug_mode', debug_mode);

lab_fcm = sp_fcm(L2);



% -------------------------------------------------------------------------
% Step 5. Gamut Decision and Post-processing
% -------------------------------------------------------------------------
decision_method = 'LABV1';  % {'CDIV2', 'LABV1'}
lab_cluster = process_labels(lab_fcm);
[gtLabels, tarLabels] = cluster_decision(img_proc, lab_proc, ...
    lab_cluster, 'method', decision_method, 'debug_mode', debug_mode);
decision_acc = nnz(ismember(tarLabels,gtLabels))/length(gtLabels);

lab_decision = ismember(lab_cluster, gtLabels);
% lab_decision = ismember(lab_cluster, tarLabels);

if debug_mode == 2
    lab_manual = uint8(lab_decision)*255;
    lab_auto = uint8(ismember(lab_cluster, tarLabels))*255;
    
    metrics1 = eval_binary_boundary(lab_proc,lab_manual,'EdgeWidth',[10,5], ...
        'TargetLabel',255,'EvalMetrics','All','debug_mode',0);
    metrics2 = eval_binary_boundary(lab_proc, lab_auto, 'EdgeWidth',[10,5], ...
        'TargetLabel',255,'EvalMetrics','All','debug_mode',0);
    
    figure;
    subplot(2,2,1), imshow(img_proc), title('Image');
    subplot(2,2,2), imshow(lab_proc), title('GT');
    
    subplot(2,2,3), imshow(lab_manual)
    xstr = sprintf("mIOU=%.2f%%,bIOU=%.2f%%,BER=%.2f%%,BLF=%.2f%%", ...
        100*metrics1.MIOU, 100*metrics1.BIOU, 100*metrics1.BER, 100*metrics1.BLF);
    xlabel(xstr), title('Manual');
    
    subplot(2,2,4), imshow(lab_auto)
    xstr = sprintf("mIOU=%.2f%%,bIOU=%.2f%%,BER=%.2f%%,BLF=%.2f%%", ...
        100*metrics2.MIOU, 100*metrics2.BIOU, 100*metrics2.BER, 100*metrics2.BLF);
    xlabel(xstr), title('Automatic');
end

% Remove small objects from the foreground and background.
P = round(numel(lab_decision) * 0.0005); conn = 4;
binary_mask = bwareaopen(lab_decision, P, conn);
binary_mask = 1 - bwareaopen(~binary_mask,P,conn);
lab_post = uint8(binary_mask) * 255;

% Label Reconstruction
label_final = zeros(size(label), 'uint8');
x = bbox(1); y = bbox(2); w = bbox(3); h = bbox(4);
if size_scale > 0
    label_roi = imresize(lab_post, [h, w], 'nearest');
else
    label_roi = lab_post;
end
label_final(y : y+h-1, x : x+w-1) = label_roi;

if debug_mode == 2
    metrics = eval_binary_boundary(label, label_final, 'EdgeWidth', [10,5], ...
        'TargetLabel',255,'EvalMetrics','all','debug_mode',debug_mode-1);
end

% end