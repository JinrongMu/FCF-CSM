function [F, V, G] = extract_features(img_uint8, featureID, varargin)
% extract_features - Get the features and coefficients of RGB images.
%
% Syntax: [F, V, G] = extract_features(img_uint8,'CGFFCM','debug_mode',1);

% Parameter Initialization    
% =========================================================
arg = inputParser; fun_name = 'extract_features';
addParameter(arg,'debug_mode',0); 
parse(arg,varargin{:});

if arg.Results.debug_mode == 1
    fprintf("\nCall functions:\t%s('%s')\n", fun_name, featureID);
end

% Method Implementation 
% =========================================================
if strcmpi(featureID, 'CGFFCM')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_homo = get_features(img_uint8, 'Homogeneity', 'v1');
    img_gabor = get_features(img_uint8, 'Texture', 'Gabor');
    img_glcm = get_features(img_uint8, 'Texture', 'GLCM');
    F = cat(3, img_lab, img_homo, img_gabor, img_glcm);
    V=[0.6 0.2 0.2]; G=[1 1 1 2 2 2 3 3];
elseif strcmpi(featureID, 'F2V1')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_pig = get_features(img_uint8, 'Skin', 'RGB');
    F = cat(3, img_lab(:,:,1), img_pig(:,:,1));
    V=[0.5 0.5]; G=[1 2];
elseif strcmpi(featureID, 'F3V1')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_pig = get_features(img_uint8, 'Skin', 'RGB');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    F = cat(3, img_lab(:,:,2), img_pig(:,:,1), img_cdi);
    V=[0.3 0.3 0.4]; G=[1 2 3];
elseif strcmpi(featureID, 'F3V2')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_pig = get_features(img_uint8, 'Skin', 'RGB');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    F = cat(3, img_lab(:,:,2), img_pig(:,:,1), img_cdi);
    V=[0.3 0.2 0.5]; G=[1 2 3];
elseif strcmpi(featureID, 'F3V3')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    F = img_lab;
    V=[1/3 1/3 1/3]; G=[1 1 1];
elseif strcmpi(featureID, 'F4V1')
    img_pig = get_features(img_uint8, 'Skin', 'RGB');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    img_hdi = get_features(img_uint8, 'Hue', 'HueProb');
    F = cat(3, img_pig(:,:,1:2), img_cdi, img_hdi);
    V=[0.4 0.6]; G=[1 1 2 2];
elseif strcmpi(featureID, 'F4V2')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_pig = get_features(img_uint8, 'Skin', 'RGB');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    img_hdi = get_features(img_uint8, 'Hue', 'HueProb');
    F = cat(3, img_lab(:,:,2), img_pig(:,:,1), img_cdi, img_hdi);
    V=[0.4 0.6]; G=[1 1 2 2];
elseif strcmpi(featureID, 'F4V3')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    F = cat(3, img_lab, img_cdi);
    V=[0.4 0.6]; G=[1 1 1 2];
elseif strcmpi(featureID, 'F4V4')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_pig = get_features(img_uint8, 'Skin', 'RGB');
    F = cat(3, img_lab(:,:,1:2), img_pig(:,:,1:2));
    V=[0.5 0.5]; G=[1 1 2 2];
elseif strcmpi(featureID, 'F4V5')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_pig = get_features(img_uint8, 'Skin', 'RGB');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    img_hdi = get_features(img_uint8, 'Hue', 'HueProb');
    F = cat(3, img_lab(:,:,2),img_pig(:,:,1), img_cdi, img_hdi);
    V = [0.2, 0.2, 0.2, 0.2]; G = [1 2 3 4];
elseif strcmpi(featureID, 'F4V6')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_pig = get_features(img_uint8, 'Skin', 'RGB');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    img_cv = get_features(img_uint8, 'ColorVolume', 'v1');
    F = cat(3, img_lab(:,:,2),img_pig(:,:,1), img_cdi, img_cv);
    V = [0.2, 0.2, 0.2, 0.2]; G = [1 2 3 4];
elseif strcmpi(featureID, 'F4V7')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_pig = get_features(img_uint8, 'Skin', 'RGB');
    img_hdi = get_features(img_uint8, 'Hue', 'HueProb');
    img_cv = get_features(img_uint8, 'ColorVolume', 'v1');
    F = cat(3, img_lab(:,:,2),img_pig(:,:,1), img_hdi, img_cv);
    V = [0.2, 0.2, 0.2, 0.2]; G = [1 2 3 4];
elseif strcmpi(featureID, 'F4V8')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    img_hdi = get_features(img_uint8, 'Hue', 'HueProb');
    img_cv = get_features(img_uint8, 'ColorVolume', 'v1');
    F = cat(3, img_lab(:,:,2),img_cdi, img_hdi, img_cv);
    V = [0.2, 0.2, 0.2, 0.2]; G = [1 2 3 4];
elseif strcmpi(featureID, 'F4V9')
    img_pig = get_features(img_uint8, 'Skin', 'RGB');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    img_hdi = get_features(img_uint8, 'Hue', 'HueProb');
    img_cv = get_features(img_uint8, 'ColorVolume', 'v1');
    F = cat(3, img_pig(:,:,1), img_cdi, img_hdi, img_cv);
    V = [0.2, 0.2, 0.2, 0.2]; G = [1 2 3 4];
elseif strcmpi(featureID, 'F5V1')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_pig = get_features(img_uint8, 'Skin', 'RGB');
    F = cat(3, img_lab, img_pig(:,:,1:2));
    V=[0.4 0.6]; G=[1 1 1 2 2];
elseif strcmpi(featureID, 'F5V2')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    img_hdi = get_features(img_uint8, 'Hue', 'HueProb');
    F = cat(3, img_lab, img_cdi, img_hdi);
    V=[0.4 0.6]; G=[1 1 1 2 2];
elseif strcmpi(featureID, 'F5V3')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_pig = get_features(img_uint8, 'Skin', 'RGB');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    img_hdi = get_features(img_uint8, 'Hue', 'HueProb');
    img_cv = get_features(img_uint8, 'ColorVolume', 'v2');
    F = cat(3, img_lab(:,:,2),img_pig(:,:,1), img_cdi, img_hdi, img_cv);
    V=[0.3 0.4, 0.3]; G=[1 1 2 2 3];
elseif strcmpi(featureID, 'F5V4')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_pig = get_features(img_uint8, 'Skin', 'RGB');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    img_hdi = get_features(img_uint8, 'Hue', 'HueProb');
    img_cv = get_features(img_uint8, 'ColorVolume', 'v1');
    F = cat(3, img_lab(:,:,2),img_pig(:,:,1), img_cdi, img_hdi, img_cv);
    V=[0.3 0.4, 0.3]; G=[1 1 2 2 3];
elseif strcmpi(featureID, 'F5V5')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_pig = get_features(img_uint8, 'Skin', 'RGB');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    img_hdi = get_features(img_uint8, 'Hue', 'HueProb');
    img_cv = get_features(img_uint8, 'ColorVolume', 'v1');
    F = cat(3, img_lab(:,:,2),img_pig(:,:,1), img_cdi, img_hdi, img_cv);
    V=[0.2, 0.2, 0.2, 0.2, 0.2]; G=[1 2 3 4 5];
elseif strcmpi(featureID, 'F5V6')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_pig = get_features(img_uint8, 'Skin', 'RGB');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    img_hdi = get_features(img_uint8, 'Hue', 'HueProb');
    img_cv = get_features(img_uint8, 'ColorVolume', 'v1');
    img_a = img_lab(:,:,2) / 50;
    img_cv = ext_normalize(img_cv, 'max');
    F = cat(3, img_a, img_pig(:,:,1), img_cdi, img_hdi, img_cv);
    V=[0.2, 0.2, 0.2, 0.2, 0.2]; G=[1 2 3 4 5];
elseif strcmpi(featureID, 'F5V7')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_pig = get_features(img_uint8, 'Skin', 'RGB');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    img_hdi = get_features(img_uint8, 'Hue', 'HueProb');
    img_cv = get_features(img_uint8, 'ColorVolume', 'v1');
    img_a = img_lab(:,:,2) / 50;
    img_cv = ext_normalize(img_cv, 'max');
    F = cat(3, img_a, img_pig(:,:,1), img_cdi, img_hdi, img_cv);
    V=[0.2, 0.2, 0.3, 0.1, 0.2]; G=[1 2 3 4 5];
elseif strcmpi(featureID, 'F5V8')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_pig = get_features(img_uint8, 'Skin', 'RGB');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    img_hsv = get_features(img_uint8, 'Color', 'HSV');
    img_cv = get_features(img_uint8, 'ColorVolume', 'v1');
    F = cat(3, img_lab(:,:,2),img_pig(:,:,1), img_cdi,img_hsv(:,:,1),img_cv);
    V=[0.2, 0.2, 0.2, 0.2, 0.2]; G=[1 2 3 4 5];
elseif strcmpi(featureID, 'F5V9')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_pig = get_features(img_uint8, 'Skin', 'RGB');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    img_hdi = get_features(img_uint8, 'Hue', 'HueProb');
    img_cv = get_features(img_uint8, 'ColorVolume', 'v2');
    F = cat(3, img_lab(:,:,2),img_pig(:,:,1), img_cdi, img_hdi, img_cv);
    V=[0.2, 0.2, 0.2, 0.2, 0.2]; G=[1 2 3 4 5];
elseif strcmpi(featureID, 'F5V10')
    img_a = get_features(img_uint8, 'ColorFit', 'AProb');
    img_pig = get_features(img_uint8, 'Skin', 'RGB');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    img_hdi = get_features(img_uint8, 'Hue', 'HueProb');
    img_cv = get_features(img_uint8, 'ColorVolume', 'v1');
    F = cat(3, img_a,img_pig(:,:,1), img_cdi, img_hdi, img_cv);
    V=[0.2, 0.2, 0.2, 0.2, 0.2]; G=[1 2 3 4 5];
elseif strcmpi(featureID, 'F6V1')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_pig = get_features(img_uint8, 'Skin', 'RGB');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    img_hdi = get_features(img_uint8, 'Hue', 'HueProb');
    % img_cv = get_features(img_uint8, 'ColorVolume', 'v2');
    F = cat(3, img_lab(:,:,1:2),img_pig(:,:,1:2), img_cdi, img_hdi);
    V=[0.3 0.3, 0.4]; G=[1 1 2 2 3 3];
elseif strcmpi(featureID, 'F6V2')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_pig = get_features(img_uint8, 'Skin', 'RGB');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    F = cat(3, img_lab, img_pig(:,:,1:2), img_cdi);
    V=[0.3 0.3 0.4]; G=[1 1 1 2 2 3];
elseif strcmpi(featureID, 'F6V3')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_pig = get_features(img_uint8, 'Skin', 'RGB');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    img_hdi = get_features(img_uint8, 'Hue', 'HueProb');
    F = cat(3, img_lab,img_pig(:,:,1), img_cdi, img_hdi);
    V=[0.3 0.3, 0.4]; G=[1 1 1 2 3 3];
elseif strcmpi(featureID, 'F6V4')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_pig = get_features(img_uint8, 'Skin', 'RGB');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    img_hdi = get_features(img_uint8, 'Hue', 'HueProb');
    img_cv = get_features(img_uint8, 'ColorVolume', 'v1');
    F = cat(3, img_lab(:,:,2),img_pig(:,:,1:2), img_cdi, img_hdi, img_cv);
    V=[1/6 1/6 1/6 1/6 1/6 1/6]; G=[1 2 3 4 5 6];
elseif strcmpi(featureID, 'F6V5')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_pig = get_features(img_uint8, 'Skin', 'RGB');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    img_hdi = get_features(img_uint8, 'Hue', 'HueProb');
    img_cv = get_features(img_uint8, 'ColorVolume', 'v1');
    F = cat(3, img_lab(:,:,1:2),img_pig(:,:,1), img_cdi, img_hdi, img_cv);
    V=[1/6 1/6 1/6 1/6 1/6 1/6]; G=[1 2 3 4 5 6];
elseif strcmpi(featureID, 'F7V1')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_pig = get_features(img_uint8, 'Skin', 'RGB');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    img_hdi = get_features(img_uint8, 'Hue', 'HueNorm');
    F = cat(3, img_lab, img_pig(:,:,1:2), img_cdi, img_hdi);
    V=[0.3 0.2 0.5]; G=[1 1 1 2 2 3 3];
elseif strcmpi(featureID, 'F7V2')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_pig = get_features(img_uint8, 'Skin', 'RGB');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    img_hdi = get_features(img_uint8, 'Hue', 'HueProb');
    F = cat(3, img_lab, img_pig(:,:,1:2), img_cdi, img_hdi);
    V=[0.3 0.2 0.5]; G=[1 1 1 2 2 3 3];
elseif strcmpi(featureID, 'F7V3')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_loc = get_features(img_uint8, 'Location', 'XY');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    img_hdi = get_features(img_uint8, 'Hue', 'HueProb');
    F = cat(3, img_lab, img_loc, img_cdi, img_hdi);
    V=[0.5 0.2 0.3]; G=[1 1 1 2 2 3 3];    
elseif strcmpi(featureID, 'F7V4')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_pig = get_features(img_uint8, 'Skin', 'RGB');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    img_hdi = get_features(img_uint8, 'Hue', 'HueProb');
    img_cv = get_features(img_uint8, 'ColorVolume', 'v1');
    F = cat(3, img_lab,img_pig(:,:,1), img_cdi, img_hdi, img_cv);
    V=[0.2, 0.2, 0.2, 0.2, 0.2]; G=[1 1 1 2 3 4 5];
elseif strcmpi(featureID, 'F7V5')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_rgb = get_features(img_uint8, 'Color', 'RGB');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    img_hdi = get_features(img_uint8, 'Hue', 'HueProb');
    img_cv = get_features(img_uint8, 'ColorVolume', 'v1');
    F = cat(3, img_lab(:,:,2),img_rgb, img_cdi, img_hdi, img_cv);
    V=[0.2, 0.2, 0.2, 0.2, 0.2]; G=[1 2 2 2 3 4 5];
elseif strcmpi(featureID, 'F7V6')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_pig1 = get_features(img_uint8, 'Skin', 'RGB');
    EI_rg = img_pig1(:,:,1);
    img_pig2 = get_features(img_uint8, 'Skin', 'LAB');
    EI_a = ext_normalize(img_pig2(:,:,2), 'range');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    img_hdi = get_features(img_uint8, 'Hue', 'HueNorm');
    F = cat(3, img_lab, EI_rg, EI_a, img_cdi, img_hdi);
    V=[0.3 0.2 0.5]; G=[1 1 1 2 2 3 3];
elseif strcmpi(featureID, 'F8V1')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_pig = get_features(img_uint8, 'Skin', 'RGB');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    img_hdi = get_features(img_uint8, 'Hue', 'HueProb');
    img_cv = get_features(img_uint8, 'ColorVolume', 'v2');
    F = cat(3, img_lab, img_pig(:,:,1:2), img_cdi, img_hdi, img_cv);
    V=[0.2 0.2 0.4 0.2]; G=[1 1 1 2 2 3 3 4];
elseif strcmpi(featureID, 'F8V2')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_pig = get_features(img_uint8, 'Skin', 'RGB');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    img_hdi = get_features(img_uint8, 'Hue', 'HueProb');
    img_gabor = get_features(img_uint8, 'Texture', 'Gabor');
    F = cat(3, img_lab, img_pig(:,:,1:2), img_cdi, img_hdi, img_gabor);
    V=[0.2 0.2 0.4 0.2]; G=[1 1 1 2 2 3 3 4];
elseif strcmpi(featureID, 'F8V3')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_pig = get_features(img_uint8, 'Skin', 'RGB');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    img_hdi = get_features(img_uint8, 'Hue', 'HueProb');
    img_glcm = get_features(img_uint8, 'Texture', 'GLCM');
    F = cat(3, img_lab, img_pig(:,:,1:2), img_cdi, img_hdi, img_glcm);
    V=[0.2 0.2 0.4 0.2]; G=[1 1 1 2 2 3 3 4];
elseif strcmpi(featureID, 'F9V1')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_pig = get_features(img_uint8, 'Skin', 'RGB');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    img_hdi = get_features(img_uint8, 'Hue', 'HueProb');
    img_gabor = get_features(img_uint8, 'Texture', 'Gabor');
    img_glcm = get_features(img_uint8, 'Texture', 'GLCM');
    F = cat(3, img_lab, img_pig(:,:,1:2), img_cdi, img_hdi,img_gabor, img_glcm);
    V=[0.3 0.2 0.4 0.1]; G=[1 1 1 2 2 3 3 4 4];
elseif strcmpi(featureID, 'F10V1')
    img_lab = get_features(img_uint8, 'Color', 'LAB');
    img_pig = get_features(img_uint8, 'Skin', 'RGB');
    img_cdi = get_features(img_uint8, 'CDI', 'CIE1931-1Sigma');
    img_hdi = get_features(img_uint8, 'Hue', 'HueProb');
    img_homo = get_features(img_uint8, 'Homogeneity', 'v2');
    F = cat(3, img_lab, img_pig(:,:,1:2), img_cdi, img_hdi,img_homo);
    V=[0.3 0.2 0.4 0.1]; G=[1 1 1 2 2 3 3 4 4 4];

end

% Debug Information
% =========================================================
if arg.Results.debug_mode == 2
    figure; 
    frows = length(V); fcols = ceil(length(G)/frows); 
    for i = 1:size(F, 3)
        subplot(frows,fcols,i); 
        f = F(:,:,i); imshow(f, []);
        xstr = sprintf('[%.2f, %.2f]', min(f(:)), max(f(:))); xlabel(xstr);
        xstr = sprintf('Group %d: %.2f', G(i), V(G(i))); title(xstr);
    end
end

end