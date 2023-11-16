% Unified access interface for image features.
% =========================================================
function img_out = get_features(img_uint8, group_type, group_key)

    available_groups = {'Color','Texture','Edge','Location','Homogeneity',...
        'Skin','CDI','Hue','ColorVolume','Saturation','ColorFit','NormPDF'};
    if nnz(strcmpi(available_groups, group_type)) == 0
        fprintf("Available 'group_type' are: \n"); disp(available_groups)
        error("Error in function '%s'.", 'getFeatures');
    end
    
    img_float = double(img_uint8)/255; group_type = lower(group_type);
    switch group_type
        case 'color'
            img_out = processColorFeature(img_float, group_key);
        case 'texture'
            img_out = processTextureFeature(img_uint8, group_key);
        case 'edge'
            img_out = processEdgeFeature(img_uint8, group_key);
        case 'location'
            img_out = processLocationFeature(img_uint8, group_key);
        case 'homogeneity'
            img_out = processLocalHomogeneity(img_float, group_key);
        case 'skin'
            img_out = processSkinIndex(img_float, group_key); 
        case 'cdi'
            img_out = processCDIFeature(img_float, group_key);
        case 'hue'
            img_out = processHueFeature(img_float, group_key);         
        case 'colorvolume'
            img_out = processColorVolume(img_float, group_key);
        case 'saturation'
            img_out = processSaturationFeature(img_float, group_key); 
        case 'colorfit'    
            img_out = processColorFittness(img_float, group_key);
        case 'normpdf'
            img_out = processNormpdfFeature(img_float, group_key);
    end

end

%% Core Function
% =========================================================
function img_out = processColorFeature(img_rgb, space)
% processColorFeature - Processing color features.

    % addpath('./Modules/Colors/')

    if nargin < 2, space = 'LAB'; end

    if isa(img_rgb, 'uint8')
        img_rgb = double(img_rgb) / 255.;
    end

    space = lower(space);
    switch space
        case 'hsy'
            img_out = rgb2hsy(img_rgb);
        case {'gray', 'xyz', 'lab', 'hsv', 'ycbcr', 'ntsc'}
            img_out = feval(strcat('rgb2', space), img_rgb);
        case {'lch', 'hsl', 'hsi', 'yuv', 'yiq', 'ydbdr', 'ypbpr'}
            img_out = colorspace(strcat('rgb->', space), img_rgb);
        otherwise % 'rgb'
            img_out = img_rgb;
    end
end

function img_out = processTextureFeature(img_uint8, texture)
% processTextureFeature - Processing texture features, ['Gabor', 'GLCM'].

    if nargin < 2, texture = 'GLCM'; end
    
    if strcmpi(texture, 'Gabor')
        img_out = processGaborFeature(img_uint8);
    else
        img_out = processGLCMFeature(img_uint8);
    end
        
end

function img_out = processEdgeFeature(img_uint8, grad_op)
% Calculate edge gradient amplitude.

    % addpath('./Functions/') 
    
    if nargin < 2
        grad_op = 'sobel'; % {'sobel', 'prewitt', 'roberts','laplacian', 'default'}
    end  
    img_gray = rgb2gray(img_uint8); 
    img_out = imgradient(img_gray, 'grad_op', grad_op);
    % img_out = ext_normalize(img_out, 'max');

end

function img_out = processLocationFeature(img_uint8, coord_type)
% Generate coordinate images.

    if nargin < 2
        coord_type = 'xy';
    end
    [numRows, numCols, ~] = size(img_uint8);
    [colsIndex, rowsIndex] = meshgrid(1:numCols, 1:numRows);
    
    if strcmpi(coord_type, 'xy')        
        img_out = cat(3, rowsIndex, colsIndex);
    end

end

function img_out = processLocalHomogeneity(img_float, version)
% Computing the local homogeneity.

    % addpath('./Functions/') 
    if nargin < 2
        version = 'v1';
    end
    
    if strcmpi(version, 'v1') % CGFFCM
        img_out = imhomogeneity(img_float, 'grad_space', 'gray',  ...
            'grad_op','default','homo_space','hsv','norm_level','image-wise');
    elseif strcmpi(version, 'v2')
        img_out = imhomogeneity(img_float, 'grad_space', 'homo', ...
            'grad_op','sobel','homo_space','hsv','norm_level', 'channel-wise');
    elseif strcmpi(version, 'v3')
        img_out = imhomogeneity(img_float, 'grad_space', 'homo', ...
            'grad_op','sobel','homo_space','lab','norm_level', 'channel-wise');
    else
        img_out = imhomogeneity(img_float, 'grad_space', 'homo', ...
            'grad_op','sobel','homo_space','rgb','norm_level', 'channel-wise');
    end

end

function img_out = processSkinIndex(img_float, space)
% Calculation of skin assessment indices EI, MI and ITA.
%
% Pigment Index: Erythema Index (EI) or Melanin Index (MI)
% ITA (Individual Typology Angle)
% B. C. K. Ly, E. B. Dyer, J. L. Feig, A. L. Chien, and S. Del Bino, 
% "Research Techniques Made Simple: Cutaneous Colorimetry: A Reliable 
% Technique for Objective Skin Color Measurement," 
% J. Invest. Dermatol., vol. 140, no. 1, pp. 3-12. e1, 2020.
% 
% space: ['RGB', 'LAB']

    if nargin < 2, space = 'RGB'; end
    
    img_lab = rgb2lab(img_float);
    ITA = atan((img_lab(:,:,1) - 50) ./ img_lab(:,:,3)) * 180/pi;

    if strcmpi(space, 'RGB')    
        Rr = img_float(:,:,1); Rr(Rr<0.1)=0.1;	% Rr(Rr>0.9)=0.9;
        Rg = img_float(:,:,2); Rg(Rg<0.1)=0.1;	% Rg(Rg>0.9)=0.9;  
        EI = log10(Rr./Rg); 
        MI = log10(1./Rr);  
    elseif  strcmpi(space, 'LAB')
        EI = img_lab(:,:,2)/100;
        MI = img_lab(:,:,1)/100;    
    end
    
    img_out = cat(3, EI, MI, ITA);
end

function img_out = processCDIFeature(img_float, ell_id)

    % addpath('./Modules/Colors/') 
    if nargin < 2, ell_id = 'CIE1931-1Sigma'; end
    
    img_out = gamutTransform(img_float,'ell_id',ell_id,'alpha',0.5);

end

function img_out = processHueFeature(img_float, hueType)
% hueType: {'Hue', 'HueNorm', 'HueProb'}

    % addpath('./Modules/Colors/') 
    if nargin < 2, hueType = 'Hue'; end
    
    [img_out, ~, ~] = rgb2hsy(img_float, 'cieType', 'CIE1931', ...
        'StartAngle', -60, 'ClockWise', true, 'hueType', hueType);
    
end

function img_out = processSaturationFeature(img_float, satType)
% hueType: {'Hue', 'HueNorm', 'HueProb'}

    % addpath('./Modules/Colors/') 
    if nargin < 2, satType = 'SatHSV'; end
    
    if strcmpi(satType, 'SatHSV')
        img_hsv = rgb2hsv(img_float);
        x = img_hsv(:,:,2);
        mu = 0.3834; sigma = 0.1003; 
    elseif strcmpi(satType, 'SatHSY')
        img_hsy = rgb2hsy(img_float);
        x = img_hsy(:,:,2);
        mu = 0.3242; sigma = 0.1129; 
    end
    img_out = exp(-(x-mu).^2/(2*sigma^2))/(sqrt(2*pi)*sigma); 
    img_out = img_out / 4;

end

function img_out = processColorVolume(img_float, version)
% Color Volume
% v1: Exploiting Color Volume and Color Difference for Salient Region Detection (TIP2019)
% v2: Color space volume and superpixel based leukocyte image segmentation (ITME2019)

    if nargin < 2, version = 'v1'; end 
    img_lab = rgb2lab(img_float);
    [L1, a1, b1] = imsplit(img_lab / 100);
    
    if strcmpi(version, 'v1')
        img_out = 4 / 3 * pi .* L1 .* a1 .* b1;
    elseif strcmpi(version, 'v2')
        [R, G, B] = imsplit(img_float);
        b2 = 2 * (a1 + B) - (R + G);
        img_out = 1 / 4 * pi .* L1 .* a1 .* b2;
    end

end

function img_out = processColorFittness(img_float, space)
    if nargin < 2, space = 'AProb'; end

    if strcmpi(space, 'AProb')
        img_lab = rgb2lab(img_float);
        mu = 18.2676; sigma_ = 5.6107; max_prob = 0.0710;
        x = img_lab(:,:,2);
        y = exp(-(x-mu).^2/(2*sigma_^2))/(sqrt(2*pi)*sigma_);
        img_out = y / max_prob;
    end

end

function img_out = processNormpdfFeature(img_float, space)

    if nargin < 2, space = 'LAB'; end

    if strcmpi(space, 'LAB')
        img_norm = rgb2lab(img_float);
        mu = [47.0634,18.2676,11.9498];
        % sigma_ = [134.0393, 0.0940, -0.2204; 0.0940,31.4803,7.8756; -0.2204,7.8756,28.7961];
        sigma_ = diag([134.0393,31.4803,28.7961]);
    elseif strcmpi(space, 'RGB')
        img_norm = img_float;
        mu = [0.5783, 0.3929, 0.3651];
        sigma_ = [0.0158,0.0129,0.0126;0.0129,0.0127,0.0124;0.0126,0.0124,0.0132];
    elseif strcmpi(space, 'XYZ')
        img_norm = rgb2xyz(img_float);
        mu = [0.2020, 0.1767, 0.1403];
        sigma_ = [0.0093,0.0087,0.0078;0.0087,0.0082,0.0075;0.0078,0.0075,0.0073];
    elseif strcmpi(space, 'Yxy')
        img_norm = rgb2Yxy(img_float, 'cieType', 'CIE1931');
        mu = [17.6745,0.3977,0.3418];
        sigma_ = [82.1648,-0.1427,-0.0125;-0.1427,0.0009,0.0002;-0.0125,0.0002,0.0002];
    elseif strcmpi(space, 'YCbCr')
        img_norm = rgb2ycbcr(img_float);
        mu = [0.4451,0.4623,0.5854];
        sigma_ = [0.0096,-0.0002,0.0004; -0.0002,0.0003,-0.0002;0.0004,-0.0002,0.0005];
    end
    arr = reshape(img_norm, [], 3);
    % p = ext_normpdf(arr, mu, sigma_);
    p = mvnpdf(arr, mu, sigma_); 
    img_out = reshape(p, [size(img_float,1), size(img_float,2)]);

end

% Auxiliary Functions
% =========================================================
function img_gabor = processGaborFeature(img_uint8)
% Gabor Filters

    img_gray = rgb2gray(img_uint8);
    
    % Design Array of Gabor Filters
    [numRows, numCols, ~] = size(img_uint8);
    
    wavelengthMin = 4/sqrt(2); 
    wavelengthMax = hypot(numRows,numCols);
    n = floor(log2(wavelengthMax/wavelengthMin));
    wavelength = 2.^(0:(n-2)) * wavelengthMin;
    
    deltaTheta = 45; 
    orientation = 0:deltaTheta:(180-deltaTheta);
    
    g = gabor(wavelength, orientation);
    
    % Extract Gabor magnitude features from source image
    gabormag = imgaborfilt(img_gray, g);
    
    % Post-process the Gabor Magnitude Images into Gabor Features.
    for i = 1:length(g)
        sigma = 0.5*g(i).Wavelength; K = 3;
        gabormag(:,:,i) = imgaussfilt(gabormag(:,:,i), K*sigma);
    end
    
    [X,Y] = meshgrid(1:numCols, 1:numRows);
    featureSet = cat(3, gabormag, X);
    featureSet = cat(3, featureSet, Y);
    
    X = reshape(featureSet, numRows*numCols, []);
    
    % Normalize the features to be zero mean, unit variance.
    X = bsxfun(@minus, X, mean(X));
    X = bsxfun(@rdivide, X, std(X));
    
    coeff = pca(X);
    
    img_gabor = reshape(X*coeff(:,1), numRows, numCols);
end

function img_glcm = processGLCMFeature(img_uint8)
% GLCM (Gray-Level Co-Occurrence Matrix)

    img_gray = rgb2gray(img_uint8);

    offsets = [0 1; -1 1; -1 0; -1 -1]; NumLevels = 9;
    [~, SI] = graycomatrix(img_gray, 'Offset', offsets, ...
        'Symmetric', true, 'NumLevels', NumLevels, 'GrayLimits', []);
    img_glcm = SI / NumLevels;
end

