% demo_superpixels.m    The superpixel demo.
% ======================================================================
%
% Input Arguments
% ======================================================================
% [1] Img: Input image
% [2] N0: Desired number of superpixels
% [3] C: Compactness factor
% 
% Output Arguments
% ======================================================================
% [1] L: Label matrix
% [2] N1: Number of superpixels computed
% 
% Notes:
% ======================================================================
% [1] number of returned superpixels may be different from the input 
% number of superpixels.
% [2] you must compile the C file using 'mex *mex.c' before using the 
% code below.
%
% Start Testing
% ======================================================================
clear 
close all
clc

availabel_methods = {'matlab', 'fuzzyslic','fuzzyslicnc','slic', 'slic0',...
    'snic', 'ers', 'scalp', 'gmmsp', 'ggmmsp'};

test_method =  'gmmsp'; 

Img = imread('bee.jpg'); N0 = 500; tic;
switch test_method
    
    case 'matlab'
        % Syntax: [L, N1] = superpixels(Img, N0, Name, Value, ...)
        % ------------------------------------------------------------
        % Img           : single | double | int16 | uint8 | uint16
        % N0            : Desired number of superpixels
        % Compactness   : Shape of superpixels, 10 (default)
        % IsInputLab    : Input image is in the L*a*b* colorspace, 
        %                 false (default) | true
        % Method        : Algorithm used to compute superpixels, 
        %                 'slic0' (default) | 'slic' 
        % NumIterations : Number of iterations used in the clustering phase 
        %                 of the algorithm, 10 (default) | numeric scalar
        [L, N1] = superpixels(Img, N0, 'Compactness', 10);
        
    case 'slic'
        % Syntax: [L, N1] = slic_mex(Img, N0, C);
        % ------------------------------------------------------------
        % Img   : 8 bit images (color or grayscale)	
        % N0    : Number of required superpixels (optional, default is 200)
        % C     : Compactness factor (optional, default is 10)
        C = 10;
        [L, N1] = slic_mex(Img, N0, C); 
        
    case 'slic0'
        % Syntax: [L, N1] = slico_mex(Img, N0, C);
        % ------------------------------------------------------------
        % How is SLICO different from SLIC?
        % ------------------------------------------------------------
        % 1. SLICO does not need compactness factor as input. It is 
        % calculated automatically
        % 2. The automatic value adapts to the content of the superpixel. 
        % So,SLICO is better suited for texture and non-texture regions
        % 3. The advantages 1 and 2 come at the cost of slightly poor 
        % boundary adherences to regions.
        % 4. This is also a very small computational overhead 
        % (but speed remains almost as fast as SLIC.
        % 5. There is a small memory overhead too w.r.t. SLIC.
        % 6. Overall, the advantages are likely to outweigh the small 
        % fdisadvantages or most applications of superpixels.
        [L, N1] = slico_mex(Img, N0); 
        
    case 'fuzzyslic'
        % Syntax: [L, N1] = fuzzyslic_mex(Img, N0, C, FC, AC);
        % ------------------------------------------------------------
        % Img   : 8 bit images (color or grayscale)
        % N0    : Number of intended superpixels (optional, default is 200)
        % C     : Compactness coefficient (optional, default is 10), used to control the color sensitive.     
        % FC    : Optional input (Use 1 to select Fuzzy SLICNC, use 0 to select Fuzzy SLIC, default is 0)
        % AC    : Optional input (Amplification coefficient, default is 0.2)
        C = 10; FC = 0; AC = 0.2;
        [L, N1] = fuzzyslic_mex(Img, N0, C, FC, AC); 
        
    case 'fuzzyslicnc'
        C = 10; FC = 1; AC = 0.2;
        [L, N1] = fuzzyslic_mex(Img, N0, C, FC, AC); 
        
    case 'snic'
        % Syntax: [L, N1] = snic_mex(Img, N0, C); 
        % ------------------------------------------------------------
        % Img   : 8 bit images (color or grayscale)	
        % N0    : Number of required superpixels (optional, default is 200)
        % C     : Compactness factor, [10,40]
        C = 20.0;
        [L, N1] = snic_mex(Img, N0, C); 
        
    case 'ers'
        % Syntax: [Labels] = ers_mex(Img,nC,lambda,sigma,conn8);
        % ------------------------------------------------------------
        % Img    : the input grey scale image or color image (double)
        % nC     : the number of desired superpixels
        % lambda : the balancing parameter (default=0.5)
        % sigma  : the kernel bandwidth (default=5.0)
        % conn8  : the flag of using 8-connected grid graph structure (default=1)
        %          if set to 0, the algorithm will use 4-connected graph instead.
        lambda = 0.5; sigma = 5.0; conn8 = 1;
        L = ers_mex(double(Img), N0, lambda, sigma, conn8); 
        N1 = length(unique(L));
        
    case 'scalp'
        % Syntax:
        % L = scalp_mex(uint8(img), K, C, single(NC));
        % L = scalp_mex(img, K, C);      % without contour prior
        % ------------------------------------------------------------
        % img	: CIE RGB color image, uint8
        % K     : Superpixel number
        % C     : Compactness parameter (default 0.075)
        % NC = double(imread('path_img_contour.png'));
        % NC = NC/max(NC(:));
        C = 0.075;  
        [L] = scalp_mex(Img, N0, C); 
        N1 = length(unique(L));
        
    case 'gmmsp'
        % Syntax:
        % L = mx_GMMSP(image, v_x, v_y);
        % L = mx_GMMSP(image, v_x, v_y, e_c, T, lambda, e_s);
        % ------------------------------------------------------------
        % v_x   : int,[8, +inf] 
        % v_y   : int,[8, +inf] 
        % e_c 	: (0, +inf) (default 8)
        % T     : int,[10, +inf] ,(default 10)
        % lambda: [8, +inf](default 8)
        % e_s   : (0, +inf) (default 2)
        v = max(8, round(sqrt(size(Img,1)*size(Img,2)/N0)));
        v_x = v; v_y = v; 
        e_c = 8; T = 10; lambda = 8; e_s = 2; 
        L = mx_GMMSP(Img, v_x, v_y, e_c, T, lambda, e_s);
        N1 = length(unique(L));
        
    case 'ggmmsp'
        % Syntax:
        % L = mx_gGMMSP(image, v_x, v_y);
        % L = mx_gGMMSP(image, v_x, v_y, e_c, T, lambda, e_s);
        % ------------------------------------------------------------
        % v_x   : int,[8, +inf] 
        % v_y   : int,[8, +inf] 
        % e_c 	: (0, +inf) (default 8)
        % T     : int,[10, +inf] ,(default 10)
        % lambda: [8, +inf](default 8)
        % e_s   : (0, +inf) (default 2)
        v = max(8, round(sqrt(size(Img,1)*size(Img,2)/N0)));
        v_x = v; v_y = v; 
        e_c = 8; T = 10; lambda = 8; e_s = 2; 
        L = mx_gGMMSP(Img, v_x, v_y, e_c, T, lambda, e_s);
        N1 = length(unique(L));
        
    otherwise
        error("Undefined Superpixel 'method': '%s'.\n", test_method)


end
runtime = toc;

% Dispaly
% ======================================================================
figure;
imagesc(L); colorbar; 
xstr = sprintf("RunTime = %.2fs", runtime);
xlabel(xstr);

figure;
BW = boundarymask(L);
imshow(imoverlay(Img,BW,'cyan'),'InitialMagnification','fit')
xstr = sprintf("NumLabels = %d, MaxLabels = %d", N1, max(L(:)));
xlabel(xstr); title(test_method)
