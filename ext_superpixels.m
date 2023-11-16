function varargout = ext_superpixels(Img, N0, varargin)
% ext_superpixels - The ensemble of superpixel algorithms.
%
% Syntax
% =================
% L = ext_superpixels(img_uint8, N0);
% [L, N1] = ext_superpixels(img_uint8, N0, 'method', 'gGMMSP');
%
% Input Arguments
% =================
% Img           8 bit images (color or grayscale).
% N0            Desired number of superpixels
%
% method        Superpixel algorithm, {}
% debug_mode    Print debug information, 0: Silent, 1: Call information, 2: Call details.
%
% Output Arguments
% =================
% out           L or {L, N1}
% L             Labels
% N1            Number of labels (final output superpixels) 
%
% Reference Papers
% =================
% [1] M. Y. Liu, O. Tuzel, S. Ramalingam, and R. Chellappa, 
% "Entropy rate superpixel segmentation." CVPR2011, pp.2097-2104.('ERS') 
% [2] R. Achanta, A. Shaji, K. Smith, A. Lucchi, P. Fua, and S. Susstrunk,
% "SLIC superpixels compared to state-of-the-art superpixel methods," ('Matab','SLIC','SLICO',)
% IEEE Trans. Pattern Anal. Mach. Intell., vol.34, no.11, pp.2274-2282,2012.
% [3] J. Shen, Y. Du, W. Wang, and X. Li,
% "Lazy random walks for superpixel segmentation," ('LRW')
% IEEE Trans. Image Process., vol. 23, no. 4, pp. 1451-1462, 2014.
% [4] Z. Ban, J. Liu, and L. Cao,
% "Superpixel segmentation using Gaussian mixture model," ('GMMSP')
% IEEE Trans. Image Process., vol. 27, no. 8, pp. 4105-4117, 2018.
% [5] Z. Ban, J. Liu, and J. Fouriaux, "GMMSP on GPU," (*'gGMMSP'*)
% Journal of Real-Time Image Processing, vol. 17, pp. 245-257, 2020.
% [6] R. Giraud, V.-T. Ta, and N. Papadakis,
% "Robust superpixels using color and contour features along linear path,"('SCALP')
% Comput. Vis. Image Und., vol. 170, pp. 1-13, 2018.
% [7] C. Wu, J. Zheng, Z. Feng, H. Zhang, L. Zhang, J. Cao, and H. Yan,
% "Fuzzy SLIC: Fuzzy simple linear iterative clustering," £¨'FuzzySLIC', 'FuzzySLICNC'£©
% IEEE Trans. Circuits Syst. Video Technol., vol. 31, no. 6, pp. 2114-2124, 2020.
% [8] R. Achanta, and S. Susstrunk,
% "Superpixels and polygons using simple non-iterative clustering," (SNIC)
% 2017 IEEE Conference on Computer Vision and Pattern Recognition (CVPR). 
% pp.4651-4660, 2017.
%
% TODO: LRW
%
% Parameter Initialization    
% =========================================================
% addpath('./Modules/Superpixels/')
% addpath('./Functions/')      % Required for the function 'process_labels'.

if ~exist('N0','var'),  N0=100; end

arg = inputParser; fun_name = 'ext_superpixels';
addParameter(arg,'method','default');
addParameter(arg,'debug_mode',0);
parse(arg,varargin{:});

% Method Implementation    
% =========================================================
method = lower(arg.Results.method); 
switch method
    
    case 'ers'
        % Syntax: 
        % L = ers_mex(Img,nC);
        % L = ers_mex(Img,nC,lambda,sigma,conn8);
        % ------------------------------------------------------------
        % Img    : the input grey scale image or color image (double)
        % nC     : the number of desired superpixels
        % lambda : the balancing parameter (default=0.5)
        % sigma  : the kernel bandwidth (default=5.0)
        % conn8  : the flag of using 8-connected grid graph structure (default=1)
        %          if set to 0, the algorithm will use 4-connected graph instead.
        lambda = 0.5; sigma = 5.0; conn8 = 1;
        L = ers_mex(double(Img), N0, lambda, sigma, conn8); 
    
    case 'slic'
        % Syntax: [L, N1] = slic_mex(Img, N0, C);
        % ------------------------------------------------------------
        % Img   : 8 bit images (color or grayscale)	
        % N0    : Number of required superpixels (optional, default is 200)
        % C     : Compactness factor (optional, default is 10)
        C = 10;
        [L, ~] = slic_mex(Img, N0, C); 
        
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
        [L, ~] = slico_mex(Img, N0); 
        
%     case 'lrw' 
%         % The execution is slow and the recommended image size is no more
%         % than 200 pixels.
%         [L, N1] = clusterLRW(Img, N0);
        
    case 'snic'
        % Syntax: [L, N1] = snic_mex(Img, N0, C); 
        % ------------------------------------------------------------
        % Img   : 8 bit images (color or grayscale)	
        % N0    : Number of required superpixels (optional, default is 200)
        % C     : Compactness factor, [10,40]
        C = 20.0;
        [L, ~] = snic_mex(Img, N0, C); 
        
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
        L = scalp_mex(Img, N0, C); 
        
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
        
    case 'fuzzyslic'
        % Syntax: [L, N1] = fuzzyslic_mex(Img, N0, C, FC, AC);
        % ------------------------------------------------------------
        % Img   : 8 bit images (color or grayscale)
        % N0    : Number of intended superpixels (optional, default is 200)
        % C     : Compactness coefficient (optional, default is 10), used to control the color sensitive.     
        % FC    : Optional input (Use 1 to select Fuzzy SLICNC, use 0 to select Fuzzy SLIC, default is 0)
        % AC    : Optional input (Amplification coefficient, default is 0.2)
        C = 10; FC = 0; AC = 0.2;
        [L, ~] = fuzzyslic_mex(Img, N0, C, FC, AC); 
        
    case 'fuzzyslicnc'
        C = 10; FC = 1; AC = 0.2;
        [L, ~] = fuzzyslic_mex(Img, N0, C, FC, AC); 
        
    otherwise  % 'matlab'
        warning("Call Matlab Function: %s.\n", 'superpixels')
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
        [L, ~] = superpixels(Img, N0, 'Compactness', 10);
        
end

[L, N1] = process_labels(L);

% Output Settings
% =========================================================
if nargout == 2
    varargout = {L, N1};
else
    varargout = {L};
end
   
% Debug Information
% =========================================================
if arg.Results.debug_mode == 1
    fprintf('\nCall Function:\t%s\n', fun_name)
    fprintf('----------------------------------------');
    fprintf('\nInput Parameters:\n');
    fprintf("\tmethod: '%s'\n\tN0: %d\n", arg.Results.method, N0);
    fprintf("Generate %d superpixels, maximum label is %d.\n", N1, max(L(:)));
   
elseif arg.Results.debug_mode == 2
    figure
    BW = boundarymask(L);
    imshow(imoverlay(Img,BW,'cyan'),'InitialMagnification','fit')
    xstr = sprintf("N0 = %d, N1 = %d, maxLabel = %d", N0, N1, max(L(:)));
    xlabel(xstr); title(arg.Results.method);

end

end