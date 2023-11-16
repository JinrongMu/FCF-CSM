% demo_supervoxels.m    The supervoxels demo.
% ======================================================================
%
% Input Arguments
% ======================================================================
% [1] Img   : 8 bit images (color or grayscale)
% [2] N0    : Number of required superpixels (optional, default is 200)
% [3] C     : Compactness factor (optional, default is 10)
% 
% Output Arguments
% ======================================================================
% [1] L     : Label matrix (in raster scan order)
% [2] N1    : Number of supervoxels computed
% 
% Notes:
% ======================================================================
% [1] The output numlabels gives number of superpixels generated.
% [2] number of returned supervoxels may be different from the input
% number of supervoxels.
% [3] you must compile the C file using:
% 
%     "mex slicsupervoxel_mex.c"
% 
% before using the Matlab code below.
%
% Example below shows how a stack may be prepared as input. It also
% shows how a desired supervoxel size may be used to obtain the
% required number of supervoxels.
% 
% You may want to try different compactness value as per you needs.
%======================================================================
% clear all;
% colorstack = 1;  % set to zero for a grayscale demo.
% im = imread('bee.jpg');
% reqdsupervoxelsize = 1000;
% if 1 == colorstack
%     stack = cat(4,im,im,im,im,im);
%     dims = size(stack);
%     numreqiredsupervoxels = dims(1)*dims(2)*dims(4)/reqdsupervoxelsize;
% else
%     im = rgb2gray(im);
%     stack = cat(3,im,im,im,im,im);
%     dims = size(stack);
%     numreqiredsupervoxels = prod(dims)/reqdsupervoxelsize;
% end
%
% compactness =  10;
% [labels, numlabels] = slicsupervoxel_mex(stack,numreqiredsupervoxels,compactness);
% figure,imagesc(labels(:,:,1));
%
% Start Testing
% ======================================================================
% 
clear 
close all
clc

test_single_slice()

test_multi_slice()


function test_single_slice()

% Load data: '*.nii.gz' or '*.nii'
filename = 'breast_slice.nii.gz';
[Img, ~] = read_nii_data(filename);

% Superpixels 
N0 = 100; C = 0.5; tic;
[L, N1] = slic_mex(Img,N0,C);
runtime = toc;

% Dispaly
figure, imagesc(L), 
xstr = sprintf("NumLabels = %d, RunTime = %.2fs", N1, runtime);
xlabel(xstr);

figure;
BW = boundarymask(L);
imshow(imoverlay(Img,BW,'cyan'),'InitialMagnification','fit')
xstr = sprintf("NumLabels = %d, RunTime = %.2fs", N1, runtime);
xlabel(xstr);

end

function test_multi_slice()

% Setup
test_method = 'slicsupervoxel';  % {'slicsupervoxel', 'superpixels3'}
debug_mode = 2;
% 1 : save nii label
% 2 : imshow + boundarymask
% 3 : imagesc + label
% 4 : volumeViewer + 3D label
% 5 : implay +  boundarymask
% 6 : implay + color label

% Load data: '*.nii.gz' or '*.nii'
filename = 'breast.nii.gz';
[Img, filename] = read_nii_data(filename);

% Superpixels 
supervoxelsize = 500; C = 0.5; 
% [m,n,d] = size(Img); N0 = round(m * n * d / supervoxelsize); 
% N0 = prod(size(Img)) / supervoxelsize;
N0 = round(numel(Img) / supervoxelsize);

tic;
if strcmpi(test_method, 'slicsupervoxel')
    [L, N1] = slicsupervoxel_mex(Img, N0, C);
else
    [L, N1] = superpixels3(Img, N0, 'Method', 'slic', 'Compactness', 0.5);
end
runtime = toc;

% Dispaly
xstr = sprintf("NumLabels = %d, RunTime = %.2fs", N1, runtime);
disp(xstr);

if debug_mode == 1
    save_nii_label(filename);
    
elseif debug_mode == 2
    for i = 1:size(Img,3)
        figure(2)
        BW = boundarymask(L(:,:,i));
        imshow(imoverlay(Img(:,:,i),BW,'cyan'),'InitialMagnification','fit'),
        title(sprintf('slice %d', i))
        pause(0.1)
    end
    
elseif debug_mode == 3
    for i = 1:size(Img,3)
        figure(3)
        imagesc(L(:,:,i));
        title(sprintf('slice %d', i))
        pause(0.1)
    end
    
elseif debug_mode == 4
    volumeViewer(L)
    
elseif debug_mode == 5
    imSize = size(Img); conn = 4;
    imPlusBoundaries = zeros(imSize(1),imSize(2),3,imSize(3),'uint8');
    % Create an RGB representation of this plane with boundary shown in cyan.
    for plane = 1:imSize(3)
        BW = boundarymask(L(:, :, plane), conn);
        imPlusBoundaries(:, :, :, plane) = imoverlay(Img(:, :, plane), BW, 'cyan');
    end
    implay(imPlusBoundaries,5)
    
elseif debug_mode == 6
    L = L + 1;
    pixelIdxList = label2idx(L);
    meanA = zeros(size(Img),'like',Img);
    for superpixel = 1:N1
        memberPixelIdx = pixelIdxList{superpixel};
        meanA(memberPixelIdx) = mean(Img(memberPixelIdx));
    end
    implay([Img meanA],5);
    
end 

end

function [data, filename] = read_nii_data(filename)
%Read nii data from '*.nii.gz' or '*.nii'.

    % filename = 'breast_slice.nii.gz';
    [~, ~, ext] = fileparts(filename);
    if strcmpi(ext, '.gz')
        files = gunzip(filename); 
        filename = files{1};
    end
    data = niftiread(filename); 

end

function save_nii_label(filename, include_info)
%Save labels to nii data '*_label.nii'.

    % filename = 'breast.nii';
    if ~exist('include_info', 'var')
        include_info = true;
    end
    
    info = niftiinfo(filename);
    info.Filename = strrep(info.Filename, '.nii', '_label.nii');
    info.Datatype = 'int32';
    if include_info 
        niftiwrite(L,info.Filename,info);
    else
        niftiwrite(L,info.Filename);
    end

end
