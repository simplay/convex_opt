% Convex Optimization - Project 1
% MICHAEL SINGLE
% 08-917-445

%% init data
clear all;
close all;
clc;

iter = 10000;
FIND_OPTIMUM = false;
im = imread('Input/fruits.png');
im = imresize(im,2);
[m,n,~] = size(im);
im = im2double(im);

%% compute demosaicing input

% bayer filter tensor
[ red_mask, green_mask, blue_mask ] = getMasks( m, n );
Omega = mat2Img(red_mask, green_mask, blue_mask);

% demosaiced img
r = im(:,:,1).*red_mask;
g = im(:,:,2).*green_mask;
b = im(:,:,3).*blue_mask;
mosaiced = mat2Img(r,g,b);

%% find best lambda
if FIND_OPTIMUM
    disp('starting finding best lambda');
    [ bestLambda ] = findBestLambda(im, mosaiced, Omega, 100 );
    disp('determined best lambda value');
else
    % determined using findBestLambda(...)
    % for Input/fruits.png at 200x266 pixels
    % bestLambda = 1899;  % 1740 for big
    bestLambda = 1740;
end

%% compute demosaiced img
demosaicedImg = demosaic(mosaiced, Omega, bestLambda, iter);

%% display results
figure;
disp = [mosaiced, (im-demosaicedImg).^2; demosaicedImg, im];
imshow(disp);