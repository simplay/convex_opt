% Convex Optimization - Project 2
% MICHAEL SINGLE
% 08-917-445

%% init data
clear all;
close all;
clc;

iter = 450; % 400;
FIND_OPTIMUM = false;
verbose = true;
im = imread('Input/fruits.png');
%im = imresize(im,2);
%im = imread('Input/woodensweater.jpg');
%im = imresize(im,0.5);
[m,n,~] = size(im);
im = im2double(im);
input = im;
%% compute demosaicing input

% bayer filter tensor
[ red_mask, green_mask, blue_mask ] = getMasks( m, n );

%red_mask = ones(size(red_mask));
%green_mask = ones(size(red_mask));
%blue_mask = ones(size(red_mask));


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
    % for Input/fruits.png at 100x133 pixels
    bestLambda = 1; %default value
    bestLambda = 2105; % found by running find best lambda script
end

%% compute demosaiced img
demosaicedImg = demosaicing_michael_single(mosaiced, Omega, bestLambda, iter, input, verbose);


%% display results
figure;
displ = [mosaiced, (im-demosaicedImg).^2; demosaicedImg, im];
imshow(displ);
disp('Image Legend From from left to the right from top to bottom:');
disp('mosaiced, ssd ground truth demosaiced, demosaiced, ground-truth');