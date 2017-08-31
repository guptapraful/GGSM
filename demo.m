clear all
clc
close all

% P. Gupta , A. K. Moorthy, R. Soundararajan and A. C. Bovik, "DIIVINE-GGSM Software Release", 
% URL: http://live.ece.utexas.edu/research/quality/diivine-ggsm.zip, 2017
% 
% Author  : Praful Gupta
% Version : 1.1
% 
% The authors are with the Laboratory for Image and Video Engineering
% (LIVE), Department of Electrical and Computer Engineering, The
% University of Texas at Austin, Austin, TX.
% 
% Kindly report any suggestions or corrections to praful_gupta@utexas.edu

im_ref = double(rgb2gray(imread('test_images/bikes.bmp')));
quality_ref = ggsm_divine(im_ref);

im_jpeg = double(rgb2gray(imread('test_images/bikes_jpeg.bmp')));
quality_jpeg = ggsm_divine(im_jpeg);