function [vx,vy,colorflow] = computeDisparity(im1,im2)

im1=im2double(im1);
im2=im2double(im2);
%% compute disparity via SIFT flow
% cellsize=3;
% gridspacing=1;

% addpath(fullfile(pwd,'mexDenseSIFT'));
% addpath(fullfile(pwd,'mexDiscreteFlow'));

% sift1 = mexDenseSIFT(im1,cellsize,gridspacing);
% sift2 = mexDenseSIFT(im2,cellsize,gridspacing);
% 
% SIFTflowpara.alpha=2*255;
% SIFTflowpara.d=40*255;
% SIFTflowpara.gamma=0.005*255;
% SIFTflowpara.nlevels=4;
% SIFTflowpara.wsize=2;
% SIFTflowpara.topwsize=10;
% SIFTflowpara.nTopIterations = 60;
% SIFTflowpara.nIterations= 30;
% 
% 
% [vx,vy,energylist]=SIFTflowc2f(sift1,sift2,SIFTflowpara);

%% compute disparity via LDOF for better performance
flowframe = mex_LDOF(im1, im2);
vx = flowframe( :, :, 2 );
vy = flowframe( :, :, 1 );
vy(:)=0;
%warpI2=warpImage(im2,vx,vy);

%% compute disparity via SGM in 'RunStereo' field

%% display flow
clear flow;
flow(:,:,1)=vx;
flow(:,:,2)=vy;
colorflow = flowToColor(flow);