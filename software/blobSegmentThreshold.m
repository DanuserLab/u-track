function maskBlobs = blobSegmentThreshold(image,minSize,plotRes,mask)
%BLOBSEGMENTTHRESHOLD segments blobs in 2D images via combination of Otsu and Rosin thresholding
%
%SYNOPSIS maskBlobs = blobSegmentThreshold(image,minSize,plotRes)
%
%INPUT  image     : 2D image to be segmented.
%       minSize   : Minimum size of a blob. 
%                   Optional. Default: 20 pixels.
%       plotRes   : 1 to plot segmentation results, 0 otherwise.
%                   Optional. Default: 0.
%       mask      : Binary mask. Optional. If not provided, the whole image
%                   domain is segmented.
%
%OUTPUT maskBlobs : Mask of blobs. 1 inside blobs, 0 outside.
%
%REMARKS While the code is in principle general, it has been extensively tested only on focal adhesions.
%
%Khuloud Jaqaman August 2009
%
% Copyright (C) 2018, Danuser Lab - UTSouthwestern 
%
% This file is part of u-track.
% 
% u-track is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% u-track is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with u-track.  If not, see <http://www.gnu.org/licenses/>.
% 
% 

%% Output
maskBlobs = [];

%% Input

%check number of input arguments
if nargin < 1
    disp('Please enter at least image to be segmented');
    return
end

%minimum blob size
if nargin < 2 || isempty(minSize)
    minSize = 20;
end

%plot results
if nargin < 3 || isempty(plotRes)
    plotRes = 0;
end

%mask
if nargin < 4 || isempty(mask)
    mask = ones(size(image));
end

if ~logical(mask)
    error('Mask must be a logical image.');
end
    
%% Segmentation

%make sure that image is in double format
image = double(image);

%remove noise by filtering image with a Gaussian whose sigma = 1 pixel
imageFiltered = filterGauss2D(image,1);

%estimate background by filtering image with a Gaussian whose sigma = 10 pixels
imageBackground = filterGauss2D(image,10);

%calculate noise-filtered and background-subtracted image
imageFilteredMinusBackground = imageFiltered - imageBackground;

%crop image
imageFilteredMinusBackground = imageFilteredMinusBackground .* mask;

%enhance features by performing a maximum filter
[sizeX,sizeY] = size(imageFilteredMinusBackground);
imageDilated = imageFilteredMinusBackground;
imageTmp(:,:,1) = imageDilated;
imageTmp(:,:,2) = [zeros(1,sizeY); imageDilated(2:end,:)];
imageTmp(:,:,3) = [imageDilated(1:end-1,:); zeros(1,sizeY)];
imageTmp(:,:,4) = [zeros(sizeX,1) imageDilated(:,2:end)];
imageTmp(:,:,5) = [imageDilated(:,1:end-1) zeros(sizeX,1)];
imageTmp(:,:,6) = [zeros(1,sizeY); [zeros(sizeX-1,1) imageDilated(2:end,2:end)]];
imageTmp(:,:,7) = [zeros(1,sizeY); [imageDilated(2:end,1:end-1) zeros(sizeX-1,1)]];
imageTmp(:,:,8) = [[zeros(sizeX-1,1) imageDilated(1:end-1,2:end)]; zeros(1,sizeY)];
imageTmp(:,:,9) = [[imageDilated(1:end-1,1:end-1) zeros(sizeX-1,1)]; zeros(1,sizeY)];
imageDilated = max(imageTmp,[],3);

% Find non zero values (due to masking)
nzInd = find(imageDilated);

%get minumum and maximum pixel values in image
minSignal = min(imageDilated(nzInd));
maxSignal = max(imageDilated(nzInd));

%normalize non zero value between 0 and 1
imageDilatedNorm = zeros(size(imageDilated));
imageDilatedNorm(nzInd) = (imageDilated(nzInd) - minSignal) / (maxSignal - minSignal);

%estimate the intensity level to use for thresholding the image
level1 = graythresh(imageDilatedNorm(nzInd)); %Otsu
[dummy, level2] = cutFirstHistMode(imageDilatedNorm(nzInd),0); %Rosin
level = 0.33333*level1 + 0.66667*level2;

%threshold the image
imageThresholded = im2bw(imageDilatedNorm,level);

%fill holes in thresholded image to make continuous blobs
imageThresholdedFilled = imfill(imageThresholded,'holes');

% go over blobs and remove those with a size smaller that minSize
labels = bwlabel(imageThresholdedFilled);
stats = regionprops(labels, 'Area'); %#ok<MRPBW>
idx = find([stats.Area] > minSize);

%output final blob mask
maskBlobs = ismember(labels, idx);

%% Plotting

if plotRes

    %get the blob edges from the final blob mask
    SE = strel('square',3);
    edgesBlobs = imdilate(maskBlobs,SE) - maskBlobs;

    %scale the original image to be between 0 and 1
    %actually scale it between the 1st and 99th percentiles
    imageScaled = (image - prctile(image(:),1)) / (prctile(image(:),99) - prctile(image(:),1));
    imageScaled(imageScaled<0) = 0;
    imageScaled(imageScaled>1) = 1;

    %plot image
    figure, hold on
    subplot(1,2,1)
    imshow(imageScaled,[])

    %give the edge pixels a value of zero in the original image
    imageScaled(edgesBlobs==1) = 0;

    %construct a 3-layered image to show blob edges on top of
    %original image
    image3Color = repmat(imageScaled,[1 1 3]);
    image3Color(:,:,1) = image3Color(:,:,1) + edgesBlobs;
    
    %plot mask edges
    subplot(1,2,2)
    imshow(image3Color,[]);

end

%% extra stuff

%otherwise, if smaller than the intermediate size ...
%     elseif sizeFACurrent < intSize
%
%         %calculate actual intensity factor (ratio of FA intensity to
%         %background)
%         intensityFactor = intensityFACurrent / intensityBG;
%
%         %discard if intensity factor is smaller than minimum allowed for
%         %this size
%         if intensityFactor < minSizeMinIntensityFactor
%             indxFA(indxFACurrent) = 0;
%         end
%
%         %otherwise, if smaller than the maximum size with extra intensity
%         %penalty
%     elseif sizeFACurrent < maxSizeIntensityPenalty
%
%         %calculate minimum allowed intensity factor
%         minIntensityFactor = yIntersept + lineSlope * sizeFACurrent;
%
%         %calculate actual intensity factor (ratio of FA intensity to
%         %background)
%         intensityFactor = intensityFACurrent / intensityBG;
%
%         %discard if intensity factor is smaller than minimum allowed for
%         %this size
%         if intensityFactor < minIntensityFactor
%             indxFA(indxFACurrent) = 0;
%         end
%
%         %otherise ...
%     else
%
%         %calculate actual intensity factor (ratio of FA intensity to
%         %background)
%         intensityFactor = intensityFACurrent / intensityBG;
%
%         %discard if intensity factor is smaller than minimum allowed
%         %overall
%         if intensityFactor < minIntensityFactorAll
%             indxFA(indxFACurrent) = 0;
%         end



%extract focal adhesion selection criteria
% minSize = selectParam.minSize;
% minSizeMinIntensityFactor = selectParam.minSizeMinIntensityFactor;
% maxSizeIntensityPenalty = selectParam.maxSizeIntensityPenalty;
% minIntensityFactorAll = selectParam.minIntensityFactorAll;
%
% %calculate straight line equation for focal adhesion selection
% intSize = minSize + (maxSizeIntensityPenalty - minSize) / 3;
% % intSize = minSize;
% x1 = intSize;
% y1 = minSizeMinIntensityFactor;
% x2 = maxSizeIntensityPenalty;
% y2 = minIntensityFactorAll;
% lineSlope = (y1 - y2) / (x1 - x2);
% yIntersept = y1 - lineSlope * x1;

