function [detectedFeatures,imageN3,errFlag] = centroidSubResFeatures2D(...
    image,cands,psfSigma,visual,bitDepth,saveResults)
%CENTROIDSUBRESFEATURES2D determines the positions and intensity amplitudes of sub-resolution features using centroid calculation
%
%SYNOPSIS [detectedFeatures,imageN3,errFlag] = centroidSubResFeatures2D(...
%    image,cands,psfSigma,visual,bitDepth,saveResults)
%
%INPUT  image      : Image being analyzed.
%       cands      : Cands structure as output from fsmCenter.
%       psfSigma   : Standard deviation of point spread function (in pixels).
%       visual     : 1 if user wants to view results; 0 otherwise.
%                    Optional. Default: 0.
%       bitDepth   : Camera bit depth. Optional. Default: 14.
%       saveResults: 1 if results are to be saved (in file 'detectedFeatures.mat'),
%                    0 otherwise. Optional. Default: 0.
%       All optional variables can be entered as [] to use default values.
%
%OUTPUT detectedFeatures: Structure with fields:
%             .xCoord    : Image coordinate system x-coordinate of detected
%                          features [x dx] (in pixels).
%             .yCoord    : Image coorsinate system y-coordinate of detected
%                          features [y dy] (in pixels).
%             .amp       : Amplitudes of PSFs fitting detected features [a da].
%       imageN3    : Image with labeled features. Blue: those from cands;
%                    Red: those from mixture-model fitting; Magenta: those
%                    from MMF which coincide with those from cands.
%                    Will be output only if visual = 1.
%       errFlag    : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, July 2011
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

detectedFeatures = [];
imageN3 = [];
errFlag = 0;

%% Input

%check whether correct number of input arguments was used
if nargin < 3
    disp('--centroidSubResFeatures2D: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%check visualization option
if nargin < 4 || isempty(visual)
    visual = 0;
else
    if visual ~= 0 && visual ~= 1
        disp('--centroidSubResFeatures2D: Variable "visual" should be 0 or 1!');
        errFlag = 1;
    end
end

%check the bit depth
if nargin < 5 || isempty(bitDepth)
    bitDepth = 14;
else
    if bitDepth <= 0 || bitDepth-floor(bitDepth) ~= 0
        disp('--centroidSubResFeatures2D: Variable "bitDepth" should be a positive integer!');
    end
end

%check whether results are to be saved
if nargin < 6 || isempty(saveResults)
    saveResults = 0;
else
    if saveResults ~= 0 && saveResults ~= 1
        disp('--detectSubResFeatures2D: Variable "saveResults" should be 0 or 1!');
    end
end

%exit if there are problems with input data
if errFlag
    disp('--centroidSubResFeatures2D: Please fix input data!');
    return
end

%get number of pixels in each direction
[numPixelsX,numPixelsY] = size(image);

%Divide image by bit depth, to normalize it between 0 and 1
image = double(image)/(2^bitDepth-1);

% %filter image to reduce noise
% imageF = filterGauss2D(image,min(1,psfSigma));

%get local maxima information from cands
status = vertcat(cands.status);
locMaxKeep = find(status==1);
locMaxPos = vertcat(cands.Lmax);
locMaxPos = locMaxPos(locMaxKeep,:);
locMaxAmp = vertcat(cands.amp);
locMaxAmp = locMaxAmp(locMaxKeep);
bgAmp = vertcat(cands.IBkg);
bgAmp = bgAmp(locMaxKeep);

%get number of local maxima
numLocMax = length(locMaxKeep);

%% Centroid calculation

%get half the PSF range in pixels
psfHalfRange = round(2*psfSigma);

%go over all local maxima
for i = 1 : numLocMax
    
    %get position
    midPixel = locMaxPos(i,:);
    
    %get part of image relevant for this local maximum
    minCoord1 = max(midPixel(1) - psfHalfRange,1);
    maxCoord1 = min(midPixel(1) + psfHalfRange,numPixelsX);
    minCoord2 = max(midPixel(2) - psfHalfRange,1);
    maxCoord2 = min(midPixel(2) + psfHalfRange,numPixelsY);
    imageLocMax = image(minCoord1:maxCoord1,minCoord2:maxCoord2);
    
    %calculate centroid in small image - notice that centroid comes back in
    %image coordinate system
    ce = centroid2D(imageLocMax);
    
    %shift to coordinates in overall image
    ce = ce(2:-1:1) + midPixel - psfHalfRange - 1;
    
    %store information in structure "detectedFeatures" - coordinates are
    %converted to image coordinate system, for consistency with
    %detectSubResFeatures2D
    detectedFeatures.xCoord(i,:) = [ce(2) 0.2];
    detectedFeatures.yCoord(i,:) = [ce(1) 0.2];
    detectedFeatures.amp(i,:) = [locMaxAmp(i)-bgAmp(i) 0];
    
end

%save output if requested
if saveResults
    save('detectedFeatures','detectedFeatures');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%visualization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if visual
    
    %make 3 layers out of original image (normalized)
    imageNorm = image/max(image(:));
    imageN3 = repmat(imageNorm,[1 1 3]);
    
    %place zeros in pixels of maxima from cands
    posL = (locMaxPos(:,2)-1)*numPixelsX + locMaxPos(:,1);
    for j=1:3
        imageN3(posL + (j-1)*numPixelsX*numPixelsY) = 0;
    end
    
    %place zeros in pixels of maxima from centroid calculation
    posC = (round(detectedFeatures.xCoord(:,1))-1)*numPixelsX ...
        + round(detectedFeatures.yCoord(:,1));
    for j=1:3
        imageN3(posC + (j-1)*numPixelsX*numPixelsY)=0;
    end
    
    %label maxima from cands in blue
    imageN3(posL + 2*numPixelsX*numPixelsY) = 1;
    
    %label maxima from mixture-model fitting in red
    %a maximum from mixture-model fitting that falls in the same pixel
    %as that from cands will appear in magenta
    imageN3(posC)=1;
    
    %plot image
    imtool(imageN3);
    
end

%%%%% ~~ the end ~~ %%%%%

