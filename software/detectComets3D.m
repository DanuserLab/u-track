function [movieInfo,featMapFinal] = detectComets3D(I,stepSize,thresh,pxlscale)
% Detect plusTip comets in a prefiltered image using a watershed method
%
% Synopsis: movieInfo = detectComets(I,stepSize,thresh)
%
% Input:
%     I - a prefiltered image
%
%     stepSize - the intensity of steps
%
%     thresh - the minimum value for detecting signal
%
%
% Ouput:
%     movieInfo   : a structure compatible with the tracker containing the
%     x/y coordinates, comet amplitude, maximum intensity and eccentricity
%
% Copyright (C) 2012 LCCB 
%
% This file is part of plusTipTracker.
% 
% plusTipTracker is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% plusTipTracker is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with plusTipTracker.  If not, see <http://www.gnu.org/licenses/>.
% 
% 
%
% Copyright (C) 2020, Danuser Lab - UTSouthwestern 
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


% Wes Legant 2014
% Ported from an initial 2D version by Kathryn Applegate, 2010 and Sebastien Besson, March 2012


% we assume each step size down the intensity profile should be on
% the order of the size of the background std; here we find how many
% steps we need and what their spacing should be. we also assume peaks
% should be taller than 3*std
nSteps = round((nanmax(I(:))-thresh)/double(stepSize));
threshList = linspace(nanmax(I(:)),thresh,nSteps);

movieInfo=struct('xCoord',[],'yCoord',[],'zCoord',[],'amp',[],'int',[]);
if nSteps<1,return; end

% Handle case where nSteps==1
if nSteps==1, slice2 = I>threshList(1); end


% compare features in z-slices startest from the highest one
for j = 1:length(threshList)-1
    % slice1 is top slice; slice2 is next slice down
    % here we generate BW masks of slices
    if j==1
        slice1 = I>threshList(j);
    else
        slice1 = slice2;
    end
    slice2 = I>threshList(j+1);
    
    % now we label them using the "bwlabeln" function from matlab which
    % labels connected components in a N-D binary image
    featMap1 = bwlabeln(slice1);
    featMap2 = bwlabeln(slice2);
    
    % get the regionproperty 'PixelIdxList' using "regionprops" function in matlab
    featProp2 = regionprops(featMap2,'PixelIdxList');
    
    % loop thru slice2 features and replace them if there are 2 or
    % more features from slice1 that contribute
    for iFeat = 1:max(featMap2(:))
        pixIdx = featProp2(iFeat,1).PixelIdxList; % pixel indices from slice2
        featIdx = unique(featMap1(pixIdx)); % feature indices from slice1 using same pixels
        featIdx(featIdx==0) = []; % 0's shouldn't count since not feature
        if length(featIdx)>1 % if two or more features contribute...
            slice2(pixIdx) = slice1(pixIdx); % replace slice2 pixels with slice1 values
        end
    end
    
end


% label slice2 again and get region properties
featProp2 = regionprops(logical(slice2),'PixelIdxList','Area');

% here we sort through features and retain only the "good" ones
% we assume the good features have area > 2 pixels
goodFeatIdx = vertcat(featProp2(:,1).Area)>2;
%    goodFeatIdxI = find(vertcat(featProp2(:,1).MaxIntensity)>2*cutOffValueInitInt);
%    goodFeatIdx = intersect(goodFeatIdxA,goodFeatIdxI);

% make new label matrix and get props
featureMap = zeros(size(I));
featureMap(vertcat(featProp2(goodFeatIdx,1).PixelIdxList)) = 1;
[featMapFinal,nFeats] = bwlabeln(featureMap);

featPropFinal = regionprops(featMapFinal,I,'PixelIdxList',...
    'Area','WeightedCentroid','MeanIntensity','MaxIntensity','PixelValues'); %'Extrema'

if nFeats==0, return; end

% centroid coordinates with 0.5 uncertainties for Khuloud's tracker
yCoord = 0.5*ones(nFeats,2);
xCoord = 0.5*ones(nFeats,2);
zCoord = 0.5*ones(nFeats,2);
temp = vertcat(featPropFinal.WeightedCentroid);
zCoord(:,1) = temp(:,3);
yCoord(:,1) = temp(:,2);
xCoord(:,1) = temp(:,1);

% SB: shoudl we use meanIntensity instead>>>
% area
featArea = vertcat(featPropFinal(:,1).Area);
amp = zeros(nFeats,2);
amp(:,1) = featArea;

% intensity
featInt = vertcat(featPropFinal(:,1).MaxIntensity);
featI = zeros(nFeats,2);
featI(:,1) = featInt;

verDate=version('-date');

% if str2double(verDate(end-3:end))>=2008 % can only calculate eccentricity
%     % if using version of matlab older than 2008
%     
%     %eccentricity
%     featEcc = vertcat(featPropFinal(:,1).Eccentricity);
%     featE = zeros(nFeats,2);
%     featE(:,1) = featEcc;
%     
% end

% make structure compatible with Khuloud's tracker
movieInfo.xCoord = xCoord*pxlscale(1);
movieInfo.yCoord = yCoord*pxlscale(2);
movieInfo.zCoord = zCoord*pxlscale(3);
movieInfo.amp = amp;
movieInfo.int = featI;
end
