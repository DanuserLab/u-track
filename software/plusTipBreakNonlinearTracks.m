function [newTrackedFeatureInfo,newTrackedFeatureIndx,newNnDistFeatures,diagnosticTrackLinearity]=plusTipBreakNonlinearTracks(trackedFeatureInfo,trackedFeatureIndx,nnDistFeatures)
% plusTipBreakNonlinearTracks splits up tracks not following uni-directional behavior
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

% this function splits up individial tracks where the displacement vectors
% created by consecutive frame-frame pairs in a track differ in direction
% by more than 45 degrees. for each instance like this, one new track is
% born. 
%
% some instances where the angle is greater than 45 degrees are retained if
% one or both of the vectors is very short (<3rd percentile), which may
% simply reflect the uncertainty in the detected position if the features
% are not moving much per frame.
%
% output has the same format as input. NOTE: nnDistFeatures will now be
% nonsensical, but these are currently not used in the gap closing cost
% matrix.  many links will now shorter than minLength, but those can be
% filtered out by post-processing step.

% this function is currently called in conditional statement during
% trackCloseGapsKalman


%get total number of tracks
[nTracks, nFrames] = size(trackedFeatureIndx);

% find track start and end frames
trackSEL = getTrackSEL(trackedFeatureInfo);

% extract x and y coordinates for all tracks for all frames (NaN will be
% used as a place holder if no linked feature exists in that frame)
% trackedFeatureInfo is such that row = track #  
% column = [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...xn,yn,zn,an,dxn,dyn,dzn]
% where n = total number of frames
% and x,y are the respective coordinates, a is the "intensity" of the
% feature (in this case area)
% and dx,dy,dz,da are the uncertainty associated with each param.  
% If no linked feature exists in that frame: NaN will be used
% as a place holder
px=trackedFeatureInfo(:,1:8:end);
py=trackedFeatureInfo(:,2:8:end);
%amp = trackedFeatureInfo(:,4:8:end); 

% get vector coordinates between linked features
vx=diff(px,1,2);
vy=diff(py,1,2);
vmag=sqrt(vx.^2+vy.^2);

% first vector matrix
v1x=vx(:,1:end-1);
v1y=vy(:,1:end-1);
v1mag=sqrt(v1x.^2+v1y.^2);

% second vector matrix
v2x=vx(:,2:end);
v2y=vy(:,2:end);
v2mag=sqrt(v2x.^2+v2y.^2);

% cos of angle between consecutive vectors (displacements)
cosV12=(v1x.*v2x+v1y.*v2y)./(v1mag.*v2mag);


%% Additional Diagnostic to assess linearity of tracking data: (MB July 2010)
% Gives the percentage of total displacement vectors
% created by consecutive frame-frame pairs which differ in direction
% by the following angles:
% 0(perfect coorelation)-45, 45-90, 90-135, 135-180(perfect anti-coorelation)

%Bin according to the angle between the displacement vectors 
% Retains Original Matrix Dimensions: Puts "1" where angle of displacement 
% vectors lies between respective values and "0" everywhere else. 
angles0to45 = cosV12 > cos(45*pi/180);
angles45to90 = (cosV12 > cos(90*pi/180)) & (cosV12 < cos(45*pi/180));
angles90to135 = (cosV12 > cos(135*pi/180)) & (cosV12 < cos(90*pi/180));
angles135to180 = cosV12 < cos(135*pi/180);

%Calculate Number of Total Displacement Vector Angles Measured For
%Percent
% Calculation: Subtract out the place holders in the cosV12 matrix which 
% are Nan (not a number) to get the total number of angles included in
% the binning 
[dim1 dim2] = size(cosV12);
totalAngles = dim1*dim2 - sum(sum(isnan(cosV12)));

%Calculate Percent Populations From Respective Bins 
diagnosticTrackLinearity.percentAngles0to45 = (sum(sum(angles0to45))/totalAngles)*100;
diagnosticTrackLinearity.percentAngles45to90 = (sum(sum(angles45to90))/totalAngles)*100;
diagnosticTrackLinearity.percentAngles90to135 = (sum(sum(angles90to135))/totalAngles)*100;
diagnosticTrackLinearity.percentAngles135to180 = (sum(sum(angles135to180))/totalAngles)*100;


 
%%
% Originally assumed max angle was 45 degrees however by setting this value
% to 180 you bypass this. If max angle is set to 180 it will not break any
% of the non-linear tracks.  Main reason for calling this function with
% this setting will be to gain the above diagnostic track linearity
% information which may be a helpful statistic. Will later incorporate this 
% information that does not call this function to save computer time. 
noBreak = 0;

if noBreak  == 1
    max_angle = 180;
    max_angle2 = 0;
else 
max_angle = 45; % set to < 180 to break tracks
max_angle2 = 135;
end 

% Note: It is potentially useful to modifiy the max angle to retain some 
% non-linear tracks depending on data set at hand: See above diagnostic 
% to help trouble-shoot this parameter
cosMax1=cos(max_angle*pi/180);
cosMax2 = cos(max_angle2*pi/180);

% lower bound displacement - if smaller than this, may just be jitter
lb=prctile(vmag(:),10); 

% keep track of where cos or displacement vectors are NaNs 
% (nanMat = matrix with value '1' where cos(angle) value is numerical and NaN at all
% other coordinates)
nanMat=swapMaskValues(isnan(cosV12) | isnan(v1mag) | isnan(v2mag),[0 1],[1 NaN]);

% these are within forward cone, so they're ok
if noBreak == 1
    okAngles = cosV12 >=cosMax1;
else 
    okAngles1 = cosV12 >= cosMax1;
    okAngles2 = cosV12 <= cosMax2; 

okAngles = okAngles1 + okAngles2;
end



% these are not in the forward cone but one or both of the vectors is shorter 
% than the nth percentile of all vectors, so maybe just jitter
okLength =   (cosV12 > cosMax2 | cosV12 < cosMax1) & (v1mag<lb | v2mag<lb); 


diagnosticTrackLinearity.percentLinksBroken = (totalAngles-sum(sum(okAngles)))/totalAngles*100;

diagnosticTrackLinearity.percentLinksSavedFromBreak = (sum(sum(okLength))/totalAngles)*100;


% if either the angle or the length criterion isn't met, then it's a bad
% link which we will break
badLinks=swapMaskValues(nanMat.*(okAngles | okLength),[0 1],[1 0]);


[badTrackIdx badTrackVecPair]=find(badLinks==1);
badTrackVecHead=badTrackVecPair+1;

doPlot=0;
if doPlot==1
    % plot the first 50 bad tracks
    b=badTrackIdx(1:50);
    figure
    plot(px(b,:)',py(b,:)')
    hold on
    x=px(b,:)'; x=x(:);
    y=py(b,:)'; y=y(:);
    scatter(x,y,'b')
    % plot break points in red
    for i=1:length(b)
        a=badTrackIdx(badTrackIdx==b(i));
        d=badTrackVecPair(badTrackIdx==b(i));
        x=px(sub2ind(size(badLinks),a,d)+length(px));
        y=py(sub2ind(size(badLinks),a,d)+length(px));

        scatter(x,y,'r')
    end
end


badLinkIdx=[badTrackIdx badTrackVecHead];
badLinkIdx=sortrows(badLinkIdx,1); % sorted indices [trackNumber headPosition]

[nBadLinks, trackIdxWithBadLink] = getMultiplicity(badLinkIdx(:,1));
% n links to break creates n+1 segments. but, since we retain the original row
% for the first segment, we only need to add n rows
nRows2add = sum(nBadLinks);

% initialize new matrices to contain new rows
newTrackedFeatureIndx = [trackedFeatureIndx; zeros(nRows2add,nFrames)];
newNnDistFeatures = [nnDistFeatures; nan(nRows2add,nFrames)];
newTrackedFeatureInfo = [trackedFeatureInfo; nan(nRows2add,8*nFrames)];

counter=nTracks+1;
for i=1:length(trackIdxWithBadLink);
    idx = badLinkIdx(badLinkIdx(:,1)==trackIdxWithBadLink(i),2);
    newSegS = [trackSEL(trackIdxWithBadLink(i),1); idx+1];  % new segment start
    newSegE = [idx; trackSEL(trackIdxWithBadLink(i),2)];    % new segment end

    for j=2:length(newSegS) % leave the first one as it is - will contain first segment
        % assign original values to new rows
        newTrackedFeatureIndx(counter,newSegS(j):newSegE(j)) = newTrackedFeatureIndx(trackIdxWithBadLink(i),newSegS(j):newSegE(j));
        newNnDistFeatures(counter,newSegS(j):newSegE(j)) = newNnDistFeatures(trackIdxWithBadLink(i),newSegS(j):newSegE(j));
        newTrackedFeatureInfo(counter,8*(newSegS(j)-1)+1:8*newSegE(j)) = newTrackedFeatureInfo(trackIdxWithBadLink(i),8*(newSegS(j)-1)+1:8*newSegE(j));

        % erase original values from original rows
        newTrackedFeatureIndx(trackIdxWithBadLink(i),newSegS(j):newSegE(j)) = zeros(size(newSegS(j):newSegE(j)));
        newNnDistFeatures(trackIdxWithBadLink(i),newSegS(j):newSegE(j)) = nan(size(newSegS(j):newSegE(j)));
        newTrackedFeatureInfo(trackIdxWithBadLink(i),8*(newSegS(j)-1)+1:8*newSegE(j)) = nan(size(8*(newSegS(j)-1)+1:8*newSegE(j)));

        counter=counter+1;
    end
end

doPlot=0;
if doPlot==1
    
    % plot broken links in red, retained tracks in blue
    px=trackedFeatureInfo(:,1:8:end)';
    py=trackedFeatureInfo(:,2:8:end)';
    
    % limit t
    px=px(5:25,:); py=py(5:25,:);   
    
    x=px(:); x(isnan(x))=[];
	y=py(:); y(isnan(y))=[];
    
    figure
    plot(px(:,:),py(:,:),'r') % original
    hold on
    scatter(x,y,'.b')
    
    px=newTrackedFeatureInfo(:,1:8:end)';
    py=newTrackedFeatureInfo(:,2:8:end)';
    px=px(5:25,:); py=py(5:25,:);
    
    plot(px(:,:),py(:,:),'b') % new

    axis equal
end

%rearrange "newTrackedFeatureIndx" such that tracks are sorted in ascending order by their
%starting point. Note that this ends up also arranging tracks starting at the
%same time in descending order from longest to shortest.

trackSEL = getTrackSEL(newTrackedFeatureInfo);
[list,indx] = sortrows(trackSEL,[1 -3]);

% rearrange data
newTrackedFeatureIndx = newTrackedFeatureIndx(indx,:);
newNnDistFeatures = newNnDistFeatures(indx,:);
newTrackedFeatureInfo = newTrackedFeatureInfo(indx,:);



