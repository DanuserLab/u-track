function [projData,M]=postProcessMTTracks(projData,tracksFinal,movieInfo,timeRange,varargin)

% Check additional input
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
ip =inputParser;
ip.addOptional('remBegEnd',true,@isscalar);
ip.addParamValue('fgapReclassScheme',1,@isscalar);
ip.addParamValue('bgapReclassScheme',1,@isscalar);
ip.parse(varargin{:})
remBegEnd = ip.Results.remBegEnd; % Flag to turn on or off the remove beginning and end 
fgapReclassScheme = ip.Results.fgapReclassScheme;
bgapReclassScheme = ip.Results.bgapReclassScheme;

% get interpolated positions for gaps and calculate velocities
[trackedFeatureInfo,trackedFeatureInfoInterp,trackInfo,trackVelocities,timeRange]=...
    plusTipGetVelocitiesFromMat(tracksFinal,movieInfo,3,timeRange);

%get number of tracks and number of time points
[nTracks,numTimePoints] = size(trackedFeatureInfo);
numTimePoints=numTimePoints/8;


% without interpolation yet
x = trackedFeatureInfo(:,1:8:end);
y = trackedFeatureInfo(:,2:8:end);

% initialize matrices for particle indices, area, and intensity
movieInfoIdx=nan(nTracks,numTimePoints);
if isfield(movieInfo,'ecc')
    featArea=nan(nTracks,numTimePoints);
    featInt =nan(nTracks,numTimePoints);
end

for iFrame=1:numTimePoints
    % these are the track numbers which exist in iFrame
    existCoordIdx=find(~isnan(x(:,iFrame)));
    % these are the corresponding xy-coordinates
    xi=x(:,iFrame); xi(isnan(xi))=[];
    yi=y(:,iFrame); yi(isnan(yi))=[];

    if ~isempty(xi)
        % distance matrix reveals where particles coincide with those recorded
        % in movieInfo
        D=createDistanceMatrix([xi,yi],[movieInfo(iFrame,1).xCoord(:,1),movieInfo(iFrame,1).yCoord(:,1)]);
        [r,c]=find(D==0); % r=track index, c=frame

        [newR,idx]=sort(r); % re-order based on track
        featIdx=c(idx); % movieInfo particle index, sorted to correspond to track indices

        % fill in movieInfoIdx with indices from particles stored in movieInfo
        movieInfoIdx(existCoordIdx,iFrame)=featIdx;
        if isfield(movieInfo,'ecc')
            % fill in particle area (pixels) at corresponding particles
            featArea(existCoordIdx,iFrame)=movieInfo(iFrame,1).amp(featIdx,1);
            % fill in particle max intensity at corresponding particles
            featInt (existCoordIdx,iFrame)=movieInfo(iFrame,1).int(featIdx,1);
        end
    end
end

% projData.trackingParameters.maxGapLength=gapCloseParam.timeWindow;
% projData.trackingParameters.minTrackLen=gapCloseParam.minTrackLen;
% projData.trackingParameters.minSearchRadius=costMatrices(1,1).parameters.minSearchRadius;
% projData.trackingParameters.maxSearchRadius=costMatrices(1,1).parameters.maxSearchRadius;
% projData.trackingParameters.maxForwardAngle=costMatrices(1,2).parameters.maxFAngle;
% projData.trackingParameters.maxBackwardAngle=costMatrices(1,2).parameters.maxBAngle;
% projData.trackingParameters.backVelMultFactor=costMatrices(1,2).parameters.backVelMultFactor;
% projData.trackingParameters.fluctRadius=costMatrices(1,2).parameters.fluctRad;

projData.nTracks = nTracks;
projData.nFrames = numTimePoints;

% figure out which frames were used in detection
m=struct2cell(movieInfo); m=m(1,:); 
detExists=find(cellfun(@(x) ~isempty(x),m)); 
sF=min(detExists); eF=max(detExists);

% frame ranges for each step
% projData.detectionFrameRange=[sF eF];
% projData.trackingFrameRange=[costMatrices(1).parameters.startFrame costMatrices(1).parameters.endFrame];
% projData.postTrackFrameRange = timeRange;

% coordinate/area/intensity info from detected particles
projData.xCoord = trackedFeatureInfoInterp(:,1:8:end);
projData.yCoord = trackedFeatureInfoInterp(:,2:8:end);
if isfield(movieInfo,'ecc')
    projData.featArea = featArea;
    projData.featInt = featInt;
end

% get frame-to-frame displacement for growth only (not forward/backward gaps)
frame2frameDispPix=sqrt(diff(x,1,2).^2+diff(y,1,2).^2);
% get rid of NaNs and linearize the vector
projData.frame2frameDispPix=frame2frameDispPix(~isnan(frame2frameDispPix(:)));

% get change in velocity between frame *pairs* for segments only
pair2pairDiffPix=diff(frame2frameDispPix,1,2);
% get rid of NaNs and linearize the vector
projData.pair2pairDiffPix=pair2pairDiffPix(~isnan(pair2pairDiffPix(:)));
% std (microns/min) of delta growthSpeed btw frames
projData.pair2pairDiffMicPerMinStd=std(pixPerFrame2umPerMin(projData.pair2pairDiffPix,projData.secPerFrame,projData.pixSizeNm));


% get all particle nearest neighbor distances from all frames in one vector
NNdist=nan(length(vertcat(movieInfo.xCoord)),1);
count=1;
for iFrame=5:length(movieInfo)
    xCoord = movieInfo(iFrame).xCoord;
    yCoord = movieInfo(iFrame).yCoord;
    
    if ~isempty(xCoord)
        xCoord=xCoord(:,1);
        yCoord=yCoord(:,1);

        D=createDistanceMatrix([xCoord yCoord],[xCoord yCoord]);
        [sD,idx]=sort(D,2);

        NNdist(count:count+length(xCoord)-1)=sD(:,2);

    end
    count=count+length(xCoord);
end

% median NN dist
projData.medNNdistWithinFramePix=nanmedian(NNdist);

% get mean displacement to median NN distance ratio
projData.meanDisp2medianNNDistRatio = mean(projData.frame2frameDispPix)/projData.medNNdistWithinFramePix;

% convert interpolated velocities to microns per minute
[projData.frame2frameVel_micPerMin]=pixPerFrame2umPerMin(trackVelocities.frame2frame,projData.secPerFrame,projData.pixSizeNm);
projData.segGapAvgVel_micPerMin=[];

% concatenate all segments and gaps into n x 4 matrices, then add info:
% [trackNum startFrame endFrame velocity seg/gapType trackLengthFrames]
segs   = vertcat(trackInfo.seg);
fgaps  = vertcat(trackInfo.fgap);
bgaps  = vertcat(trackInfo.bgap);
ugaps  = vertcat(trackInfo.ugap);

compositeMatrix = [];
if ~isempty(segs)
    segs =  [segs   1*ones(size(segs,1),1)];
    compositeMatrix = [compositeMatrix; segs];
end
if ~isempty(fgaps)
    fgaps = [fgaps  2*ones(size(fgaps,1),1)];
    compositeMatrix = [compositeMatrix; fgaps];
end
if ~isempty(bgaps)
    bgaps = [bgaps  3*ones(size(bgaps,1),1)];
    compositeMatrix = [compositeMatrix; bgaps];
end
if ~isempty(ugaps)
    ugaps = [ugaps  4*ones(size(ugaps,1),1)];
    compositeMatrix = [compositeMatrix; ugaps];
end



% put segs/gaps into one matrix and sort to see track profiles in order
aT=sortrows(compositeMatrix,[1 2]);

% lifetime is subtrack length in frames
lifeTimes=aT(:,3)-aT(:,2);
% get total distance traveled over all subtracks
totalDispPix=aT(:,4).*lifeTimes;
% add lifetime and total distplacement to matrix
aT=[aT lifeTimes totalDispPix];

% convert pix/frame to micron/min velocities
[aT(:,4)] = pixPerFrame2umPerMin(aT(:,4),projData.secPerFrame,projData.pixSizeNm);

% aT will now contain consolidated rows (we will further use this one to 
% calculate the stats), while aTreclass will be stored in projData.

[aT,aTreclass,dummy,projData]=plusTipMergeSubtracks(projData,aT, ...
    fgapReclassScheme, bgapReclassScheme);

% assign the matrix retaining where growth fgaps are indicated with
% trackType=5 (This structure will be read into plotting functions so one
% can visualize which tracks have been reclassified)
projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix=aTreclass;


% recalculate segment average speeds to reflect consolidation
projData.segGapAvgVel_micPerMin=nan(size(projData.frame2frameVel_micPerMin));
for iSub=1:size(aT,1)
    projData.segGapAvgVel_micPerMin(aT(iSub,1),aT(iSub,2):aT(iSub,3)-1)=aT(iSub,4);
end

% get track numbers that contain an fgap or bgap
projData.tracksWithFgap = unique(aT(aT(:,5)==2,1));
projData.tracksWithBgap = unique(aT(aT(:,5)==3,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add 2 more columns corresponding to extra subTrack information
% Though obvious here when partition these subtracks based on 
% regional criteria it can become ambigious. 
% Therefore it is useful to mark these here while the dataStruct for the 
% whole cell compound tracks are in tact.
% Column 8: Nuc Event 1= yes 0 = no
% Column 9: Growth Before Term Event = 1, Growth Before Fgap = 2, 
% Growth before Bgap = 3, Growth before undefined gap = 4

[dummy nucEventsIdx dummy] =  unique(aT(:,1),'first');
[dummy termEventsIdx dummy] = unique(aT(:,1),'last'); 

aT(nucEventsIdx,8) = 1; 
allIdx = 1:length(aT(:,1)); 
nonNucIdx = setdiff(allIdx,nucEventsIdx); 

aT(nonNucIdx,8) = 0; 

aT(termEventsIdx,9) = 1;

fIdx = find(aT(:,5) == 2); 
bIdx = find(aT(:,5) == 3); 
uIdx = find(aT(:,5) == 4); 

beforeFgapIdx=fIdx-1;
beforeBgapIdx = bIdx-1; 
beforeUgapIdx = uIdx-1; 

aT(beforeFgapIdx,9) = 2; 
aT(beforeBgapIdx,9) = 3; 
aT(beforeUgapIdx,9) = 4; 
aT([fIdx;bIdx;uIdx],9)= 0; 


% perform lifetime and displacement unit conversions
aT(:,6)=aT(:,6).* projData.secPerFrame; % convert lifetimes to seconds
aT(:,7)=aT(:,7).*(projData.pixSizeNm/1000); % convert displacements to microns

% mark if part of compound versus non-compound track 
   gapIdx = sort([fIdx;bIdx]); 
   compIdx= unique(sort([gapIdx ; (gapIdx+1) ; (gapIdx -1)]));
   % set those part of a compound track to in column 10 to 1 
   % so marked for later partitioning
   aT(compIdx,10) = 1;
   
   compDataMat = aT(compIdx,:); 
 
   
    
   
 % Segregate Tracks That Are Exclusively From Single Tracks
   singleDataMat = aT;
   singleDataMat(compIdx,:) = [];
   
   
   % remove uIdx 
   uIdxSingleMat = find(singleDataMat(:,5) == 4);
   
   toRemove = sort([uIdxSingleMat;uIdxSingleMat+1;uIdxSingleMat-1]); 
   
   singleDataMat(toRemove(toRemove~=0),:) = []; 
   

   if remBegEnd == 1 
      compDataMat = plusTipRemBegEnd(compDataMat,projData,1); 
      singleDataMat = plusTipRemBegEnd(singleDataMat,projData,1);  
   end 
       
     projData.compDataMat = compDataMat;  
     projData.singleDataMat = singleDataMat; 
   
  
% save this dataStruct for subRoi 
% partitioning and stat calculations from pooled data.  
projData.mergedDataMatAllSubTracksConverted = aT; 


   
   


if remBegEnd == 1
% Remove growth subtracks that start in the first frame and end in the
% last frame (as well as flanking fgap and bgaps) 

[dataMatCrpSecMic projData] = plusTipRemBegEnd(aT,projData);
projData.remBegEnd = 'yes';
else 
    dataMatCrpSecMic = aT;
    projData.remBegEnd = 'no'; 
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate stats using the matrix beginning/end data has been
% . M records speeds (microns/min), lifetimes (sec), and
% displacements (microns) for growths, fgaps,and bgaps.
% note the input 0,0 just tells the program that it is NOT 
% calling from subRoi or poolGroupData
[projData,M]=plusTipDynamParam(dataMatCrpSecMic,projData,0,0);
