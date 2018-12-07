function [costMat,nonlinkMarker,indxMerge,numMerge,indxSplit,numSplit,...
    errFlag] = costMatRandomDirectedSwitchingMotionCloseGaps(trackedFeatInfo,...
    trackedFeatIndx,trackStartTime,trackEndTime,costMatParam,gapCloseParam,...
    kalmanFilterInfo,nnDistLinkedFeat,probDim,movieInfo)
%costMatRandomDirectedSwitchingMotionCloseGaps provides a cost matrix for closing gaps and capturing merges/splits using Kalman filter information
%
%SYNOPSIS [costMat,nonlinkMarker,indxMerge,numMerge,indxSplit,numSplit,...
%    errFlag] = costMatRandomDirectedSwitchingMotionCloseGaps(trackedFeatInfo,...
%    trackedFeatIndx,trackStartTime,trackEndTime,costMatParam,gapCloseParam,...
%    kalmanFilterInfo,nnDistLinkedFeat,probDim,movieInfo)
%
%INPUT  trackedFeatInfo: The positions and amplitudes of the tracked
%                        features from linkFeaturesKalman.
%                        Number of rows = number of tracks.
%                        Number of columns = 8*number of frames.
%                        Each row consists of
%                        [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                        in image coordinate system (coordinates in
%                        pixels). NaN is used to indicate time points
%                        where the track does not exist.
%       trackedFeatIndx: Connectivity matrix of features between frames.
%                        Rows indicate continuous tracks, while columns
%                        indicate frames. A track that ends before the
%                        last time point is followed by zeros, and a track
%                        that starts at a time after the first time point
%                        is preceded by zeros.
%       trackStartTime : Starting time of all tracks.
%       trackEndTime   : Ending time of all tracks.
%       costMatParam   : Structure containing variables needed for cost
%                        calculation. Contains the fields:
%             .linearMotion   : 0 for only random motion;
%                               1 for random + directed motion;
%                               2 for random + directed motion with the
%                               possibility of instantaneous switching to
%                               opposite direction (but same speed),
%                               i.e. something like 1D diffusion.
%             .minSearchRadius: Minimum allowed search radius (in pixels).
%             .maxSearchRadius: Maximum allowed search radius (in pixels).
%                               This value is the maximum search radius
%                               between two consecutive frames as when
%                               linking between consecutive frames. It will
%                               be calcualted for different time gaps
%                               based on the scaling factor of Brownian
%                               motion (expanding it will make use of the
%                               fields .brownScaling and .timeReachConfB).
%             .brownStdMult   : Factor multiplying Brownian
%                               displacement std to get search radius.
%                               Vector with number of entries equal to
%                               gapCloseParam.timeWindow (defined below).
%             .linStdMult     : Factor multiplying linear motion std to get
%                               search radius. Vector with number of entries
%                               equal to gapCloseParam.timeWindow (defined
%                               below).
%             .brownScaling   : Power with which the Brownian part of the
%                               search radius scales with time. It has 2
%                               elements, the first indicating the power
%                               before timeReachConfB (see below) and the
%                               second indicating the power after
%                               timeReachConfB.
%             .linScaling     : Power with which the Linear part of the
%                               search radius scales with time. It has 2
%                               elements, the first indicating the power
%                               before timeReachConfL (see below) and the
%                               second indicating the power after
%                               timeReachConfL
%             .timeReachConfB : Time gap for reaching confinement for
%                               2D Brownian motion. For smaller time gaps,
%                               expected displacement increases with
%                               (time gap)^brownScaling. For larger time gaps,
%                               expected displacement increases slowly, with
%                               (time gap)^0.01.
%             .timeReachConfL : Time gap for reaching confinement for
%                               linear motion. Time scaling similar to
%                               timeReachConfB above.
%             .lenForClassify : Minimum length of a track to classify it as
%                               directed or Brownian.
%             .maxAngleVV     : Maximum allowed angle between two
%                               directions of motion for potential linking
%                               (in degrees).
%             .useLocalDensity: 1 if local density of features is used to expand
%                               their search radius if possible, 0 otherwise.
%             .nnWindow       : Time window to be used in estimating the
%                               nearest neighbor distance of a track at its start
%                               and end.
%           Optional fields ...
%             .ampRatioLimit  : Minimum and maximum allowed ratio between
%                               the amplitude of a merged feature and the
%                               sum of the amplitude of the two features
%                               making it.
%                               Default: [], in which case the amplitude is 
%                               not used for cost calculation.
%             .lftCdf         : Lifetime cumulative density function.
%                               Column vector, specifying cdf for
%                               lifetime = 0 to movie length.
%                               Default: [], in which case cdf is not used.
%             .gapPenalty     : Penalty for increasing temporary
%                               disappearance time, to be used in gap
%                               closing cost. Disappearing for n frames,
%                               i.e. closing a gap of n+1 frames,
%                               gets a penalty of gapPenalty^n.
%                               Default: 1, i.e. no penalty.
%             .resLimit       : Resolution limit. Used to expand search 
%                               radius for merging and splitting, if
%                               motion- and density-based search radius is
%                               smaller than resolution limit.
%                               Default: 0, i.e. resolution limit is not
%                               used for search radius expansion.
%             .gapExcludeMS   : 1 to impose that the possibility of gap closing
%                               prohibits merging/splitting for an
%                               end/start, 0 to allow the two to compete.
%                               In other words, if 1, then gaps are given
%                               absolute priority.
%                               Default: 0, i.e. no prioritization.
%             .strategyBD     : Strategy for calculating the birth and
%                               death cost. 
%                               If a positive number, then this is taken as
%                               the percentile of the gap closing and merging
%                               and splitting cost distribution to use for
%                               calculating the birth and death cost.
%                               If 0, then the percentile is calculated
%                               from the structure (in essence the
%                               over-connectivity) of the cost matrix.
%                               If -1, then the birth and death cost is not
%                               calculated as a percentile but rather its
%                               value is determined by the definitions of
%                               the various factors contributing to the gap
%                               closing and merging and splitting costs.
%                               Default: 0 (middle option above).
%       gapCloseParam  : Structure containing variables needed for gap closing.
%                        Contains the fields:
%             .timeWindow : Largest time gap between the end of a track and the
%                           beginning of another that could be connected to it.
%             .tolerance  : Relative change in number of tracks in two
%                           consecutive gap closing steps below which
%                           iteration stops.
%             .mergeSplit : Logical variable with value 1 if the merging
%                           and splitting of trajectories are to be consided;
%                           0 if merging and splitting are not allowed, 2
%                           if only merging is allowed, and 3 if only
%                           splitting is allowed.
%       kalmanFilterInfo:Structure array with number of entries equal to
%                        number of frames in movie. Contains the fields:
%             .stateVec   : Kalman filter state vector for each
%                           feature in frame.
%             .stateCov   : Kalman filter state covariance matrix
%                           for each feature in frame.
%             .noiseVar   : Variance of state noise for each
%                           feature in frame.
%             .stateNoise : Estimated state noise for each feature in
%                           frame.
%             .scheme     : 1st column: propagation scheme connecting
%                           feature to previous feature. 2nd column:
%                           propagation scheme connecting feature to
%                           next feature.
%       nnDistLinkedFeat:Matrix indicating the nearest neighbor
%                        distances of features linked together within
%                        tracks.
%       probDim        : Problem dimensionality. 2 (for 2D) or 3 (for 3D).
%       movieInfo      : movieInfo as input to trackCloseGapsKalman. Not
%                        really used in this code, but needed for
%                        compatibility with other cost functions.
%
%OUTPUT costMat       : Cost matrix.
%       nonlinkMarker : Value indicating that a link is not allowed.
%       indxMerge     : Index of tracks that have possibly merged with
%                       tracks that end before the last time points.
%       numMerge      : Number of such tracks.
%       indxSplit     : Index of tracks from which tracks that begin after
%                       the first time point might have split.
%       numSplit      : Number of such tracks.
%       errFlag       : 0 if function executes normally, 1 otherwise.
%
%REMARKS
%
%Khuloud Jaqaman, April 2007
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

costMat = [];
nonlinkMarker = [];
indxMerge = [];
numMerge = [];
indxSplit = [];
numSplit = [];
errFlag = [];

%% Input

%check whether correct number of input arguments was used
if nargin ~= nargin('costMatRandomDirectedSwitchingMotionCloseGaps')
    disp('--costMatRandomDirectedSwitchingMotionCloseGaps: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%get cost matrix parameters
linearMotion = costMatParam.linearMotion;
minSearchRadius = costMatParam.minSearchRadius;
maxSearchRadius = costMatParam.maxSearchRadius;
brownStdMult = costMatParam.brownStdMult;
brownScaling = costMatParam.brownScaling;
timeReachConfB = costMatParam.timeReachConfB;
lenForClassify = costMatParam.lenForClassify;
useLocalDensity = costMatParam.useLocalDensity;
linStdMult   = costMatParam.linStdMult;
linScaling = costMatParam.linScaling;
timeReachConfL = costMatParam.timeReachConfL;
sin2AngleMax = (sin(costMatParam.maxAngleVV*pi/180))^2;
sin2AngleMaxVD = 0.5;
nnWindow = costMatParam.nnWindow;
if useLocalDensity
    closestDistScale = 2;
    maxStdMult = 100;
else
    closestDistScale = [];
    maxStdMult = [];
end
if isfield(costMatParam,'ampRatioLimit') && ~isempty(costMatParam.ampRatioLimit)
    minAmpRatio = costMatParam.ampRatioLimit(1);
    maxAmpRatio = costMatParam.ampRatioLimit(2);
    useAmp = 1;
else
    minAmpRatio = 0;
    maxAmpRatio = Inf;
    useAmp = 0;
end
if isfield(costMatParam,'lftCdf') && ~isempty(costMatParam.lftCdf)
    lftCdf = costMatParam.lftCdf;
    oneMinusLftCdf = 1 - lftCdf;
else
    lftCdf = [];
end
if isfield(costMatParam,'gapPenalty') && ~isempty(costMatParam.gapPenalty)
    gapPenalty = costMatParam.gapPenalty;
else
    gapPenalty = 1;
end
if isfield(costMatParam,'resLimit') && ~isempty(costMatParam.resLimit)
    resLimit = costMatParam.resLimit;
else
    resLimit = 0;
end
if isfield(costMatParam,'gapExcludeMS') && ~isempty(costMatParam.gapExcludeMS)
    gapExcludeMS = costMatParam.gapExcludeMS;
else
    gapExcludeMS = 0;
end
if isfield(costMatParam,'strategyBD') && ~isempty(costMatParam.strategyBD)
    strategyBD = costMatParam.strategyBD;
else
    strategyBD = 0;
end

%get gap closing/merging & splitting parameters
timeWindow = gapCloseParam.timeWindow;
mergeSplit = gapCloseParam.mergeSplit;

%make sure that timeReachConfB and timeReachConfL are <= timeWindow
timeReachConfB = min(timeReachConfB,timeWindow);
timeReachConfL = min(timeReachConfL,timeWindow);

%find the number of tracks to be linked and the number of frames in the movie
[numTracks,numFrames] = size(trackedFeatInfo);
numFrames = numFrames / 8;

%list the tracks that start and end in each frame
tracksPerFrame = repmat(struct('starts',[],'ends',[]),numFrames,1);
for iFrame = 1 : numFrames    
    tracksPerFrame(iFrame).starts = find(trackStartTime == iFrame); %starts
    tracksPerFrame(iFrame).ends = find(trackEndTime == iFrame); %ends
end

%% Pre-processing

%get the x,y-coordinates and amplitudes at the starts of tracks
coordStart = zeros(numTracks,probDim);
ampStart   = zeros(numTracks,1);
for iTrack = 1 : numTracks
    coordStart(iTrack,:) = full(trackedFeatInfo(iTrack,...
        (trackStartTime(iTrack)-1)*8+1:(trackStartTime(iTrack)-1)*8+probDim));
    ampStart(iTrack) = full(trackedFeatInfo(iTrack,(trackStartTime(iTrack)-1)*8+4));
end

%get the x,y-coordinates and amplitudes at the ends of tracks
coordEnd = zeros(numTracks,probDim);
ampEnd   = zeros(numTracks,1);
for iTrack = 1 : numTracks
    coordEnd(iTrack,:) = full(trackedFeatInfo(iTrack,...
        (trackEndTime(iTrack)-1)*8+1:(trackEndTime(iTrack)-1)*8+probDim));
    ampEnd(iTrack) = full(trackedFeatInfo(iTrack,(trackEndTime(iTrack)-1)*8+4));
end

%determine the types, velocities, noise stds, centers and mean
%displacements of all tracks
[trackType,xyzVelS,xyzVelE,noiseStd,trackCenter,trackMeanDispMag] = estimTrackTypeParamRDS(...
    trackedFeatIndx,trackedFeatInfo,kalmanFilterInfo,lenForClassify,probDim);

%if by chance some tracks are labeled linear when linearMotion=0, make them
%not linear
if linearMotion == 0
    trackType(trackType==1) = 0;
end

%find the 10th percentile of the noise standard deviation in order to use
%that for undetermined tracks whose mean displacement cannot be used to
%estimate their noise standard deviation
%(after removing std = 1 which indicates the simple initialization conditions)
%use 10% to be quite strict - basically, unless such a short track falls in
%the search area of a longer track, it won't get linked to anything
noiseStdAll = noiseStd(noiseStd ~= 1);
undetBrownStd = prctile(noiseStdAll,10);

%for undetermined tracks that have a mean displacement estimate (i.e. all
%tracks longer than 1 frame), use the mean displacement estimate to assign
%a noiseStd value (instead of the Kalman filter)
indx = find(noiseStd==1 & ~isnan(trackMeanDispMag));
noiseStd(indx) = trackMeanDispMag(indx)/sqrt(2);

%calculate the average mean displacement for all tracks, to assign to
%tracks that have no mean displacement estimate
meanDispAllTracks = nanmean(trackMeanDispMag);

%determine the search areas of all tracks
[longVecSAll,longVecEAll,...
    shortVecSAll,shortVecEAll,...
    shortVecS3DAll,shortVecE3DAll,...
    longVecSAllMS,longVecEAllMS,...
    shortVecSAllMS,shortVecEAllMS,...
    shortVecS3DAllMS,shortVecE3DAllMS,...
    longRedVecSAll,longRedVecEAll,...
    longRedVecSAllMS,longRedVecEAllMS] = ...
    ...
    getSearchRegionRDS(...
    ...
    xyzVelS,xyzVelE,noiseStd,trackType,undetBrownStd,timeWindow,...
    brownStdMult,linStdMult,timeReachConfB,timeReachConfL,...
    minSearchRadius,maxSearchRadius,useLocalDensity,closestDistScale,...
    maxStdMult,nnDistLinkedFeat,nnWindow,trackStartTime,trackEndTime,...
    probDim,resLimit,brownScaling,linScaling,linearMotion);

%% Gap closing

%find all pairs of ends and starts that can potentially be linked
%determine this by looking at time gaps between ends and starts
%and by looking at the distance between pairs
indxEnd2 = [];
indxStart2 = [];

%get the absolute upper limit of acceptable displacements in one frame
%as the maximum of (maximum velocity multiplied by probDim*linStdMult(1),
%maxSearchRadiusCG)
maxDispAllowed = max( max( abs([xyzVelS(:);xyzVelE(:)]) ) * probDim * linStdMult(1) * 3, ...
    maxSearchRadius );
% maxDispAllowed10 = 10 * maxDispAllowed;

%go over all frames until the one before last
for iFrame = 1 : numFrames - 1
    
    %find tracks that end in this frame
    endsToConsider = tracksPerFrame(iFrame).ends;
    
    for jFrame = iFrame + 1 : min(iFrame+timeWindow,numFrames)
        
        %find tracks that start in this frame
        startsToConsider = tracksPerFrame(jFrame).starts;
        
        %calculate the distance between ends and starts
        dispMat2 = createDistanceMatrix(coordEnd(endsToConsider,:),...
            coordStart(startsToConsider,:));
        
        %find possible pairs
        %         if (jFrame-iFrame) <= 2
        %             [indxEnd3,indxStart3] = find(dispMat2 <= (maxDispAllowed10 * sqrt(jFrame-iFrame)));
        %         else
        tmpFrame = jFrame-iFrame;
        [indxEnd3,indxStart3] = find(dispMat2 <= (maxDispAllowed * tmpFrame));
        %         end
        if size(indxEnd3,1) == 1
            indxEnd3 = indxEnd3';
            indxStart3 = indxStart3';
        end

        %add them to the list of possible pairs
        indxEnd2 = [indxEnd2; endsToConsider(indxEnd3)];
        indxStart2 = [indxStart2; startsToConsider(indxStart3)];
                
    end %(for jFrame = iFrame + 1 : iFrame + timeWindow)
    
end %(for iFrame = 1 : numFrames)

%get total number of pairs
numPairs = length(indxEnd2);

%clear variables from memory
clear dispMat2 maxDispAllowed

%reserve memory for cost matrix vectors
indx1 = zeros(numPairs,1); %row number in cost matrix
indx2 = zeros(numPairs,1); %column number in cost matrix
cost  = zeros(numPairs,1); %cost value

%put time scaling of linear motion in a vector
% timeScalingLin = ones(timeWindow,1);
timeScalingLin = [(1:timeReachConfL).^linScaling(1) ...
    (timeReachConfL)^linScaling(1) * (2:timeWindow-timeReachConfL+1).^linScaling(2)];

%put time scaling of Brownian motion in a vector
% timeScalingBrown = ones(timeWindow,1);
timeScalingBrown = [(1:timeReachConfB).^brownScaling(1) ...
    (timeReachConfB)^brownScaling(1) * (2:timeWindow-timeReachConfB+1).^brownScaling(2)];

% timeGapAll = [];
% timeGapAll2 = [];

%go over all possible pairs of starts and ends
for iPair = 1 : numPairs
    
    %get indices of starts and ends
    iStart = indxStart2(iPair);
    iEnd = indxEnd2(iPair);
    
    %determine the time gap between them
    timeGap = trackStartTime(iStart) - trackEndTime(iEnd);
    %     timeGapAll2 = [timeGapAll2; timeGap];

    %get the types of the two tracks
    trackTypeS = trackType(iStart);
    trackTypeE = trackType(iEnd);
    
    %calculate the vector connecting the end of track iEnd to the
    %start of track iStart and compute its magnitude
    dispVec = coordStart(iStart,:) - coordEnd(iEnd,:);
    dispVecMag = norm(dispVec);

    %determine whether the connecting vector is parallel or anti-parallel
    %to the tracks' directions of motion
    parallelToS = (dispVec * xyzVelS(iStart,:,1)') > 0;
    parallelToE = (dispVec * xyzVelE(iEnd,:,1)') > 0;
    
    %determine the search area of track iStart
    if linearMotion == 1 && ~parallelToS
        longVecS = longRedVecSAll(:,timeGap,iStart);
    else
        longVecS = longVecSAll(:,timeGap,iStart);
    end
    shortVecS = shortVecSAll(:,timeGap,iStart);
        
    %determine the search area of track iEnd
    if linearMotion == 1 && ~parallelToE
        longVecE = longRedVecEAll(:,timeGap,iEnd);
    else
        longVecE = longVecEAll(:,timeGap,iEnd);
    end
    shortVecE = shortVecEAll(:,timeGap,iEnd);

    %calculate the magnitudes of the long and short search vectors
    %of both start and end
    longVecMagS = norm(longVecS);
    shortVecMagS = norm(shortVecS);
    longVecMagE = norm(longVecE);
    shortVecMagE = norm(shortVecE);

    %project the connecting vector onto the long and short vectors
    %of track iStart and take absolute value
    projStartLong = abs(dispVec * longVecS) / longVecMagS;
    projStartShort = abs(dispVec * shortVecS) / shortVecMagS;
    
    %project the connecting vector onto the long and short vectors
    %of track iEnd and take absolute value
    projEndLong = abs(dispVec * longVecE) / longVecMagE;
    projEndShort = abs(dispVec * shortVecE) / shortVecMagE;
    
    %get second short vector and project along it if problem is 3D
    if probDim == 3
        shortVecS3D = shortVecS3DAll(:,timeGap,iStart); %second short vectors
        shortVecE3D = shortVecE3DAll(:,timeGap,iEnd);
        shortVecMagS3D = sqrt(shortVecS3D' * shortVecS3D); %their magnitudes
        shortVecMagE3D = sqrt(shortVecE3D' * shortVecE3D);
        projStartShort3D = abs(dispVec * shortVecS3D) / shortVecMagS3D; %projection
        projEndShort3D = abs(dispVec * shortVecE3D) / shortVecMagE3D;
    else %if problem is 2D, make values zero
        shortVecMagS3D = 0;
        shortVecMagE3D = 0;
        projStartShort3D = 0;
        projEndShort3D = 0;
    end
    
    %calculate the vector connecting the centers of the two tracks
    %     cen2cenVec = trackCenter(iStart,:) - trackCenter(iEnd,:);
    cen2cenVec = dispVec;
    cen2cenVecMag = sqrt(cen2cenVec * cen2cenVec');
    
    %decide whether this is a possible link based on the types of
    %the two tracks
    switch trackTypeE
        
        case 1 %if end is directed
            
            switch trackTypeS
                
                case 1 %if start is directed
                    
                    %calculate the cosine of the angle between velocity
                    %vectors
                    cosAngle = longVecE' * longVecS / (longVecMagE * longVecMagS);
                    
                    %calculate the square sine of the angle between velocity vectors
                    sin2Angle = 1 - cosAngle^2;
                    
                    %calculate the square sine of the angle between each
                    %motion direction vector and the center-to-center vector
                    sin2AngleE = 1 - (cen2cenVec * longVecE / ...
                        (longVecMagE * cen2cenVecMag))^2;
                    sin2AngleS = 1 - (cen2cenVec * longVecS / ...
                        (longVecMagS * cen2cenVecMag))^2;

                    %check whether 
                    %(1) the end of track iEnd is within the search
                    %rectangle of the start of track iStart,
                    %(2) the start of track iStart is within the search
                    %rectangle of the end of track iEnd,
                    %(3) the angle between the two directions of motion
                    %is within acceptable bounds, and 
                    %(4) the angle between directions of motion and vector
                    %connecting end and start is within acceptable bounds
                    possibleLink = ...
                        ((projEndLong <= longVecMagE) && ...
                        (projEndShort <= shortVecMagE) && ...
                        (projEndShort3D <= shortVecMagE3D)) && ...
                        ...
                        ((projStartLong <= longVecMagS) && ...
                        (projStartShort <= shortVecMagS) && ...
                        (projStartShort3D <= shortVecMagS3D)) && ...
                        ...
                        (sin2Angle <= sin2AngleMax) && ...
                        ...
                        ((sin2AngleE <= sin2AngleMaxVD) && ...
                        (sin2AngleS <= sin2AngleMaxVD));
                    
                    %for directed motion without switching, only allow
                    %links between tracks moving in the same direction
                    if linearMotion == 1
                        possibleLink = possibleLink && (cosAngle >= 0);
                    end
                    
                case 0 %if start is Brownian

                    %calculate the square sine of the angle between the
                    %end's motion direction vector and the center-to-center
                    %vector
                    sin2AngleE = 1 - (cen2cenVec * longVecE / ...
                        (longVecMagE * cen2cenVecMag))^2;

                    %check whether
                    %(1) the start of track iStart is within the
                    %search rectangle of the end of track iEnd
                    %(2) the end of track iEnd is within the search disc
                    %of the start of track iStart, and
                    %(3) the angle between end's direction of motion and
                    %vector connecting end and start is within acceptable
                    %bounds
                    possibleLink = ...
                        ((projEndLong <= longVecMagE) && ...
                        (projEndShort <= shortVecMagE) && ...
                        (projEndShort3D <= shortVecMagE3D)) && ...
                        ...
                        (dispVecMag <= longVecMagS) && ...
                        ...
                        (sin2AngleE <= sin2AngleMaxVD);

                otherwise %if start is undetermined

                    %calculate the square sine of the angle between the
                    %end's motion direction vector and the center-to-center
                    %vector
                    sin2AngleE = 1 - (cen2cenVec * longVecE / ...
                        (longVecMagE * cen2cenVecMag))^2;

                    %check whether
                    %(1) the start of track iStart is within the search
                    %rectangle of the end of track iEnd,
                    %(2) the end of track iEnd is within the search disc
                    %of the start of track iStart, and -- NO MORE
                    %(3) the angle between end's direction of motion and
                    %vector connecting end and start is within acceptable
                    %bounds
                    possibleLink = ...
                        ((projEndLong <= longVecMagE) && ...
                        (projEndShort <= shortVecMagE) && ...
                        (projEndShort3D <= shortVecMagE3D)) && ...
                        ...
                        (sin2AngleE <= sin2AngleMaxVD);
                    %                         ...
                    %                         (dispVecMag <= longVecMagS) && ...
                    
            end %(switch trackTypeS)
            
        case 0 %if end is Brownian
            
            switch trackTypeS
                
                case 1 %if start is directed
                    
                    %calculate the square sine of the angle between the
                    %start's motion direction vector and the
                    %center-to-center vector
                    sin2AngleS = 1 - (cen2cenVec * longVecS / ...
                        (longVecMagS * cen2cenVecMag))^2;
                    
                    %check whether
                    %(1) the end of track iEnd is within the search
                    %rectangle of the start of track iStart,
                    %(2) the start of track iStart is within the search
                    %disc of the end of track iEnd, and
                    %(3) the angle between start's direction of motion and
                    %vector connecting end and start is within acceptable
                    %bounds
                    possibleLink = ...
                        (dispVecMag <= longVecMagE) && ...
                        ...
                        ((projStartLong <= longVecMagS) && ...
                        (projStartShort <= shortVecMagS) && ...
                        (projStartShort3D <= shortVecMagS3D)) && ...
                        ...
                        (sin2AngleS <= sin2AngleMaxVD);
                    
                case 0 %if start is Brownian
                    
                    %check whether the end of track iEnd is within the search
                    %disc of the start of track iStart and vice versa
                    possibleLink = ...
                        (dispVecMag <= longVecMagE) && ...
                        (dispVecMag <= longVecMagS);
                    
                otherwise %if start is undetermined
                    
                    %check whether the end of track iEnd is within the search
                    %disc of the start of track iStart and vice versa -- NO
                    %MORE START
                    possibleLink = ...
                        (dispVecMag <= longVecMagE); % && ...
                    %                         (dispVecMag <= longVecMagS);
                    
            end %(switch trackTypeS)
            
        otherwise %if end is undetermined
            
            switch trackTypeS
                
                case 1 %if start is directed
                    
                    %calculate the square sine of the angle between the
                    %start's motion direction vector and the
                    %center-to-center vector
                    sin2AngleS = 1 - (cen2cenVec * longVecS / ...
                        (longVecMagS * cen2cenVecMag))^2;
                    
                    %check whether
                    %(1) the end of track iEnd is within the search
                    %rectangle of the start of track iStart,
                    %(2) the start of track iStart is within the search
                    %disc of the end of track iEnd, and -- NO MORE
                    %(3) the angle between start's direction of motion and
                    %vector connecting end and start is within acceptable
                    %bounds
                    possibleLink = ...
                        ((projStartLong <= longVecMagS) && ...
                        (projStartShort <= shortVecMagS) && ...
                        (projStartShort3D <= shortVecMagS3D)) && ...
                        ...
                        (sin2AngleS <= sin2AngleMaxVD);
                    %                         (dispVecMag <= longVecMagE) && ...
                    %                         ...
                    
                case 0 %if start is Brownian
                    
                    %check whether the end of track iEnd is within the search
                    %disc of the start of track iStart and vice versa -- NO
                    %MORE END
                    possibleLink = ...
                        (dispVecMag <= longVecMagS);
                    %                         (dispVecMag <= longVecMagE) && ...
                    
                otherwise %if start is undetermined
                    
                    %check whether the end of track iEnd is within the search
                    %disc of the start of track iStart and vice versa
                    possibleLink = ...
                        (dispVecMag <= longVecMagE) && ...
                        (dispVecMag <= longVecMagS);
                    
            end %(switch trackTypeS)
            
    end %(switch trackTypeE)
    
    %if this is a possible link ...
    if possibleLink
        
        %calculate the average displacement for the two tracks combined
        meanDispTrack1 = trackMeanDispMag(iStart);
        meanDispTrack1(isnan(meanDispTrack1)) = meanDispAllTracks;
        meanDispTrack2 = trackMeanDispMag(iEnd);
        meanDispTrack2(isnan(meanDispTrack2)) = meanDispAllTracks;
        meanDisp2Tracks = mean([meanDispTrack1 meanDispTrack2]);
        
        %calculate the cost of linking
        dispVecMag2 = dispVecMag ^ 2;
        if trackTypeE == 1 && trackTypeS == 1
            cost12 = dispVecMag2 * (1 + mean([sin2Angle sin2AngleE sin2AngleS])) ...
                / (timeScalingLin(timeGap) * meanDisp2Tracks)^2;
            %                 / (meanDisp2Tracks)^2;
        elseif trackTypeE == 1
            cost12 = dispVecMag2 * (1 + sin2AngleE) ...
                / (mean([timeScalingLin(timeGap)*meanDispTrack2 ...
                timeScalingBrown(timeGap)*meanDispTrack1]))^2;
            %                 / (meanDisp2Tracks)^2;
        elseif trackTypeS == 1
            cost12 = dispVecMag2 * (1 + sin2AngleS) ...
                / (mean([timeScalingLin(timeGap)*meanDispTrack1 ...
                timeScalingBrown(timeGap)*meanDispTrack2]))^2;
            %                 / (meanDisp2Tracks)^2;
        else
            cost12 = dispVecMag2 ...
                / (timeScalingBrown(timeGap) * meanDisp2Tracks)^2;
            %                 / (meanDisp2Tracks)^2;
        end
        
        %penalize cost for lifetime considerations
        if ~isempty(lftCdf)
            cost12 = cost12 / oneMinusLftCdf(trackEndTime(iStart)-trackStartTime(iEnd)+2);
        end
        
        %if the lifetime consideration does not make this link impossible
        if isfinite(cost12)
            
            %penalize cost for gap length considerations
            %by using timeGap - 2, there is no penalty for things that disappear for one frame and then come back
            %this puts those gaps at the same cost level as merging and splitting
            cost12 = cost12 * gapPenalty^(timeGap-2);
            
            %store time gap
            %             timeGapAll = [timeGapAll; timeGap];
            
            %add this cost to the list of costs
            cost(iPair) = cost12;
            
            %specify the location of this pair in the cost matrix
            indx1(iPair) = iEnd; %row number
            indx2(iPair) = iStart; %column number
            
        end
        
    end %(if possibleLink)
    
end %(for iPair = 1 : numPairs)

%keep only pairs that turned out to be possible
possiblePairs = find(indx1 ~= 0);
indx1 = indx1(possiblePairs);
indx2 = indx2(possiblePairs);
cost  = cost(possiblePairs);

clear possiblePairs

%% Merging and splitting

%define some merging and splitting variables
numMerge  =  0; %index counting merging events
indxMerge = []; %vector storing merging track number
altCostMerge = []; %vector storing alternative costs of not merging
numSplit  =  0; %index counting splitting events
indxSplit = []; %vector storing splitting track number
altCostSplit = []; %vector storing alternative costs of not splitting

%if merging and/or splitting are to be considered ...
if mergeSplit > 0

    %get the absolute upper limit of acceptable displacements in one frame
    %as the maximum of (maximum velocity multiplied by probDim*linStdMult(1),
    %maxSearchRadiusCG)
    maxDispAllowed = max( max( abs([xyzVelS(:);xyzVelE(:)]) ) * probDim * linStdMult(1), ...
        maxSearchRadius );
    maxDispAllowed = max(maxDispAllowed,resLimit);
    
    %define number of time points before/after current one to use in
    %calculating mean amplitude to evaluate merges and splits (1 = this
    %point + 1 before/after, 2 = this point + 2 before/after, etc.)
    nTpMS = 2;

    %costs of merging
    if mergeSplit == 1 || mergeSplit == 2

        %go over all track end times
        for endTime = 1 : numFrames-1

            %find tracks that end in this frame
            endsToConsider = tracksPerFrame(endTime).ends;
            
            %if requested, remove tracks that have a gap closing 
            %possibility - no merging allowed in this case
            if gapExcludeMS
                endsToConsider = setdiff(endsToConsider,indx1);
            end

            %find tracks that start before or in this frame and end after this
            %frame
            mergesToConsider = intersect(vertcat(tracksPerFrame(1:endTime).starts),...
                vertcat(tracksPerFrame(endTime+1:end).ends));

            %get index indicating frame of merging
            timeIndx  = endTime*8;

            %calculate displacement between track ends and other tracks in the
            %next frame
            dispMat2 = createDistanceMatrix(coordEnd(endsToConsider,:), ...
                full(trackedFeatInfo(mergesToConsider,timeIndx+1:timeIndx+probDim)));

            %find possible pairs
            [indxEnd2,indxMerge2] = find(dispMat2 <= maxDispAllowed);
            numPairs = length(indxEnd2);

            %clear memory
            clear dispMat2

            %map from indices to track indices
            indxEnd2 = endsToConsider(indxEnd2);
            indxMerge2 = mergesToConsider(indxMerge2);

            %reserve memory for cost vectors and related vectors
            indx1MS   = zeros(numPairs,1);
            indx2MS   = zeros(numPairs,1);
            costMS    = zeros(numPairs,1);
            altCostMS = zeros(numPairs,1);
            indxMSMS  = zeros(numPairs,1);

            %go over all possible pairs
            for iPair = 1 : numPairs

                %get indices of ending track and track it might merge with
                iEnd = indxEnd2(iPair);
                iMerge = indxMerge2(iPair);

                %calculate the vector connecting the end of track iEnd to the
                %point of merging and compute its magnitude
                dispVec = full(trackedFeatInfo(iMerge,timeIndx+1:timeIndx+probDim)) ...
                    - coordEnd(iEnd,:);
                dispVecMag = sqrt(dispVec * dispVec');

                %determine whether the connecting vector is parallel or anti-parallel
                %to the ending track's direction of motion
                parallelToE = (dispVec * xyzVelE(iEnd,:,1)') > 0;
                
                %determine the search area of track iEnd
                if linearMotion == 1 && ~parallelToE
                    longVecE = longRedVecEAllMS(:,1,iEnd);
                else
                    longVecE = longVecEAllMS(:,1,iEnd);
                end
                shortVecE = shortVecEAllMS(:,1,iEnd);

                %calculate the magnitudes of the long and short search vectors
                %of the end
                longVecMagE = sqrt(longVecE' * longVecE);
                shortVecMagE = sqrt(shortVecE' * shortVecE);

                %project the connecting vector onto the long and short vectors
                %of track iEnd and take absolute value
                projEndLong = abs(dispVec * longVecE) / longVecMagE;
                projEndShort = abs(dispVec * shortVecE) / shortVecMagE;

                %get the amplitude of track iEnd before its end - take the
                %last nTpMS+1 points
                indxBefore = 8*(endTime-1)+4 - 8*(0:nTpMS);
                indxBefore = indxBefore(indxBefore > 1);
                ampE = full(trackedFeatInfo(iEnd,indxBefore));
                ampE = mean(ampE(ampE~=0));

                %get the amplitude of the merging track before and after
                %merging - take nTpMS+1 points on each side
                ampM1 = full(trackedFeatInfo(iMerge,indxBefore)); %before merging
                ampM1 = mean(ampM1(ampM1~=0));
                indxAfter = 8*endTime+4 + 8*(0:nTpMS);
                indxAfter = indxAfter(indxAfter < 8*numFrames);
                ampM = full(trackedFeatInfo(iMerge,indxAfter)); %after merging
                ampM = mean(ampM(ampM~=0));

                %calculate the ratio of the amplitude after merging to the sum
                %of the amplitudes before merging
                ampRatio = ampM / (ampE + ampM1);
                
                %calculate the individual ratios
                ampRatioIndME = ampM / ampE;
                ampRatioIndMM1 = ampM / ampM1;
                
                %if amplitude is not to be used, give amplitude-related
                %variables dummy values
                if ~useAmp
                    [ampRatio,ampM,ampM1] = deal(1);
                    [ampRatioIndME,ampRatioIndMM1] = deal(1.1);
                end

                %decide whether this is a possible link based on displacement,
                %directionality and amplitude ratio
                if trackType(iEnd) == 1 %if ending track is linear

                    %get second short vector and project along it if problem is 3D
                    if probDim == 3
                        shortVecE3D = shortVecE3DAllMS(:,1,iEnd);
                        shortVecMagE3D = sqrt(shortVecE3D' * shortVecE3D);
                        projEndShort3D = abs(dispVec * shortVecE3D) / shortVecMagE3D;
                    else %if problem is 2D, make values zero
                        shortVecMagE3D = 0;
                        projEndShort3D = 0;
                    end

                    %calculate the vector connecting the centers of the two tracks
                    %                     cen2cenVec = trackCenter(iEnd,:) - trackCenter(iMerge,:);
                    cen2cenVec = dispVec;
                    cen2cenVecMag = sqrt(cen2cenVec * cen2cenVec');

                    %calculate the square sine of the angle between the
                    %end's motion direction vector and the center-to-center vector
                    sin2AngleE = 1 - (cen2cenVec * longVecE / ...
                        (longVecMagE * cen2cenVecMag))^2;

                    %check whether ...
                    %(1) the feature to be merged with is within the search
                    %region of the end of track iEnd
                    %(2) the center-to-center vector is reasonably well
                    %aligned with the directionality of track iEnd
                    %(3) the amplitudes satisfy the following: 
                    %   (i) ampRatio is within acceptable limits,
                    %   (ii) the intensity after merging is larger than
                    %the individual intensities before merging,
                    %   (iii) ampRatio is closer to 1 than ampRatioIndMM1
                    %(i.e. the ratio were the track to be merged with left
                    %alone)
                    possibleLink = projEndLong <= longVecMagE && ...
                        projEndShort <= shortVecMagE && ...
                        projEndShort3D <= shortVecMagE3D && ...
                        sin2AngleE <= sin2AngleMaxVD && ...
                        ampRatio >= minAmpRatio && ampRatio <= maxAmpRatio && ...
                        ampRatioIndME > 1 && ampRatioIndMM1 > 1 && ...
                        abs(ampRatio-1) < abs(ampRatioIndMM1-1);

                else %if ending track is Brownian or undetermined

                    %assign the dummy value of zero to sin2AngleE
                    sin2AngleE = 0;

                    %look at displacement and amplitudes only (no
                    %directionality)
                    possibleLink = dispVecMag <= longVecMagE && ...
                        ampRatio >= minAmpRatio && ampRatio <= maxAmpRatio && ...
                        ampRatioIndME > 1 && ampRatioIndMM1 > 1 && ...
                        abs(ampRatio-1) < abs(ampRatioIndMM1-1);

                end

                %if this is a possible link ...
                if possibleLink

                    %calculate the cost of linking
                    dispVecMag2 = dispVecMag ^ 2; %due to displacement
                    ampCost = ampRatio; %due to amplitude
                    ampCost(ampCost<1) = ampCost(ampCost<1) ^ (-2); %punishment harsher when intensity of merged feature < sum of intensities of merging features
                    meanDisp2Tracks = trackMeanDispMag(iEnd); %for displacement scaling
                    if isnan(meanDisp2Tracks)
                        meanDisp2Tracks = meanDispAllTracks;
                    end
                    cost12 = dispVecMag2 * ampCost * (1 + sin2AngleE) ...
                        / (meanDisp2Tracks^2); %cost

                    %penalize cost for lifetime considerations
                    if ~isempty(lftCdf)
                        cost12 = cost12 / oneMinusLftCdf(trackEndTime(iMerge)-trackStartTime(iEnd)+2);
                    end
                    
                    %if the lifetime consideration does not make this link impossible
                    if ~isinf(cost12)
                        
                        %add this cost to the list of costs
                        costMS(iPair) = cost12;
                        
                        %check whether the track being merged with has had
                        %something possibly merging with it in this same frame
                        prevAppearance = find(indxMSMS == iMerge);

                        %if this track in this frame did not appear before ...
                        %THIS SECTION IS OUTDATED
                        if isempty(prevAppearance)

                            %increase the "merge index" by one
                            numMerge = numMerge + 1;

                            %save the merging track's index
                            indxMSMS(iPair) = iMerge;

                            %store the location of this pair in the cost matrix
                            indx1MS(iPair) = iEnd; %row number
                            indx2MS(iPair) = numMerge+numTracks; %column number

                            %calculate the alternative cost of not merging for the
                            %track that the end is possibly merging with

                            %get the average square displacement in this track
                            trackCoord = trackedFeatInfo(indxMSMS(iPair),:);
                            trackCoord = reshape(trackCoord',8,[]);
                            if issparse(trackCoord)
                                trackCoord = full(trackCoord);
                                trackCoord(trackCoord==0) = NaN;
                                if probDim == 2
                                    trackCoord(3,:) = 0;
                                    trackCoord(7,:) = 0;
                                end
                            end
                            dispVecMag2 = (diff(trackCoord,1,2)).^2;
                            dispVecMag2 = nanmean(dispVecMag2,2);
                            dispVecMag2 = sum(dispVecMag2(1:probDim));
                            
                            %if the average square displacement is smaller
                            %than resLimit^2, then expand it
                            dispVecMag2 = max([dispVecMag2 resLimit^2]);

                            %calculate intensity cost if no merge happens
                            ampCost = ampM / ampM1;
                            ampCost(ampCost<1) = ampCost(ampCost<1) ^ (-2);
                            
                            %get track's mean displacement
                            meanDisp1Track = trackMeanDispMag(indxMSMS(iPair));
                            meanDisp1Track(isnan(meanDisp1Track)) = ...
                                meanDispAllTracks;
                            
                            %calculate alternative cost
                            cost12 = dispVecMag2 * ampCost ...
                                / (meanDisp1Track^2);

                            %although alternative cost is still calculated,
                            %it is actually not used any more in the end
                            
                            %add this cost to the list of alternative costs
                            altCostMS(iPair) = cost12;

                        else %if this track in this frame appeared before

                            %do not increase the "merge index" or save the merging
                            %track's index (they are already saved)

                            %store the location of this pair in the cost matrix
                            indx1MS(iPair) = iEnd; %row number
                            indx2MS(iPair) = indx2MS(prevAppearance); %column number

                            %no need to calculate and save the alternative cost
                            %since that is already saved from previous encounter

                        end %(if isempty(prevAppearance))

                    end %(if ~isinf(cost12))

                end %(if possibleLink)

            end %(for iPair = 1 : numPairs)

            %keep only pairs that turned out to be possible
            possiblePairs = find(indx1MS ~= 0);
            indx1MS   = indx1MS(possiblePairs);
            indx2MS   = indx2MS(possiblePairs);
            costMS    = costMS(possiblePairs);
            possibleMerges = find(indxMSMS ~= 0);
            indxMSMS  = indxMSMS(possibleMerges);
            altCostMS = altCostMS(possibleMerges);
            clear possiblePairs possibleMerges

            %append these vectors to overall cost vector and related vectors
            indx1 = [indx1; indx1MS];
            indx2 = [indx2; indx2MS];
            cost  = [cost; costMS];
            altCostMerge = [altCostMerge; altCostMS];
            indxMerge = [indxMerge; indxMSMS];

        end %(for endTime = 1 : numFrames-1)
        
    end %(if mergeSplit == 1 || mergeSplit == 2)

    %costs of splitting
    if mergeSplit == 1 || mergeSplit == 3

        %go over all track starting times
        for startTime = 2 : numFrames

            %find tracks that start in this frame
            startsToConsider = tracksPerFrame(startTime).starts;
            
            %if requested, remove tracks that have a gap closing 
            %possibility - no splits allowed in this case
            if gapExcludeMS
                startsToConsider = setdiff(startsToConsider,indx2);
            end

            %find tracks that start before this frame and end after or in this frame
            splitsToConsider = intersect(vertcat(tracksPerFrame(1:startTime-1).starts),...
                vertcat(tracksPerFrame(startTime:end).ends));

            %get index indicating time of splitting
            timeIndx  = (startTime-2)*8;

            %calculate displacement between track starts and other tracks in the
            %previous frame
            dispMat2 = createDistanceMatrix(coordStart(startsToConsider,:), ...
                full(trackedFeatInfo(splitsToConsider,timeIndx+1:timeIndx+probDim)));

            %find possible pairs
            [indxStart2,indxSplit2] = find(dispMat2 <= maxDispAllowed);
            numPairs = length(indxStart2);

            %clear memory
            clear dispMat2

            %map from indices to track indices
            indxStart2 = startsToConsider(indxStart2);
            indxSplit2 = splitsToConsider(indxSplit2);

            %reserve memory for cost vectors and related vectors
            indx1MS   = zeros(numPairs,1);
            indx2MS   = zeros(numPairs,1);
            costMS    = zeros(numPairs,1);
            altCostMS = zeros(numPairs,1);
            indxMSMS  = zeros(numPairs,1);

            %go over all possible pairs
            for iPair = 1 : numPairs

                %get indices of starting track and track it might have split from
                iStart = indxStart2(iPair);
                iSplit = indxSplit2(iPair);

                %calculate the vector connecting the end of track iStart to the
                %point of splitting and compute its magnitude
                dispVec = coordStart(iStart,:) - full(trackedFeatInfo(iSplit,...
                    timeIndx+1:timeIndx+probDim));
                dispVecMag = sqrt(dispVec * dispVec');

                %determine whether the connecting vector is parallel or anti-parallel
                %to the starting track's direction of motion
                parallelToS = (dispVec * xyzVelS(iStart,:,1)') > 0;
                
                %determine the search area of track iStart
                if linearMotion == 1 && ~parallelToS
                    longVecS = longRedVecSAllMS(:,1,iStart);
                else
                    longVecS = longVecSAllMS(:,1,iStart);
                end
                shortVecS = shortVecSAllMS(:,1,iStart);

                %calculate the magnitudes of the long and short search vectors
                %of the start
                longVecMagS = sqrt(longVecS' * longVecS);
                shortVecMagS = sqrt(shortVecS' * shortVecS);

                %project the connecting vector onto the long and short vectors
                %of track iStart and take absolute value
                projStartLong = abs(dispVec * longVecS) / longVecMagS;
                projStartShort = abs(dispVec * shortVecS) / shortVecMagS;
                
                %get the amplitude of track iStart after its start - take
                %the first nTpMS+1 points
                indxAfter = 8*(startTime-1)+4 + 8*(0:nTpMS);
                indxAfter = indxAfter(indxAfter < 8*numFrames);
                ampS = full(trackedFeatInfo(iStart,indxAfter));
                ampS = mean(ampS(ampS~=0));

                %get the amplitude of the splitting track after and before
                %splitting - take nTpMS+1 points on each side
                ampSp1 = full(trackedFeatInfo(iSplit,indxAfter)); %after splitting
                ampSp1 = mean(ampSp1(ampSp1~=0));
                indxBefore = 8*(startTime-2)+4 - 8*(0:nTpMS);
                indxBefore = indxBefore(indxBefore > 1);
                ampSp = full(trackedFeatInfo(iSplit,indxBefore)); %before splitting
                ampSp = mean(ampSp(ampSp~=0));

                %calculate the ratio of the amplitude before splitting to the sum
                %of the amplitudes after splitting
                ampRatio = ampSp / (ampS + ampSp1);

                %calculate the individual ratios
                ampRatioIndSpS = ampSp / ampS;
                ampRatioIndSpSp1 = ampSp / ampSp1;
                
                %if amplitude is not to be used, give amplitude-related
                %variables dummy values
                if ~useAmp
                    [ampRatio,ampSp,ampSp1] = deal(1);
                    [ampRatioIndSpS,ampRatioIndSpSp1] = deal(1.1);
                end

                %decide whether this is a possible link based on displacement,
                %directionality and amplitude ratio
                if trackType(iStart) == 1

                    %get second short vector and project along it if problem is 3D
                    if probDim == 3
                        shortVecS3D = shortVecS3DAllMS(:,1,iStart); %second short vectors
                        shortVecMagS3D = sqrt(shortVecS3D' * shortVecS3D); %their magnitudes
                        projStartShort3D = abs(dispVec * shortVecS3D) / shortVecMagS3D; %projection
                    else %if problem is 2D, make values zero
                        shortVecMagS3D = 0;
                        projStartShort3D = 0;
                    end

                    %calculate the vector connecting the centers of the two tracks
                    %                     cen2cenVec = trackCenter(iStart,:) - trackCenter(iSplit,:);
                    cen2cenVec = dispVec;
                    cen2cenVecMag = sqrt(cen2cenVec * cen2cenVec');

                    %calculate the square sine of the angle between the
                    %end's motion direction vector and the center-to-center vector
                    sin2AngleS = 1 - (cen2cenVec * longVecS / ...
                        (longVecMagS * cen2cenVecMag))^2;

                    %check whether ...
                    %(1) the feature to be split from is within the search
                    %region of the start of track iStart
                    %(2) the center-to-center vector is reasonably well
                    %aligned with the directionality of track iStart
                    %(3) the amplitudes satisfy the following: 
                    %   (i) ampRatio is within acceptable limits,
                    %   (ii) the intensity before splitting is larger than
                    %the individual intensities after splitting,
                    %   (iii) ampRatio is closer to 1 than ampRatioIndSpSp1
                    %(i.e. the ratio were the track to be splitted from left
                    %alone)
                    possibleLink = projStartLong <= longVecMagS && ...
                        projStartShort <= shortVecMagS && ...
                        projStartShort3D <= shortVecMagS3D && ...
                        sin2AngleS <= sin2AngleMaxVD && ...
                        ampRatio >= minAmpRatio && ampRatio <= maxAmpRatio && ...
                        ampRatioIndSpS > 1 && ampRatioIndSpSp1 > 1 && ...
                        abs(ampRatio-1) < abs(ampRatioIndSpSp1-1);
                    
                else %if starting track is Brownian or undetermined

                    %assign the dummy value of zero to sin2AngleS
                    sin2AngleS = 0;

                    %look at displacement and amplitude ratio only (no
                    %directionality)
                    possibleLink = dispVecMag <= longVecMagS && ...
                        ampRatio >= minAmpRatio && ampRatio <= maxAmpRatio && ...
                        ampRatioIndSpS > 1 && ampRatioIndSpSp1 > 1 && ...
                        abs(ampRatio-1) < abs(ampRatioIndSpSp1-1);

                end

                %if this is a possible link ...
                if possibleLink

                    %calculate the cost of linking
                    dispVecMag2 = dispVecMag ^ 2; %due to displacement
                    ampCost = ampRatio; %due to amplitude
                    ampCost(ampCost<1) = ampCost(ampCost<1) ^ (-2); %punishment harsher when intensity of splitting feature < sum of intensities of features after splitting
                    meanDisp2Tracks = trackMeanDispMag(iStart); %for displacement scaling
                    if isnan(meanDisp2Tracks)
                        meanDisp2Tracks = meanDispAllTracks;
                    end
                    cost12 = dispVecMag2 * ampCost * (1 + sin2AngleS) ...
                        / (meanDisp2Tracks^2);

                    %penalize cost for lifetime considerations
                    if ~isempty(lftCdf)
                        cost12 = cost12 / oneMinusLftCdf(trackEndTime(iStart)-trackStartTime(iSplit)+2);
                    end
                    
                    %if the lifetime consideration does not make this link impossible
                    if ~isinf(cost12)
                        
                        %add this cost to the list of costs
                        costMS(iPair) = cost12;
                        
                        %check whether the track being split from has had something
                        %possibly splitting from it in this same frame
                        prevAppearance = find(indxMSMS == iSplit);

                        %if this track in this frame did not appear before ...
                        %THIS SECTION IS OUTDATED
                        if isempty(prevAppearance)
                            
                            %increase the "split index" by one
                            numSplit = numSplit + 1;
                            
                            %save the splitting track's number
                            indxMSMS(iPair) = iSplit;

                            %store the location of this pair in the cost matrix
                            indx1MS(iPair) = numSplit+numTracks; %row number
                            indx2MS(iPair) = iStart; %column number

                            %calculate the alternative cost of not splitting for the
                            %track that the start is possibly splitting from

                            %get the average square displacement in this track
                            trackCoord = trackedFeatInfo(indxMSMS(iPair),:);
                            trackCoord = reshape(trackCoord',8,[]);
                            if issparse(trackCoord)
                                trackCoord = full(trackCoord);
                                trackCoord(trackCoord==0) = NaN;
                                if probDim == 2
                                    trackCoord(3,:) = 0;
                                    trackCoord(7,:) = 0;
                                end
                            end
                            dispVecMag2 = (diff(trackCoord,1,2)).^2;
                            dispVecMag2 = nanmean(dispVecMag2,2);
                            dispVecMag2 = sum(dispVecMag2(1:probDim));

                            %if the average square displacement is smaller
                            %than resLimit^2, then expand it
                            dispVecMag2 = max([dispVecMag2 resLimit^2]);

                            %calculate intensity cost if no split happens
                            ampCost = ampSp / ampSp1;
                            ampCost(ampCost<1) = ampCost(ampCost<1) ^ (-2);

                            %get track's mean displacement
                            meanDisp1Track = trackMeanDispMag(indxMSMS(iPair));
                            meanDisp1Track(isnan(meanDisp1Track)) = ...
                                meanDispAllTracks;
                            
                            %calculate alternative cost
                            cost12 = dispVecMag2 * ampCost ...
                                / (meanDisp1Track^2);

                            %although alternative cost is still calculated,
                            %it is actually not used any more in the end
                            
                            %add this cost to the list of alternative costs
                            altCostMS(iPair) = cost12;

                        else %if this track in this frame appeared before

                            %do not increase the "split index" or save the
                            %splitting track's index (they are already saved)

                            %store the location of this pair in the cost matrix
                            indx1MS(iPair) = indx1MS(prevAppearance); %row number
                            indx2MS(iPair) = iStart; %column number

                            %no need to calculate and save the alternative cost
                            %since that is already saved from previous appearance

                        end %(if isempty(prevAppearance))

                    end %(if ~isinf(cost12))

                end %(if possibleLink)

            end %(for for iPair = 1 : numPairs)

            %keep only pairs that turned out to be possible
            possiblePairs = find(indx1MS ~= 0);
            indx1MS   = indx1MS(possiblePairs);
            indx2MS   = indx2MS(possiblePairs);
            costMS    = costMS(possiblePairs);
            possibleSplits = find(indxMSMS ~= 0);
            altCostMS = altCostMS(possibleSplits);
            indxMSMS  = indxMSMS(possibleSplits);
            clear possiblePairs possibleSplits
            
            %append these vectors to overall cost and related vectors
            indx1 = [indx1; indx1MS];
            indx2 = [indx2; indx2MS];
            cost  = [cost; costMS];
            altCostSplit = [altCostSplit; altCostMS];
            indxSplit = [indxSplit; indxMSMS];

        end %(for startTime = 2 : numFrames)
        
    end %(if mergeSplit == 1 || mergeSplit == 3)

end %(if mergeSplit)

%create cost matrix without births and deaths
numEndSplit = numTracks + numSplit;
numStartMerge = numTracks + numMerge;
costMat = sparse(indx1,indx2,cost,numEndSplit,numStartMerge);

%% Append cost matrix to allow births and deaths ...

%determine the cost of birth and death
if strategyBD == -1 %value determined from cost definitions
    
    costBD = 9 ... %maximum expected contribution to cost from [displacement/mean(displacement)]^2
        * 2 ... %maximum expected contribution to cost from angle penalty
        * gapPenalty^(timeWindow-2); %maximum expected contribution to cost from gap penalty
    
    if mergeSplit > 0 && useAmp %maximum contribution to cost from amplitude ratio when merging & splitting
        if minAmpRatio > 0 && maxAmpRatio < Inf
            costBDTmp = 9 * 2 * max(maxAmpRatio,minAmpRatio^-2);
        else
            disp('With chosen parameters, there are no bounds on amplitude ratio for merging and splitting. Thus, birth and death cost cannot be calculated from cost definitions, and will be calculated using an automatically-determine percentile.');
            strategyBD = 0;
            costBDTmp = NaN;
        end
        costBD = max(costBD,costBDTmp);
    end
    
end

if strategyBD == 0 %automatically-determined percentile
    
    tmp = (costMat~=0);
    numPotAssignRow = full(sum(tmp,2));
    numPotAssignCol = full(sum(tmp)');
    numPotAssignColAll = sum(numPotAssignCol);
    numPotAssignRowAll = sum(numPotAssignRow);
    numPartCol = length(numPotAssignCol) * 2;
    extraCol = (numPotAssignColAll-numPartCol)/numPotAssignColAll;
    numPartRow = length(numPotAssignRow) * 2;
    extraRow = (numPotAssignRowAll-numPartRow)/numPotAssignRowAll;
    prctile2use = min(100, 100 - mean([extraRow extraCol])*100);
    costBD = 1.05*prctile(cost(:),prctile2use);
    
elseif strategyBD > 0 %user-defined percentile
    
    costBD = 1.05*prctile(cost(:),strategyBD);
        
end

%get the cost for the lower right block
% costLR = min(min(min(costMat))-1,-1);
costLR = costBD;

%create cost matrix that allows for births and deaths
% costMat = [costMat ... %costs for links (gap closing + merge/split)
%     spdiags([costBD*ones(numTracks,1); altCostSplit],0,numEndSplit,numEndSplit); ... %costs for death
%     spdiags([costBD*ones(numTracks,1); altCostMerge],0,numStartMerge,numStartMerge) ...  %costs for birth
%     sparse(indx2,indx1,costLR*ones(length(indx1),1),numStartMerge,numEndSplit)]; %dummy costs to complete the cost matrix

costMat = [costMat ... %costs for links (gap closing + merge/split)
    spdiags(costBD*ones(numTracks+numSplit,1),0,numEndSplit,numEndSplit); ... %costs for death
    spdiags(costBD*ones(numTracks+numMerge,1),0,numStartMerge,numStartMerge) ...  %costs for birth
    sparse(indx2,indx1,costLR*ones(length(indx1),1),numStartMerge,numEndSplit)]; %dummy costs to complete the cost matrix

%determine the nonlinkMarker
nonlinkMarker = min(floor(full(min(min(costMat))))-5,-5);


%% ~~~ the end ~~~
