function [costMat,nonlinkMarker,indxMerge,numMerge,indxSplit,numSplit,...
    errFlag] = plusTipCostMatCloseGaps(trackedFeatInfo,...
    trackedFeatIndx,trackStartTime,trackEndTime,costMatParam,gapCloseParam,...
    kalmanFilterInfo,nnDistLinkedFeat,probDim,movieInfo)
%COSTMATLINEARMOTIONCLOSEGAPS provides a cost matrix for closing gaps and capturing merges/splits using Kalman filter information
%
%SYNOPSIS [costMat,nonlinkMarker,indxMerge,numMerge,indxSplit,numSplit,...
%    errFlag] = costMatLinearMotionCloseGaps(trackedFeatInfo,...
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
%          .fluctRad             : Size in pixels of tube radius around track
%                                trajectory. The search region in the backward
%                                direction will expand out from the track at
%                                the final point.  This value also determines
%                                the search radius around the final point
%                                wherein any candidates will be considered for
%                                forward linking, even if they fall slightly
%                                behind the point.  This ensures that tracks
%                                starting from a fluctuation during a pause
%                                will still be picked up as candidates for
%                                pause.
%          .maxFAngle          : Max angle in degrees allowed between the end
%                                track's final velocity vector and the
%                                displacement vector between end and start.
%                                Also the max angle between the end and start
%                                tracks themselves.
%          .maxBAngle          : Angle in degrees used to expand backward
%                                search region, giving a distance-dependent
%                                criterion for how far a start track
%                                can be from the lattice to be considered a
%                                candidate for linking. THIS IS CURRENTLY A
%                                HARDWIRED PARAMETER
%          .backVelMultFactor : Muliplication factor of max growth speed used
%                               to define candidate search area in the
%                               backward direction.
%       gapCloseParam  : Structure containing variables needed for gap closing.
%                        Contains the fields:
%             .timeWindow : Largest time gap between the end of a track and the
%                           beginning of another that could be connected to
%                           it.
%             .mergeSplit : Logical variable with value 1 if the merging
%                           and splitting of trajectories are to be consided;
%                           and 0 if merging and splitting are not allowed.
%                           For MT tracking, there are no merges/splits, so
%                           this should be 0.
%       kalmanFilterInfo: Structure array with number of entries equal to
%                         number of frames in movie. Contains the fields:
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
%Kathryn Applegate, 2009
%Philippe Roudot 2014: 3D support
%Adapted from costMatCloseGaps.m by Khuloud Jaqaman, April 2007
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
if nargin ~= nargin('plusTipCostMatCloseGaps')
    disp('--plusTipCostMatCloseGaps: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

doTest=0;
doPlot=0;

% get user-set parameters
fluctRad=costMatParam.fluctRad;
maxFAngle = costMatParam.maxFAngle*pi/180;
maxBAngle = costMatParam.maxBAngle*pi/180;
backVelMultFactor = costMatParam.backVelMultFactor;
tMax=gapCloseParam.timeWindow;

%find the number of tracks to be linked and the number of frames in the movie
[nTracks,nFrames] = size(trackedFeatInfo);
nFrames = nFrames / 8;

%list the tracks that start and end in each frame
tracksPerFrame = repmat(struct('starts',[],'ends',[]),nFrames,1);
for iFrame = 1 : nFrames
    tracksPerFrame(iFrame).starts = find(trackStartTime == iFrame); %starts
    tracksPerFrame(iFrame).ends = find(trackEndTime == iFrame); %ends
end

%% Gap closing

% extract feature positions and velocity components
px=trackedFeatInfo(:,1:8:end);
py=trackedFeatInfo(:,2:8:end);
pz=trackedFeatInfo(:,3:8:end);

vx=diff(px,1,2);
vy=diff(py,1,2);
vz=diff(pz,1,2);

% all the instantaneous velocities from all tracks in one vector
velInst=sqrt(vx.^2+vy.^2+vz.^2);
velInst=velInst(:);
velInst(isnan(velInst))=[];

% TRACK STARTS
trackStartPxyzVxyz = zeros(nTracks,6);
% x and y coordinates of the track's first point
trackStartPxyzVxyz(:,1)=cell2mat(arrayfun(@(i) px(i,find(~isnan(px(i,:)),1,'first')),[1:nTracks]','UniformOutput',0));
trackStartPxyzVxyz(:,2)=cell2mat(arrayfun(@(i) py(i,find(~isnan(py(i,:)),1,'first')),[1:nTracks]','UniformOutput',0));
trackStartPxyzVxyz(:,3)=cell2mat(arrayfun(@(i) pz(i,find(~isnan(pz(i,:)),1,'first')),[1:nTracks]','UniformOutput',0));

% average of first three velocity vectors (made from last 4 points
% on track, if that many exist), x and y components
trackStartPxyzVxyz(:,4)=cell2mat(arrayfun(@(i) mean(vx(i,find(~isnan(vx(i,:)),3,'first'))),[1:nTracks]','UniformOutput',0));
trackStartPxyzVxyz(:,5)=cell2mat(arrayfun(@(i) mean(vy(i,find(~isnan(vy(i,:)),3,'first'))),[1:nTracks]','UniformOutput',0));
trackStartPxyzVxyz(:,6)=cell2mat(arrayfun(@(i) mean(vz(i,find(~isnan(vz(i,:)),3,'first'))),[1:nTracks]','UniformOutput',0));

% TRACK ENDS
trackEndPxyzVxyz = zeros(nTracks,6);
% x and y coordinates of the track's last point
trackEndPxyzVxyz(:,1)=cell2mat(arrayfun(@(i) px(i,find(~isnan(px(i,:)),1,'last')),[1:nTracks]','UniformOutput',0));
trackEndPxyzVxyz(:,2)=cell2mat(arrayfun(@(i) py(i,find(~isnan(py(i,:)),1,'last')),[1:nTracks]','UniformOutput',0));
trackEndPxyzVxyz(:,3)=cell2mat(arrayfun(@(i) pz(i,find(~isnan(pz(i,:)),1,'last')),[1:nTracks]','UniformOutput',0));

% average of last three velocity vectors (made from last 4 points
% on track, if that many exist), x and y components
trackEndPxyzVxyz(:,4)=cell2mat(arrayfun(@(i) mean(vx(i,find(~isnan(vx(i,:)),3,'last'))),[1:nTracks]','UniformOutput',0));
trackEndPxyzVxyz(:,5)=cell2mat(arrayfun(@(i) mean(vy(i,find(~isnan(vy(i,:)),3,'last'))),[1:nTracks]','UniformOutput',0));
trackEndPxyzVxyz(:,6)=cell2mat(arrayfun(@(i) mean(vz(i,find(~isnan(vz(i,:)),3,'last'))),[1:nTracks]','UniformOutput',0));

% get velocity components for each track from kalman filter (very similar to trackEndVxy)
xyzVel=cell2mat(arrayfun(@(iTrack) kalmanFilterInfo(trackEndTime(iTrack))...
    .stateVec(trackedFeatIndx(iTrack,trackEndTime(iTrack)),probDim+1:2*probDim),...
    [1:nTracks]','UniformOutput',0));


trackEndSpeed=sqrt(sum(xyzVel.^2,2));
vMax=prctile(trackEndSpeed,95);
vMed=median(trackEndSpeed);

% get start and end frames for each track
sFrameAll=zeros(nTracks,1);
eFrameAll=zeros(nTracks,1);
for iFrame=1:nFrames
    sFrameAll(tracksPerFrame(iFrame).starts)=iFrame;
    eFrameAll(tracksPerFrame(iFrame).ends)=iFrame;
end

% initialize matrices for pair indices and cost components
indx1 = zeros(10*nTracks,1);
indx2 = zeros(10*nTracks,1);
costComponents  = zeros(10*nTracks,5);

linkCount = 1;
for iFrame = 1:nFrames-1
    %find tracks that end in this frame
    endsToConsider = tracksPerFrame(iFrame).ends;

    if isempty(endsToConsider)
        continue
    end

    % these are the frames to consider for finding starts
    jFrame=iFrame+1:min(iFrame+tMax,nFrames);

    % find tracks that start in possible frame range
    startsToConsider=arrayfun(@(x) tracksPerFrame(x).starts,jFrame,'uniformoutput',0);
    % get number of starts in each frame
    nStarts=cellfun(@(x) length(x),startsToConsider,'uniformoutput',0);

    if isempty(startsToConsider)
        continue
    end

    % n-vector of time gaps allowable
    tGap = jFrame - iFrame;

    % forward and backward cutoff distances - based on velocity and the
    % current time gap
    cutoffDistFwd = vMax*min(sqrt(tMax),tGap);
    cutoffDistBwd = min(vMed*tMax,backVelMultFactor*vMax*tGap);

    % plot the cutoff distances as a function of time gap
    if doPlot==1 && iFrame==1
        figure;
        plot(cutoffDistFwd);
        hold on;
        plot(cutoffDistBwd,'r');
        legend('forward links','backward links','location','best')
        title('cutoff distance')
        xlabel('time gap (frames)')
        ylabel('cutoff distance (pixels)')
    end

    % make vectors containing forward/backward cutoffs for each start track candidate
    cutFwdPerVec=cell2mat(arrayfun(@(i,j) repmat(i,[j,1]),cutoffDistFwd,cell2mat(nStarts),'uniformoutput',0)');
    cutBwdPerVec=cell2mat(arrayfun(@(i,j) repmat(i,[j,1]),cutoffDistBwd,cell2mat(nStarts),'uniformoutput',0)');

    startsToConsider=cell2mat(startsToConsider');

    nStarts = length(startsToConsider);
    nEnds = length(endsToConsider);

    % end track last pt
    epX=repmat(trackEndPxyzVxyz(endsToConsider,1),[1 nStarts]);
    epY=repmat(trackEndPxyzVxyz(endsToConsider,2),[1 nStarts]);
    epZ=repmat(trackEndPxyzVxyz(endsToConsider,3),[1 nStarts]);
    % start track first pt
    spX=repmat(trackStartPxyzVxyz(startsToConsider,1)',[nEnds 1]);
    spY=repmat(trackStartPxyzVxyz(startsToConsider,2)',[nEnds 1]);
    spZ=repmat(trackStartPxyzVxyz(startsToConsider,3)',[nEnds 1]);
    % nEnds x nStarts distance matrix containing displacement vector magnitude
    dispMag=sqrt((epX-spX).^2+(epY-spY).^2+(epZ-spZ).^2);

    % we will only consider end/start pairs where the distance from end to start is less than
    % the max(forwardCutoff,backwardCutoff)
    maxCut=max(repmat(cutFwdPerVec',[nEnds 1]),repmat(cutBwdPerVec',[nEnds 1]));

    % initialize with larger-than-needed matrix
    dPerp=zeros(nStarts*nEnds,1); % distance from start track's first point to end track's lattice
    dPara=zeros(nStarts*nEnds,1); % distance from the dPerp point on end track lattice to end track's last point
    evYc=zeros(nStarts*nEnds,1);  % x component of end track's instantaneous velocity at dPerp point
    evXc=zeros(nStarts*nEnds,1);  % y component of end track's instantaneous velocity at dPerp point
    evZc=zeros(nStarts*nEnds,1);  % y component of end track's instantaneous velocity at dPerp poi
    endLinkIdx=zeros(nStarts*nEnds,1);   % end track's index for candidate pair
    startLinkIdx=zeros(nStarts*nEnds,1); % start track's index for candidate pair
    sAll=zeros(nStarts*nEnds,1);         % local start index that feeds into startCandidateIdx
    sXall=zeros(nStarts*nEnds,1);        % x-coordinates of all start tracks
    sYall=zeros(nStarts*nEnds,1);        % y-coordinates of all start tracks
    sZall=zeros(nStarts*nEnds,1);        % y-coordinates of all start tracks
    endCounter=zeros(nStarts*nEnds,1);   % keep track of the current iEnd loop index

    if doTest==1
        endRange=1:nEnds;
    else
        endRange=1:nEnds;
    end

    count=1;
    for iEnd=endRange
        % indices correspoinding to current set of starts and ends, respectively
        startCandidateIdx=find(dispMag(iEnd,:)<maxCut(iEnd,:))';
        endCandidateIdx=repmat(iEnd,[length(startCandidateIdx) 1]);

        if isempty(startCandidateIdx)
            continue
        end

        % coordinates of the first point in each startsToConsider track
        sX=trackStartPxyzVxyz(startsToConsider(startCandidateIdx),1);
        sY=trackStartPxyzVxyz(startsToConsider(startCandidateIdx),2);
        sZ=trackStartPxyzVxyz(startsToConsider(startCandidateIdx),3);

        if doTest==1
            % re-assign points to be the current end track
            xtemp=px(endsToConsider(iEnd),:)';
            xtemp(isnan(xtemp))=[];
            ytemp=py(endsToConsider(iEnd),:)';
            ytemp(isnan(ytemp))=[];
            ztemp=pz(endsToConsider(iEnd),:)';
            ztemp(isnan(ztemp))=[];

            % here only look at the end tracks >=10 frames long
            if length(xtemp)<25
                continue
            end

            % store track for later
            xTemp{iEnd,1}=xtemp;
            yTemp{iEnd,1}=ytemp;
            zTemp{iEnd,1}=ztemp;

            % add fake points around the end track's last point at each
            % pixel center +/- halfWidth pixels from this point
            halfWidth=round(max([max(cutoffDistFwd); max(cutoffDistBwd)]));
            [sX sY sZ]= meshgrid([round(xtemp(end))-halfWidth:round(xtemp(end))+halfWidth],...
                [round(ytemp(end))-halfWidth:round(ytemp(end))+halfWidth],...
                [round(ztemp(end))-halfWidth:round(ztemp(end))+halfWidth]);
            sX=sX(:);
            sY=sY(:);
            sZ=sZ(:);
            negIdx=find(sX<=0 | sY<=0 | sZ<=0 );
            sX(negIdx)=[];
            sY(negIdx)=[];
            sZ(negIdx)=[];

            sXall(count:count+length(sX)-1)=sX;
            sYall(count:count+length(sY)-1)=sY;
            sZall(count:count+length(sZ)-1)=sZ;


            % assume all fake points are within the fwd and bwd cutoff distances
            cutFwdPerVec(count:count+length(sX)-1)=max(cutoffDistFwd);
            cutBwdPerVec(count:count+length(sX)-1)=max(cutoffDistBwd);

        end
        sXall(count:count+length(sX)-1)=sX;
        sYall(count:count+length(sX)-1)=sY;
        sZall(count:count+length(sZ)-1)=sZ;
        endCounter(count:count+length(sX)-1)=iEnd;

        % call subfunction to calculate magnitude of the components of
        % the vector pointing from end to start, as well as the components
        % of the local velocity along the end track at its closest point to
        % the start track

        [dPerpTemp,dParaTemp,evYcTemp,evXcTemp,evZcTemp]=pt2segDist...
            ([pz(endsToConsider(iEnd),:)',py(endsToConsider(iEnd),:)',px(endsToConsider(iEnd),:)'],...
            trackStartPxyzVxyz(endsToConsider(iEnd),:),[sZ, sY,sX],cutoffDistBwd,0);

        % dPerp is the component perpendicular to the end track
        dPerp(count:count+length(sX)-1)=dPerpTemp;
        % dPara is the component parallel to the end track
        dPara(count:count+length(sX)-1)=dParaTemp;

        % evX/Yc are the velocity components of the end track at
        % the point closest (c) each startsToConsider track starts
        evXc(count:count+length(sX)-1)=evXcTemp;
        evYc(count:count+length(sX)-1)=evYcTemp;
        evZc(count:count+length(sX)-1)=evZcTemp;

        if doTest==1
            % redefine these based on test points
            endCandidateIdx=repmat(iEnd,[length(sX) 1]);
            endLinkIdx  (count:count+length(sX)-1) = endsToConsider(endCandidateIdx);
            startLinkIdx(count:count+length(sX)-1)=1:length(sX);
            sAll(count:count+length(sX)-1)=1:length(sX);
        else
            % starts/endsToConsider indices of the ones checked
            endLinkIdx  (count:count+length(sX)-1) =   endsToConsider(endCandidateIdx);
            startLinkIdx(count:count+length(sX)-1) = startsToConsider(startCandidateIdx);
            sAll(count:count+length(sX)-1)=startCandidateIdx;
        end

        count=count+length(sX);

    end % end iterating thru track ends

    if count==1
        % there aren't any points to continue with
        continue
    end

    % trim down the vectors
    dPerp(count:end)=[];
    dPara(count:end)=[];
    evYc(count:end)=[];
    evXc(count:end)=[];
    evZc(count:end)=[];
    endLinkIdx(count:end)=[];
    startLinkIdx(count:end)=[];
    sAll(count:end)=[];

    % velocity at starts of startsToConsider tracks
    if doTest==1
        % for test, assume they all point in same direction as end track's last pt
        svX = evXc; %trackEndPxyVxy(endLinkIdx,3);
        svY = evYc; %trackEndPxyVxy(endLinkIdx,4);
        svZ = evZc; %trackEndPxyVxy(endLinkIdx,4);

    else
        svX = trackStartPxyzVxyz(startLinkIdx,4);
        svY = trackStartPxyzVxyz(startLinkIdx,5);
        svZ = trackStartPxyzVxyz(startLinkIdx,6);
    end
    svMag = sqrt(svX.^2 + svY.^2 + svZ.^2);

    % cos of angle between start track beginning and direction of end
    % track at closest point to start
    evMagC=sqrt(evXc.^2+evYc.^2+evZc.^2);
    cosTheta = (evXc.*svX + evYc.*svY+ evZc.*svZ)./(evMagC.*svMag);

    % velocity at final point (f) of endsToConsider tracks
    evXf = trackEndPxyzVxyz(endLinkIdx,4);
    evYf = trackEndPxyzVxyz(endLinkIdx,5);
    evZf = trackEndPxyzVxyz(endLinkIdx,6);

    evMagF = sqrt(evXf.^2 + evYf.^2 +  evZf.^2);

    % displacement vector (start minus end)
    if doTest==1
        dispX = sXall-trackEndPxyzVxyz(endLinkIdx,1);
        dispY = sYall-trackEndPxyzVxyz(endLinkIdx,2);
        dispZ = sZall-trackEndPxyzVxyz(endLinkIdx,3);
    else
        dispX = trackStartPxyzVxyz(startLinkIdx,1)-trackEndPxyzVxyz(endLinkIdx,1);
        dispY = trackStartPxyzVxyz(startLinkIdx,2)-trackEndPxyzVxyz(endLinkIdx,2);
        dispZ = trackStartPxyzVxyz(startLinkIdx,3)-trackEndPxyzVxyz(endLinkIdx,3);
    end
    dispMag = sqrt(dispX.^2 + dispY.^2 + dispZ.^2);

    % cos angle between end track's end and start track's start
    cosEF_SF = (evXf.*svX + evYf.*svY+evZf.*svZ)./(evMagF.*svMag); % cos(alpha)

    % cos angle between end track's end and displacement vector
    cosEF_D  = (evXf.*dispX + evYf.*dispY+ evZf.*dispZ)./(evMagF.*dispMag); % cos(beta)

    % criteria for backward linking:
    % perp dist (dPerp) must be smaller than user-set fluctRad
    % nearest pt needs to not be the end of the endTrack and parallel dist
    % should be smaller than backward cutoff
    % angle between tracks should be less than max forward angle
    bwdIdx=find(dPerp<=(fluctRad+dPara*tan(maxBAngle)) & (dPara>0 & dPara<=cutBwdPerVec(sAll)) & cosTheta>=cos(maxFAngle));

    if ~isempty(bwdIdx)
        % record indices and parts of cost for forward links
        indx1(linkCount:linkCount+length(bwdIdx)-1) = endLinkIdx(bwdIdx);
        indx2(linkCount:linkCount+length(bwdIdx)-1) = startLinkIdx(bwdIdx);

        % cost - keep several pieces of data here for now
        % [dPerp dPara cosTheta 2 (for backward)]
        costComponents(linkCount:linkCount+length(bwdIdx)-1,1:4) = [dPerp(bwdIdx) dPara(bwdIdx) cosTheta(bwdIdx) 2*ones(length(bwdIdx),1)];
        linkCount = linkCount+length(bwdIdx);
    end

    % criteria for forward linking:
    % parallel dist (dPara) must be 0 (indicates closest pt is the end pt)
    % end-start dist must be smaller than forward cutoff
    % end-displacement angle must be smaller than max forward angle
    % angle between tracks should be less than max forward angle
    fwdIdx1=find(dPerp<=cutFwdPerVec(sAll) & dPara==0 & cosEF_D>=cos(maxFAngle) & cosTheta>=cos(maxFAngle));

    % for forward links, currently cosTheta=cosEF_SF and dPerp=dispMag
    % reassign dPerp and dPara with components of displacement vector
    if ~isempty(fwdIdx1)
        dPara(fwdIdx1)=dPerp(fwdIdx1).*cosEF_D(fwdIdx1);
        dPerp(fwdIdx1)=sqrt(dPerp(fwdIdx1).^2-dPara(fwdIdx1).^2);
    end

    % but also count those tracks falling within fluctRad of the track
    % end in the backward direction as forward links. here we calculate
    % dispMagApparent, which is the hypotenuse length of the triangle
    % formed by dPara and dPerp (as opposed to the Euclidean distance
    % between them, which is captured by dispMag)
    dispMagApparent = sqrt(dPara.^2+dPerp.^2);
    fwdIdx2=setdiff(find(dispMagApparent<=fluctRad & cosTheta>=cos(maxFAngle)),fwdIdx1);
    
    % combine them
    fwdIdx=[fwdIdx1; fwdIdx2];
    
    if ~isempty(fwdIdx)
        % record indices and parts of cost for forward links
        indx1(linkCount:linkCount+length(fwdIdx)-1) = endLinkIdx(fwdIdx);
        indx2(linkCount:linkCount+length(fwdIdx)-1) = startLinkIdx(fwdIdx);

        % cost - keep several pieces of data here for now
        % [dPerp dPara cosTheta 1 (for forward)]
        costComponents(linkCount:linkCount+length(fwdIdx)-1,1:4) = [dPerp(fwdIdx) dPara(fwdIdx) cosTheta(fwdIdx) ones(length(fwdIdx),1)];
        linkCount = linkCount+length(fwdIdx);
    end

    % some candidates will be in both fwd and bwd lists because the bwd
    % criterion specifies dPara>0, which includes some in the fluctation
    % radius captured by the fwdIdx.  we leave those in the fwdIdx list and
    % get rid of them from the bwdIdx list.
    bwdIdx=setdiff(bwdIdx,fwdIdx);

    if doTest==1

        for iEnd=endRange
            % get only those indices corresponding to the current end
            idxTemp=find(endCounter==iEnd);

            % if there are some, make an image where the forward cone is
            % green, the backward cone is red, and other pixels in blue
            if ~isempty(idxTemp)

                if isempty(intersect(idxTemp,fwdIdx)) && isempty(intersect(idxTemp,bwdIdx))
                    continue
                end
                img=zeros(max(sYall(idxTemp)),max(sXall(idxTemp)),max(sZall(idxTemp)));
                img(sub2ind(size(img),sYall(intersect(idxTemp,fwdIdx)),sXall(intersect(idxTemp,fwdIdx)),sZall(intersect(idxTemp,fwdIdx))))=1;
                img(sub2ind(size(img),sYall(intersect(idxTemp,bwdIdx)),sXall(intersect(idxTemp,bwdIdx)),sZall(intersect(idxTemp,fwdIdx))))=2;


                figure
                imagesc(img);
                hold on
                x=xTemp{iEnd,1};
                y=yTemp{iEnd,1};
                plot(x,y,'linewidth',2)
                axis equal
                scatter(x,y,'b.')
                scatter(x(end),y(end),'b*')

            end
        end

    end
end



indx1(linkCount:end) =[];
indx2(linkCount:end) =[];
costComponents(linkCount:end,:)=[];
costComponents(:,5)=sFrameAll(indx2)-eFrameAll(indx1);

% type is 1 for forward, 2 for backward
type=costComponents(:,4);

% calculate the cost
costPerpPerc=prctile(costComponents(:,1),99);
costParaPerc=prctile(costComponents(:,2),99);
if((costPerpPerc==0)||(costParaPerc==0)) costPerpPerc=1;costParaPerc=1; end;

costPerp = costComponents(:,1)./costPerpPerc; % normalized dperp only
costPara = costComponents(:,2)./costParaPerc; % normalized dpara only
costAngle = 1-costComponents(:,3);    % 1-cos(angle) only
costPerpParaAngle=costPerp+costPara+costAngle;
cost = 1.1.^costComponents(:,5).*costPerpParaAngle; % for now, this seems to be the best cost

% plot histograms of costs for forward and backward
doPlot=0;
if doPlot==1
    % define populations for forward/backward linkings costs (all tracks)
    pop1=cost(type==1); % forward
    pop2=cost(type==2); % backward

    % put them into a matrix
    M=nan(max([length(pop1) length(pop2)]),2);
    M(1:length(pop1),1)=pop1;
    M(1:length(pop2),2)=pop2;

    % create x-axis bins spanning all costs in sample
    n=linspace(min([pop1;pop2]),max([pop1;pop2]),25);

    % bin the samples
    [x1,nbins1] = histc(pop1,n); % forward
    [x2,nbins2] = histc(pop2,n); % backward

    % make new matrix of binned samples
    M=nan(max([length(x1) length(x2)]),2);
    M(1:length(x1),1)=x1;
    M(1:length(x2),2)=x2;

    % make the plot
    figure
    bar(n,M,'stack')
    colormap([1 0 0; 0 0 1])
    legend('Forward Costs','Shrinkage Costs','Location','best')
    title('Cost for all tracks')
    xlabel('cost');
    ylabel('number of track pairs');
    hold on
    deathCost=prctile(cost,90);
    plot([deathCost;deathCost],[0,max([x1+x2])])
    text(deathCost,max([x1+x2])/2,['\leftarrow ' 'death cost'])
    %saveas(gcf,[histDir filesep 'allTracksCostStackedHist.fig'])
    %saveas(gcf,[histDir filesep 'allTracksCostStackedHist.tif'])

    % get indices from endsToConsider with only 1 potential start link
    [num, u] = getMultiplicity(indx1);
    only1 = u(num==1);

    costSingles=cost(cell2mat(arrayfun(@(x) find(indx1==x),only1,'uniformoutput',0)));
    typeSingles=type(cell2mat(arrayfun(@(x) find(indx1==x),only1,'uniformoutput',0)));

    % define populations for forward/backward linkings costs (tracks with only one link)
    pop1=costSingles(typeSingles==1); % forward
    pop2=costSingles(typeSingles==2); % backward

    % put them into a matrix
    M=nan(max([length(pop1) length(pop2)]),2);
    M(1:length(pop1),1)=pop1;
    M(1:length(pop2),2)=pop2;

    % create x-axis bins spanning all costs in sample
    %n=linspace(min([pop1;pop2]),max([pop1;pop2]),25);
    % here let's use the same span for all costs

    % bin the samples
    [x1,nbins1] = histc(pop1,n); % forward
    [x2,nbins2] = histc(pop2,n); % backward

    % make new matrix of binned samples
    M=nan(max([length(x1) length(x2)]),2);
    M(1:length(x1),1)=x1;
    M(1:length(x2),2)=x2;

    % make the plot
    figure
    bar(n,M,'stack')
    colormap([1 0 0; 0 0 1])
    legend('Forward Costs','Shrinkage Costs','Location','best')
    title('Cost for tracks with only one potential link')
    xlabel('cost');
    ylabel('number of track pairs');
    hold on
    deathCost=prctile(cost,90);
    plot([deathCost;deathCost],[0,max([x1+x2])])
    text(deathCost,max([x1+x2])/2,['\leftarrow ' 'death cost'])
    %saveas(gcf,[histDir filesep 'singleTracksCostStackedHist.fig'])
    %saveas(gcf,[histDir filesep 'singleTracksCostStackedHist.tif'])

end

% plot those tracks with 3 potential connections and their costs
% doPlot=1;
if doPlot==1
    figure
    imagesc(.75*zeros(round(max(py(:))),round(max(px(:))))); colormap gray
    hold on

    [num, u] = getMultiplicity(indx1);
    c = u(num==3);
    for j=1:5 %length(c);
        b=find(indx1==c(j));
        [indx1(b) indx2(b) type(b) cost(b)]
        for i=1:length(b)
            idx=b(i);

            %get current end track's coordinates
            currentTrackE = [px(indx1(idx),:); py(indx1(idx),:)]';
            currentTrackE = currentTrackE(trackStartTime(indx1(idx)):trackEndTime(indx1(idx)),:);

            %get current start track's coordinates
            currentTrackS = [px(indx2(idx),:); py(indx2(idx),:)]';
            currentTrackS = currentTrackS(trackStartTime(indx2(idx)):trackEndTime(indx2(idx)),:);

            % plot the tracks in blue
            plot(currentTrackE(:,1),currentTrackE(:,2),'g')
            plot(currentTrackS(:,1),currentTrackS(:,2),'r')

            % plot points along tracks in red (ends) or green (starts)
            scatter(currentTrackE(:,1),currentTrackE(:,2),'b.')
            scatter(currentTrackS(:,1),currentTrackS(:,2),'b.')

            % plot possible connections
            if type(idx)==1
                x=[currentTrackE(end,1);currentTrackS(1,1)];
                y=[currentTrackE(end,2);currentTrackS(1,2)];
                plot(x,y,'c')
                text(mean(x),mean(y),[' \leftarrow ' sprintf('%3.2f',cost(idx))],'color','c');
            else
                x=[currentTrackE(end,1);currentTrackS(1,1)];
                y=[currentTrackE(end,2);currentTrackS(1,2)];
                plot(x,y,'y')
                text(mean(x),mean(y),[sprintf('%3.2f',cost(idx)) '\rightarrow '],'color','y','horizontalAlignment','right');
            end

            % end track end vectors
            %quiver(currentTrackE(end,1),currentTrackE(end,2),xyzVel(indx1(iEnd),1),xyzVel(indx1(iEnd),2),'r')
            quiver(trackEndPxyzVxyz(indx1(idx),1),trackEndPxyzVxyz(indx1(idx),2),trackEndPxyzVxyz(indx1(idx),3),trackEndPxyzVxyz(indx1(idx),4),'b')
            % start track end vectors
            %quiver(currentTrackS(end,1),currentTrackS(end,2),xyzVel(indx2(iStart),1),xyzVel(indx2(iStart),2),'r')
            quiver(trackEndPxyzVxyz(indx2(idx),1),trackEndPxyzVxyz(indx2(idx),2),trackEndPxyzVxyz(indx2(idx),3),trackEndPxyzVxyz(indx2(idx),4),'b')

        end

    end
    axis equal
end


%% Merging and splitting

%define some merging and splitting variables
numMerge  =  0; %index counting merging events
indxMerge = []; %vector storing merging track number
altCostMerge = []; %vector storing alternative costs of not merging
numSplit  =  0; %index counting splitting events
indxSplit = []; %vector storing splitting track number
altCostSplit = []; %vector storing alternative costs of not splitting

%create cost matrix without births and deaths
numEndSplit = nTracks;
numStartMerge = nTracks;
costMat = sparse(indx1,indx2,cost,numEndSplit,numStartMerge);

%% Append cost matrix to allow births and deaths ...

%determine the cost of birth and death
costBD = prctile(cost,90);

%get the cost for the lower right block
costLR = min(min(min(costMat))-1,-1);

% create cost matrix that allows for births and deaths
costMat = [costMat ... %costs for links (gap closing + merge/split)
    spdiags([costBD*ones(nTracks,1); altCostSplit],0,numEndSplit,numEndSplit); ... %costs for death
    spdiags([costBD*ones(nTracks,1); altCostMerge],0,numStartMerge,numStartMerge) ...  %costs for birth
    sparse(indx2,indx1,costLR*ones(length(indx1),1),numStartMerge,numEndSplit)]; %dummy costs to complete the cost matrix

%determine the nonlinkMarker
nonlinkMarker = min(floor(full(min(min(costMat))))-5,-5);


%% ~~~ the end ~~~


function [dPerp,dPara,evY,evX,evZ]=pt2segDist(segZYX,endTrackStartPxyVxy,ptZYX,maxCutoff,doPlot)
% find nearest point on the directed segment B-C to point P
% see http://www.geometrictools.com/Documentation/DistancePointLine.pdf for
% inspiration

% the point of this subfunction is to find the perpendicular and parallel
% distances of one or more points (stored in ptYX), which represent a
% candidate track's starting position, to a line segment (segYX), which
% represents the end track itself.  we also want to know the
% velocity components of the point on segment nearest the
% track start.  because the MT tarck may not be straight, dPara, the
% component parallel to the track, is here the actual distance along the
% lattice, not the euclidean distance from track end to the candidate start.
% dPerp is perpendicular to the line segment unless of course
% the point of interest is nearest one of the two segment ends; then it is
% defined as the Euclidean distance.  to avoid dPerp being zero in the
% backward direction (along the lattice), we add "phantom points" before
% the first point of segYX, extending the lattice backward in the direction
% of the end track's start velocity.  this extrapolation is not implemented
% for the forward direction, as dPerp = 0 is one of the criterion for a
% point to be considered in the forward direction.

if doPlot==1
    figure
end

segZYXorig=segZYX;

% find how many phantom pts are needed to extend bwd vector past max bwd cutoff
% endTrackStartPxyVxy is not the whole matrix - only for this particular end
% track
endMag=sqrt(sum(endTrackStartPxyVxy(1,4:6).^2));
if (endMag~=0)
    nRepPhantom=ceil(max(maxCutoff)/endMag);
else
    nRepPhantom=0;
end 

% create the phantom points extending back from end track's first pt,
% pointing towards it with the end track's starting velocity
phantomZYX=repmat(endTrackStartPxyVxy(1,3:-1:1),[nRepPhantom,1])-...
    repmat([nRepPhantom:-1:1]',[1,3]).*repmat(endTrackStartPxyVxy(1,6:-1:4),[nRepPhantom,1]);
segZYX=[phantomZYX; segZYX];

% treat every consecutive pair of points in segYX as a line segment BC
bZYX=segZYX(1:end-1,:); % here's a vector containing the first pt of each line segment (B)
bcZYX=diff(segZYX); % velocity components
BC=sqrt(sum(bcZYX.^2,2)); % BC length

% keep track of which entries don't exist
nanEntries=double(isnan(bcZYX(:,1)));
nanEntries=swapMaskValues(nanEntries,[0,1],[1,nan]);
nPts=size(ptZYX,1);
dPerp=zeros(nPts,1);
dPara=zeros(nPts,1);
evZ=zeros(nPts,1);
evY=zeros(nPts,1);
evX=zeros(nPts,1);
for i=1:size(ptZYX,1)
    % velocity components for vectors from all points on seg BC to P
    temp=repmat(ptZYX(i,:),[size(segZYX,1),1])-segZYX;

    bpZYX=temp(1:end-1,:); % velocity components for P-B
    BP=sqrt(sum(bpZYX.^2,2)); % distance from P to B

    cpZYX=temp(2:end,:); % velocity components for P-C
    CP=sqrt(sum(cpZYX.^2,2)); % distance from P to C

    % get fraction mag(vector from B to point on line closest to P)/mag(BC)
    % if t0<0, closest to b; if t0>0 then closet to c
    t0=(bcZYX(:,1).*bpZYX(:,1)+bcZYX(:,2).*bpZYX(:,2)+bcZYX(:,3).*bpZYX(:,3))./(BC.^2);

    D=zeros(length(t0),1);
    extraPt=zeros(length(t0),3);

    % P falls outside segment and is closest to B
    idx=find(t0<=0);
    D(idx)=BP(idx);
    extraPt(idx,:)=segZYX(idx,:); % just duplicate the first point

    % P falls outside segment and is closest to C
    idx=find(t0>=1);
    D(idx)=CP(idx);
    extraPt(idx,:)=segZYX(idx+1,:); % duplicate the last point

    % P falls within BC segment
    idx=find(t0>0 & t0<1);
    pZYX=repmat(ptZYX(i,:),[length(idx),1]);
    b_plus_t0M=bZYX(idx,:)+repmat(t0(idx),[1 3]).*bcZYX(idx,:); % location of perp point along segment
    nearVec=pZYX-b_plus_t0M;
    D(idx)=sqrt(sum(nearVec.^2,2));
    extraPt(idx,:)=b_plus_t0M; % we'll have to insert this point into the list of pts in segYX

    % don't consider where track didn't exist
    D=D.*nanEntries;

    % this is the segment index with the lowest distance
    d1Idx=find(D==nanmin(D),1,'first'); % this is where we will insert the extra point
    dPerp(i)=D(d1Idx); % distance from P to nearest point on BC

    % add in the extra point corresponding to where the dPerp vector falls on
    % the BC line
    temp=[segZYX; 0 0 0];
    temp(d1Idx+2:end,:)=temp(d1Idx+1:end-1,:);
    temp(d1Idx+1,:)=extraPt(d1Idx,:);

    % get local segment velocity at the extra point. do this by finding the
    % pt behind up and the three pts ahead of the extra pt (closest pt on
    % the line to the start pt of interest). since there is a duplicated
    % pt, this could throw off the velocity calculation if we just did a
    % mean of the differences between pts. instead, sum the differences and
    % divide by the number of *segments* (nPts-1), minus 1 for the zero
    % difference from one point to its duplicate. this should give a good
    % local estimate of the instantaneous velocity
    pts2getLocalVel = temp(max(1,d1Idx-1):min(size(temp,1),d1Idx+3),:);
    pts2getLocalVel(isnan(pts2getLocalVel(:,1)),:)=[];
    velZYX=sum(diff(pts2getLocalVel))./(size(pts2getLocalVel,1)-2);
    evZ(i)=velZYX(1);
    evY(i)=velZYX(2);
    evX(i)=velZYX(3);

    % calculate pt-to-pt displacements towards the track end and sum them
    % this is the total shrinkage distance
    % if dPara=0, P is nearest the track end (no shrinkage)
    % if dPara=trackLength, P is nearest the track start (complete shrinkage)
    dPara(i)=nansum(sqrt(sum(diff(temp(d1Idx+1:end,:)).^2,2)));

    if doPlot==1
        % show all the pts to check as blue dot
        scatter(ptZYX(:,2),ptZYX(:,1),'.')

        % plot end track as a blue line
        plot(segZYX(:,2),segZYX(:,1))
        plot(segZYXorig(:,2),segZYXorig(:,1),'LineWidth',2)
        hold on;
        % add blue dots for detection events
        scatter(temp(:,2),temp(:,1),'b.')
        % add magenta line showing start velocity vector
        quiver(endTrackStartPxyVxy(1),endTrackStartPxyVxy(2),endTrackStartPxyVxy(3),endTrackStartPxyVxy(4),0,'m')
        % plot candidate start point in red
        scatter(ptZYX(i,2),ptZYX(i,1),'r.')
        % show dPerp,dPara distances next to point
        text(ptZYX(i,2),ptZYX(i,1),['\leftarrow ' sprintf('%3.2f',dPerp(i)) ', ' sprintf('%3.2f',dPara(i))])
        % show newly-created extra pt as green circle
        scatter(extraPt(d1Idx,2),extraPt(d1Idx,1),'g')
        % show end-track calculated instantaneous velocity at new pt
        quiver(extraPt(d1Idx,2),extraPt(d1Idx,1),velZYX(2)+.1,velZYX(1)+.1,0,'r')
        % connect new pt to candidate start pt with green line
        plot([ptZYX(i,2); extraPt(d1Idx,2)],[ptZYX(i,1); extraPt(d1Idx,1)],'g')
        axis equal
    end
end