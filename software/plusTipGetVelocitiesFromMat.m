function [trackedFeatureInfo,trackedFeatureInfoInterp,trackInfo,trackVelocities,timeRange] = plusTipGetVelocitiesFromMat(trackedFeatureInfo,movieInfo,minTrackLen,timeRange)
% plusTipGetVelocitiesFromMat fills forward and backward gaps in EB3 trajectory and calculates segment velocities
%
% INPUT: trackedFeatureInfo : nTracks x 8*nFrames matrix resulting from
%                             tracking or simulated data from simEBtracks.m
%        movieInfo          : Output of feature tracking
%        minTrackLen        : minimum number of frames a track segment
%                             should be to get considered. segments shorter
%                             than minTrackLen will be connected with the
%                             adjacent gaps to make a bigger gap
%        timeRange          : frame range over which to run post-processing
%
% OUTPUT: trackedFeatureInfoInterp : same size as trackedFeatureInfo but
%                                    has gaps filled in with interpolated
%                                    positions
%         trackInfo                : nTrack x 1 structure with the
%                                    following fields:
%                                    .seg (segment)
%                                    .fgap (forward gap)
%                                    .bgap (backward gap)
%                                    .ugap (unclassifiable gap)
%                                    each field is an nSeg(nGap) x 4 matrix
%                                    containing the track number, start
%                                    frame, end frame, and average velocity
%         trackVelocities          : structure with fields .frame2frame and
%                                    .segmentAvgs, where each is an
%                                    nTracks x nFrames-1 matrix. entry ij
%                                    contains the velocity of the ith
%                                    microtubule between the j and j+1
%                                    frames. the former gives instantaneous
%                                    velocities whereas the latter gives
%                                    averages over the whole segment or gap
%         timeRange                : frame range over which post-processing
%                                    ran - needed for post-processing fn
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


% if tracksFinal is the input, convert to matrix format
if isstruct(trackedFeatureInfo)
    [trackedFeatureInfo,trackedFeatureIndx] = convStruct2MatNoMS(trackedFeatureInfo,movieInfo);
    clear trackedFeatureIndx
end
nFrames = size(trackedFeatureInfo,2)/8;


if nargin<4 || isempty(timeRange)
    startFrame=1;
    endFrame=nFrames;
elseif isequal(unique(size(timeRange)),[1 2])
    if timeRange(1)<=timeRange(2) && timeRange(2)<=nFrames
        startFrame = timeRange(1);
        endFrame = timeRange(2);
    else
        startFrame = 1;
        endFrame = nFrames;
    end

else
    error('--plusTipGetVelocitiesFromMat: timeRange should be [startFrame endFrame] or [] for all frames')
end

temp=nan(size(trackedFeatureInfo));
temp(:,8*(startFrame-1)+1:8*endFrame)=trackedFeatureInfo(:,8*(startFrame-1)+1:8*endFrame);
trackedFeatureInfo=temp;

% find whole rows full of nans and remove them
nanRowIdx=find(sum(isnan(trackedFeatureInfo),2)==size(trackedFeatureInfo,2));
if ~isempty(nanRowIdx)
    trackedFeatureInfo(nanRowIdx,:)=[];
end

% since we need a few frames to calculate segment directions for the dot
% product calculation (to determine whether forward (growing) or backward
% (shrinking) velocity, tracks should be at least 3 frames long
minTrackLen = max(3,minTrackLen);

% number of tracks and frames
nTracks = size(trackedFeatureInfo,1);
nFrames = size(trackedFeatureInfo,2)/8;

% initialize matrices to hold velocity and position info
perFrameVelocities=nan*zeros(nTracks,nFrames-1);
trackedFeatureInfoInterp = nan.*(trackedFeatureInfo);

% initialize structure to keep track of segments, forward gaps, backward
% gaps, and unclassifiable gaps
trackInfo(nTracks,1).seg  = [];
trackInfo(nTracks,1).fgap = [];
trackInfo(nTracks,1).bgap = [];
trackInfo(nTracks,1).ugap = [];


% find and fill gaps using linear interpolation of known positions
for iTrack=1:nTracks
    % get coordinates of the detected feature for this track
    px=trackedFeatureInfo(iTrack,1:8:end);
    py=trackedFeatureInfo(iTrack,2:8:end);

    % find start/end frames for each segment of the track
    trackIdx = find(~isnan(px));
    temp = diff(~isnan(px));
    trackStarts = unique([trackIdx(1) find(temp == 1)+1]);
    trackEnds   = unique([trackIdx(end) find(temp == -1)]);
    trackLengths = trackEnds - trackStarts + 1;

    % remove segments that are less than minTrackLen
    tooShort = find(trackLengths < minTrackLen);
    for iShortSeg=1:length(tooShort)
        px(trackStarts(tooShort(iShortSeg)):trackEnds(tooShort(iShortSeg)))=nan;
        py(trackStarts(tooShort(iShortSeg)):trackEnds(tooShort(iShortSeg)))=nan;
    end

    % again find start/end frames for each segment of the track
    trackIdx = find(~isnan(px));
    temp=(diff(~isnan(px)));
    if ~isempty(trackIdx) % after removal, we may have an empty row
        trackStarts = unique([trackIdx(1) find(temp == 1)+1]);
        trackEnds = unique([trackIdx(end) find(temp == -1)]);
        
        % fill in segment info [trackIndex trackStart trackEnd 0], where the 0
        % will be filled with the segment velocity
        trackInfo(iTrack,1).seg = [iTrack*ones(size(trackStarts')) trackStarts' trackEnds' zeros(size(trackStarts'))];
    else
        trackStarts=-1;
    end



    % pxFill/pyFill are the coordinates of the track where we have
    % estimated the position of the feature (using a constant velocity) to fill the gap.
    pxFill=px; pyFill=py;

    % fill in gaps and get gap info
    % if n segs then n-1 gaps
    if length(trackStarts) > 1 % then there's a gap
        gapStarts = trackEnds(1:end-1);
        gapEnds = trackStarts(2:end);
        gapLengths = gapEnds-gapStarts+1;

        for iGap=1:length(gapStarts)

            % get xy-vectors showing direction/speed-per-frame of the parts of the tracks
            % directly flanking the gap
            beforeGapDirXY =[(px(gapStarts(iGap))-px(gapStarts(iGap)-2))/2; (py(gapStarts(iGap))-py(gapStarts(iGap)-2))/2];
            afterGapDirXY = [(px(gapEnds(iGap)+2)-px(gapEnds(iGap)))/2; (py(gapEnds(iGap)+2)-py(gapEnds(iGap)))/2];

            % get xy-vector showing direction/speed-per-frame of the particle during the gap
            gapDirXY = [px(gapEnds(iGap))-px(gapStarts(iGap)); py(gapEnds(iGap))-py(gapStarts(iGap))];

            % take dot product of before/after and gap vectors and determine if
            % a forward or backward gap
            dotProd = [dot(beforeGapDirXY,gapDirXY); dot(afterGapDirXY,gapDirXY)];

            % classify as a forward or backward gap (or unclassified, if dot
            % poduct with both segments doesn't have the same sign
            if all(dotProd>0) % forward gap
                trackInfo(iTrack,1).fgap = [trackInfo(iTrack,1).fgap; iTrack gapStarts(iGap) gapEnds(iGap) 0];

            elseif all(dotProd<0) % backward gap
                trackInfo(iTrack,1).bgap = [trackInfo(iTrack,1).bgap; iTrack gapStarts(iGap) gapEnds(iGap) 0];

            else % one of dot products is pos, one neg...can't classify
                trackInfo(iTrack,1).ugap = [trackInfo(iTrack,1).ugap; iTrack gapStarts(iGap) gapEnds(iGap) 0];
            end

            pxFill(gapStarts(iGap):gapEnds(iGap))=linspace(px(gapStarts(iGap)),px(gapEnds(iGap)),gapLengths(iGap));
            pyFill(gapStarts(iGap):gapEnds(iGap))=linspace(py(gapStarts(iGap)),py(gapEnds(iGap)),gapLengths(iGap));

        end % end filling gaps, getting info
    end % end if there's a gap

    % find average speed for the track from detected positions (those in
    % input matrix) - this corresponds to average speed in forward
    % direction. we will give negative sign to backward gaps in the next
    % step.
    perFrameVelocities(iTrack,:)=sqrt(diff(pxFill,[],2).^2+diff(pyFill,[],2).^2);

    for iBackGap=1:size(trackInfo(iTrack,1).bgap,1)
        % get gap start and end frame
        gapS=trackInfo(iTrack,1).bgap(iBackGap,2);
        gapE=trackInfo(iTrack,1).bgap(iBackGap,3);

        % make velocities negative for backward gaps
        perFrameVelocities(iTrack,gapS:gapE-1)= -perFrameVelocities(iTrack,gapS:gapE-1);

        % record average per-frame velocity over backward gap
        trackInfo(iTrack,1).bgap(iBackGap,4) = mean(perFrameVelocities(iTrack,gapS:gapE-1));
    end

    for iForwGap=1:size(trackInfo(iTrack,1).fgap,1)
        % get gap start and end frame
        gapS=trackInfo(iTrack,1).fgap(iForwGap,2);
        gapE=trackInfo(iTrack,1).fgap(iForwGap,3);

        % record average per-frame velocity over forward gap
        trackInfo(iTrack,1).fgap(iForwGap,4) = mean(perFrameVelocities(iTrack,gapS:gapE-1));
    end

    for iUnclassGap=1:size(trackInfo(iTrack,1).ugap,1)
        % get gap start and end frame
        gapS=trackInfo(iTrack,1).ugap(iUnclassGap,2);
        gapE=trackInfo(iTrack,1).ugap(iUnclassGap,3);

        % record average per-frame velocity over backward gap
        trackInfo(iTrack,1).ugap(iUnclassGap,4) = mean(perFrameVelocities(iTrack,gapS:gapE-1));
    end

    for iSeg=1:size(trackInfo(iTrack,1).seg,1)
        % get segment start and end frame
        segS=trackInfo(iTrack,1).seg(iSeg,2);
        segE=trackInfo(iTrack,1).seg(iSeg,3);

        % record average per-frame velocity over backward gap
        trackInfo(iTrack,1).seg(iSeg,4) = mean(perFrameVelocities(iTrack,segS:segE-1));
    end

    % now fill in the output matrix
    for iFrame = 1:nFrames
        if ~isnan(pxFill(iFrame))
            trackedFeatureInfoInterp(iTrack,8*(iFrame-1)+1:8*(iFrame-1)+8)=[pxFill(iFrame) pyFill(iFrame) zeros(1,6)];
        end
    end

end
% get rows that don't contain a track
nanRowIdx=find(sum(isnan(trackedFeatureInfoInterp),2)==size(trackedFeatureInfoInterp,2));
% remove these rows from relevant data
if ~isempty(nanRowIdx)
    trackedFeatureInfo(nanRowIdx,:)=[];
    trackedFeatureInfoInterp(nanRowIdx,:)=[];
    perFrameVelocities(nanRowIdx,:)=[];
end
nTracks = size(trackedFeatureInfo,1);


% removing NaN-data from trackInfo is a bit trickier...
fNames=fieldnames(trackInfo);
for iNan=1:length(nanRowIdx)
    for iTrack=nanRowIdx(iNan):length(trackInfo)
        for iField=1:length(fNames)
            if ~isempty(trackInfo(iTrack).(fNames{iField}))
                trackInfo(iTrack).(fNames{iField})(:,1)= trackInfo(iTrack).(fNames{iField})(:,1)-1;
            end
        end
    end
end
trackInfo(nanRowIdx)=[];

% concatenate the segment and gap data
binnedData = [vertcat(trackInfo.seg); vertcat(trackInfo.fgap); vertcat(trackInfo.bgap); vertcat(trackInfo.ugap)];
avgSegVel = zeros(nTracks,nFrames-1);
for i=1:size(binnedData,1)
    avgSegVel(binnedData(i,1),binnedData(i,2):binnedData(i,3)-1) = binnedData(i,4);
end

trackVelocities.frame2frame = perFrameVelocities;
trackVelocities.segmentAvgs = avgSegVel;
