function trackSEL = getTrackSEL(trackedFeatureInfo,segmentSEL)
%GETTRACKSEL outputs track start times, end times and lifetimes
%
%SYNOPSIS trackSEL = getTrackSEL(trackedFeatureInfo);
%
%INPUT  trackedFeatureInfo: -- EITHER -- 
%                           Output of trackWithGapClosing:
%                           Matrix indicating the positions and amplitudes 
%                           of the tracked features to be plotted. Number 
%                           of rows = number of tracks, while number of 
%                           columns = 8*number of time points. Each row 
%                           consists of 
%                           [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                           in image coordinate system (coordinates in
%                           pixels). 
%                           If matrix is full, NaN is used to indicate time 
%                           points where the track does not exist.
%                           If matrix is sparse, NaN indicates gaps but
%                           times before and after a track are indicated by
%                           0 (they are basically unfilled elements).
%                           -- OR -- 
%                           Output of trackCloseGapsKalman:
%                           Structure array with number of entries equal to
%                           the number of tracks (or compound tracks when
%                           merging/splitting are considered). Contains the
%                           fields:
%           .tracksFeatIndxCG: Connectivity matrix of features between
%                              frames, after gap closing. Number of rows
%                              = number of track segments in compound
%                              track. Number of columns = number of frames
%                              the compound track spans. Zeros indicate
%                              frames where track segments do not exist
%                              (either because those frames are before the
%                              segment starts or after it ends, or because
%                              of losing parts of a segment.
%           .tracksCoordAmpCG: The positions and amplitudes of the tracked
%                              features, after gap closing. Number of rows
%                              = number of track segments in compound
%                              track. Number of columns = 8 * number of 
%                              frames the compound track spans. Each row
%                              consists of 
%                              [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                              NaN indicates frames where track segments do
%                              not exist.
%           .seqOfEvents     : Matrix with number of rows equal to number
%                              of events happening in a track and 4
%                              columns:
%                              1st: Frame where event happens;
%                              2nd: 1 - start of track, 2 - end of track;
%                              3rd: Index of track segment that ends or starts;
%                              4th: NaN - start is a birth and end is a death,
%                                   number - start is due to a split, end
%                                   is due to a merge, number is the index
%                                   of track segment for the merge/split.
%       segmentSEL        : Relevant only for tracks in structure format. 1
%                           to get the times for the segments inside the
%                           compound tracks, 0 to get the times for each
%                           overall compound track.
%                           Optional. Default: 0.
%
%OUTPUT trackSEL          : An array with 3 columns and number of rows equal
%                           to number of (compound) tracks. 1st column 
%                           indicates track start times, 2nd column
%                           indicates track end times and 3rd column
%                           indicates track lifetimes.
%
%REMARKS The details of compound tracks are currently ignored. The start
%times, end time and life times are for the overall compound tracks, not
%their branches.
%
%Khuloud Jaqaman, August 2006
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trackSEL = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin < 1
    disp('--getTrackSEL: Incorrect number of input arguments!');
    return
end

if nargin < 2 || isempty(segmentSEL)
    segmentSEL = 0;
end
if segmentSEL == 1 && ~isstruct(trackedFeatureInfo)
    segmentSEL = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Track information extraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get number of tracks
if isstruct(trackedFeatureInfo)
    numTracks = length(trackedFeatureInfo);
else
    numTracks = size(trackedFeatureInfo,1);
end

%alocate memory for output
if segmentSEL
    trackSEL = NaN(10*numTracks,3);
else
    trackSEL = NaN(numTracks,3);
end

%if input is in structure format
if isstruct(trackedFeatureInfo)
    
    if segmentSEL
        
        %initialize global segment index
        iSegGlob = 0;
        
        %go over all compound tracks
        for i=1:numTracks
            
            %get track's sequence of events
            seqOfEvents = trackedFeatureInfo(i).seqOfEvents;
            
            %get number of segments
            numSeg = size(seqOfEvents,1) / 2;
            
            %get the track segments' start times
            trackST = seqOfEvents(seqOfEvents(:,2)==1,1);
            
            %get rows storing segment end information
            segEndRows = seqOfEvents(seqOfEvents(:,2)==2,:);
            
            %sort the rows in ascending order of segment number
            [dummy,indxSort] = sort(segEndRows(:,3));
            segEndRows = segEndRows(indxSort,:);
            
            %extract the track segments' end times
            trackET = segEndRows(:,1);
            
            %for ends resulting from merges, subtract 1 from end time to
            %get the real track segment end time
            trackET(~isnan(segEndRows(:,4))) = ...
                trackET(~isnan(segEndRows(:,4))) - 1;
                        
            %store segment start and end times
            trackSEL(iSegGlob+1:iSegGlob+numSeg,1) = trackST;
            trackSEL(iSegGlob+1:iSegGlob+numSeg,2) = trackET;
            
            %update global segment index
            iSegGlob = iSegGlob + numSeg;
            
        end
        
        %remove all extra rows in trackSEL
        lastIndxKeep = find(~isnan(trackSEL(:,1)),1,'last');
        trackSEL = trackSEL(1:lastIndxKeep,:);
        
    else
        
        
        % vertically concatenate the r x 4 seqOfEvents
        seqOfEventsMatrix = vertcat(trackedFeatureInfo.seqOfEvents);
        
        %find track end times
        % find the indices of the last row of each compound track
        seqEndIdx = cumsum(cellfun('size',{trackedFeatureInfo.seqOfEvents},1));
        
        %find track start times
        % the first row comes after each last row except the last
        seqBegIdx = [1 seqEndIdx(1:end-1)+1];
        trackSEL(:,1:2) = [seqOfEventsMatrix(seqBegIdx,1) ...
                           seqOfEventsMatrix(seqEndIdx,1)];
        
    end
    
else %if input is a matrix
    
    %make new matrix which contains only one column per time point
    trackedFeatureInfo = trackedFeatureInfo(:,1:8:end);
    
    %if matrix is in sparse format
    if issparse(trackedFeatureInfo)
        
        %find non-empty tracks
        indxGood = find((max(abs(trackedFeatureInfo),[],2))~=0);
        
        %find track start times
        for i = indxGood'
            trackSEL(i,1) = find(trackedFeatureInfo(i,:)~=0,1,'first');
        end
        
        %find track end times
        for i = indxGood'
            trackSEL(i,2) = find(trackedFeatureInfo(i,:)~=0,1,'last');
        end
        
    else %if in full format
        
        %find non-empty tracks
        indxGood = find(~isnan(max(trackedFeatureInfo,[],2)));
        
        %find track start times
        for i = indxGood'
            trackSEL(i,1) = find(~isnan(trackedFeatureInfo(i,:)),1,'first');
        end
        
        %find track end times
        for i = indxGood'
            trackSEL(i,2) = find(~isnan(trackedFeatureInfo(i,:)),1,'last');
        end
        
    end
    
end %(if isstruct(trackedFeatureInfo))

%calculate track lifetimes
trackSEL(:,3) = trackSEL(:,2) - trackSEL(:,1) + 1;


%%%%% ~~ the end ~~ %%%%%

