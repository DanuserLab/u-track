function [trackedFeatureInfo,trackedFeatureIndx,trackStartRow,numSegments,aggregStateMat] = ...
    convStruct2MatIgnoreMS(tracksFinal,startFirstFrame)
%CONVSTRUCT2MATIGNOREMS converts tracks from structure format to matrix format, ignoring merges/splits.
%
%SYNPOSIS [trackedFeatureInfo,trackedFeatureIndx,trackStartRow,numSegments] = ...
%    convStruct2MatIgnoreMS(tracksFinal)
%
%INPUT  tracksFinal: Output of trackCloseGapsKalman, when run with
%                    gapCloseParam.mergeSplit = 1.
%       startFirstFrame: 1 to start martix at first frame where there are
%                    tracks, 0 otherwise. The point is to avoid many empty
%                    columns before any tracks start.
%                    Optional. Default: 0.
%OUTPUT trackedFeatureInfo, trackedFeatureIndx: Output of trackWithGapClosing.
%                    Every segment in tracksFinal becomes a separate track.
%       trackStartRow: Row where each compound track starts in
%                      the output matrices.
%       numSegments: Number of segments in each compound track.
%       aggregStateMat: Matrix of aggregation states, same dimensions as
%                       trackedFeatureIndx. Calculated only if aggregation
%                       state is supplied in tracksFinal.
%
%Khuloud Jaqaman, February 2008
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

%% input - output
if isfield(tracksFinal,'aggregState')
    calcAggregation = 1;
else
    calcAggregation = 0;
    aggregStateMat = [];
end

if nargin < 2 || isempty(startFirstFrame)
    startFirstFrame = 0;
end

%% conversion

%get number of tracks
numTracks = length(tracksFinal);

%get number of time points
tmp = vertcat(tracksFinal.seqOfEvents);
numTimePoints = max(tmp(:,1));

%get first frame where there are tracks and modify number of time points if
%startFirstFrame = 1
if startFirstFrame
    firstFrame = min(tmp(:,1));
    numTimePoints = numTimePoints - firstFrame + 1;
else
    firstFrame = 1;
end

%get number of segments making each track
numSegments = zeros(numTracks,1);
for iTrack = 1 : numTracks
    numSegments(iTrack) = size(tracksFinal(iTrack).tracksCoordAmpCG,1);
end

%locate the row of the first track of each compound track in the
%big matrix of all tracks (to be constructed in the next step)
trackStartRow = ones(numTracks,1);
for iTrack = 2 : numTracks
    trackStartRow(iTrack) = trackStartRow(iTrack-1) + numSegments(iTrack-1);
end

%reserve memory for matrix of tracks
trackedFeatureInfo = NaN(trackStartRow(end)+numSegments(end)-1,8*numTimePoints);

%reserve memory for matrix of aggregation state
if calcAggregation
    aggregStateMat = NaN(trackStartRow(end)+numSegments(end)-1,numTimePoints);
end

%put all tracks together in a matrix
for iTrack = 1 : numTracks
    startTime = tracksFinal(iTrack).seqOfEvents(1,1) - firstFrame + 1;
    endTime   = tracksFinal(iTrack).seqOfEvents(end,1) - firstFrame + 1;
    trackedFeatureInfo(trackStartRow(iTrack):trackStartRow(iTrack)+...
        numSegments(iTrack)-1,8*(startTime-1)+1:8*endTime) = ...
        tracksFinal(iTrack).tracksCoordAmpCG;
    if calcAggregation
        aggregStateMat(trackStartRow(iTrack):trackStartRow(iTrack)+...
            numSegments(iTrack)-1,startTime:endTime) = ...
            tracksFinal(iTrack).aggregState;
    end
end

if nargout >= 2
    
    %reserve memory for matrix of feature indices
    trackedFeatureIndx = zeros(trackStartRow(end)+numSegments(end)-1,numTimePoints);
    
    %put indices in matrix
    for iTrack = 1 : numTracks
        startTime = tracksFinal(iTrack).seqOfEvents(1,1) - firstFrame + 1;
        endTime   = tracksFinal(iTrack).seqOfEvents(end,1) - firstFrame + 1;
        trackedFeatureIndx(trackStartRow(iTrack):trackStartRow(iTrack)+...
            numSegments(iTrack)-1,startTime:endTime) = ...
            tracksFinal(iTrack).tracksFeatIndxCG;
    end
    
end

%% ~~~ the end ~~~
