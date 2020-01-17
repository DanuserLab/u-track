function [trackIndx,errFlag] = chooseTracks(trackedFeatureInfo,criteria,probDim)
%CHOOSETRACKS outputs the indices of tracks that satisfy the input criteria
%
%SYNOPSIS [trackIndx,errFlag] = chooseTracks(trackedFeatureInfo,criteria);
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
%                           pixels). NaN is used to indicate time points 
%                           where the track does not exist.
%                           -- OR -- 
%                           Output of trackCloseGapsKalman:
%                           Structure array with number of entries equal to
%                           the number of tracks (or compound tracks when
%                           merging/splitting are considered). Contains the
%                           fields:
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
%       criteria          : Structure with fields:
%           .lifeTime        :Structure with fields:
%               .min            :minimum lifetime.
%               .max            :maximum lifetime.
%           .startTime       :Structure with fields:
%               .min            :minimum start time.
%               .max            :maximum start time.
%           .endTime         :Structure with fields:
%               .min            :minimum end time.
%               .max            :maximum end time.
%           .initialAmp      :Structure with fields:
%               .min            :minimum initial amplitude.
%               .max            :maximum initial amplitude.
%           .finalAmp        :Structure with fields:
%               .min            :minimum initial amplitude.
%               .max            :maximum initial amplitude.
%           .initialXCoord   :Structure with fields:
%               .min            :minimum initial x-coordinate.
%               .max            :maximum initial x-coordinate.
%           .initialYCoord   :Structure with fields:
%               .min            :minimum initial y-coordinate.
%               .max            :maximum initial y-coordinate.
%           .initialZCoord   :Structure with fields:
%               .min            :minimum initial z-coordinate.
%               .max            :maximum initial z-coordinate.
%           .numSegments     : Structure with fields:
%               .min            :minimum number of segments making a
%                                compound track.
%               .max            :maximum number of segments making a
%                                compound track.
%           .trackType       :0/1 to choose random/linear tracks.
%                           All criteria are optional. Leave out or give as
%                           [] if not of interest.
%       probDim           : 2 for 2D, 3 for 3D. Optional. Default: 2.
%
%OUTPUT trackIndx         : Indices of (compound) tracks that satisfy the input criteria.
%       errFlag           : 0 if function executes normally, 1 otherwise
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

%% Output

trackIndx = [];
errFlag = 0;

%% Input

%check whether correct number of input arguments was used
if nargin < 2
    disp('--chooseTracks: Incorrect number of input arguments!');
    errFlag = 1;
    return
end

if nargin < 3 || isempty(probDim)
    probDim = 2;
end

%% Preamble

%get number of tracks
if isstruct(trackedFeatureInfo)
    numTracks = length(trackedFeatureInfo);
else
    numTracks = size(trackedFeatureInfo,1);
end

%get track start, end and life times
trackSEL = getTrackSEL(trackedFeatureInfo);

%initialize vector of track indices
trackIndx = ones(numTracks,1);

%% Lifetime
if isfield(criteria,'lifeTime') && ~isempty(criteria.lifeTime)
    
    %find track lifetimes
    comparisonVec = trackSEL(:,3);

    %get minimum lifetime
    if isfield(criteria.lifeTime,'min') && ~isempty(criteria.lifeTime.min)
        minCrit = criteria.lifeTime.min;
    else
        minCrit = min(comparisonVec) - 1;
    end
    
    %get maximum lifetime
    if isfield(criteria.lifeTime,'max') && ~isempty(criteria.lifeTime.max)
        maxCrit = criteria.lifeTime.max;
    else
        maxCrit = max(comparisonVec) + 1;
    end
    
    %assign one to tracks that satisfy this criterion plus all criteria
    %above
    trackIndx = trackIndx & comparisonVec >= minCrit & comparisonVec <= maxCrit;
end

%% Start time
if isfield(criteria,'startTime') && ~isempty(criteria.startTime)
    
    %find track start times
    comparisonVec = trackSEL(:,1);

    %get minimum start time
    if isfield(criteria.startTime,'min') && ~isempty(criteria.startTime.min)
        minCrit = criteria.startTime.min;
    else
        minCrit = min(comparisonVec) - 1;
    end
    
    %get maximum start time
    if isfield(criteria.startTime,'max') && ~isempty(criteria.startTime.max)
        maxCrit = criteria.startTime.max;
    else
        maxCrit = max(comparisonVec) + 1;
    end
    
    %assign one to tracks that satisfy this criterion plus all criteria
    %above
    trackIndx = trackIndx & comparisonVec >= minCrit & comparisonVec <= maxCrit;
end

%% End time
if isfield(criteria,'endTime') && ~isempty(criteria.endTime)
    
    %find track end times
    comparisonVec = trackSEL(:,2);

    %get minimum end time
    if isfield(criteria.endTime,'min') && ~isempty(criteria.endTime.min)
        minCrit = criteria.endTime.min;
    else
        minCrit = min(comparisonVec) - 1;
    end
    
    %get maximum end time
    if isfield(criteria.endTime,'max') && ~isempty(criteria.endTime.max)
        maxCrit = criteria.endTime.max;
    else
        maxCrit = max(comparisonVec) + 1;
    end
    
    %assign one to tracks that satisfy this criterion plus all criteria
    %above
    trackIndx = trackIndx & comparisonVec >= minCrit & comparisonVec <= maxCrit;
end

%% Initial amplitude
if isfield(criteria,'initialAmp') && ~isempty(criteria.initialAmp)
    
    %find initial amplitudes
    comparisonVec = zeros(numTracks,1);
    if isstruct(trackedFeatureInfo)
        for i=1:numTracks
            comparisonVec(i) = trackedFeatureInfo(i).tracksCoordAmpCG(1,4);
        end
    else
        for i=1:numTracks
            comparisonVec(i) = trackedFeatureInfo(i,(trackSEL(i,1)-1)*8+4);
        end
    end
    
    %get minimum initial amplitude
    if isfield(criteria.initialAmp,'min') && ~isempty(criteria.initialAmp.min)
        minCrit = criteria.initialAmp.min;
    else
        minCrit = min(comparisonVec) - 1;
    end
    
    %get maximum initial amplitude
    if isfield(criteria.initialAmp,'max') && ~isempty(criteria.initialAmp.max)
        maxCrit = criteria.initialAmp.max;
    else
        maxCrit = max(comparisonVec) + 1;
    end
    
    %assign one to tracks that satisfy this criterion plus all criteria
    %above
    trackIndx = trackIndx & comparisonVec >= minCrit & comparisonVec <= maxCrit;
end

%% Final amplitude
if isfield(criteria,'finalAmp') && ~isempty(criteria.finalAmp)
    
    %find final amplitudes
    comparisonVec = zeros(numTracks,1);
    if isstruct(trackedFeatureInfo)
        for i=1:numTracks
            comparisonVec(i) = max(trackedFeatureInfo(i).tracksCoordAmpCG(:,end-4));
        end
    else
        for i=1:numTracks
            comparisonVec(i) = trackedFeatureInfo(i,(trackSEL(i,2)-1)*8+4);
        end
    end
    
    %get minimum final amplitude
    if isfield(criteria.finalAmp,'min') && ~isempty(criteria.finalAmp.min)
        minCrit = criteria.finalAmp.min;
    else
        minCrit = min(comparisonVec) - 1;
    end
    
    %get maximum final amplitude
    if isfield(criteria.finalAmp,'max') && ~isempty(criteria.finalAmp.max)
        maxCrit = criteria.finalAmp.max;
    else
        maxCrit = max(comparisonVec) + 1;
    end
    
    %assign one to tracks that satisfy this criterion plus all criteria
    %above
    trackIndx = trackIndx & comparisonVec >= minCrit & comparisonVec <= maxCrit;
end

%% Initial x-coordinate
if isfield(criteria,'initialXCoord') && ~isempty(criteria.initialXCoord)
    
    %find initial x-coordinates
    comparisonVec = zeros(numTracks,1);
    if isstruct(trackedFeatureInfo)
        for i=1:numTracks
            comparisonVec(i) = trackedFeatureInfo(i).tracksCoordAmpCG(1,1);
        end
    else
        for i=1:numTracks
            comparisonVec(i) = trackedFeatureInfo(i,(trackSEL(i,1)-1)*8+1);
        end
    end

    %get minimum initial x-coordinate
    if isfield(criteria.initialXCoord,'min') && ~isempty(criteria.initialXCoord.min)
        minCrit = criteria.initialXCoord.min;
    else
        minCrit = min(comparisonVec) - 1;
    end
    
    %get maximum initial x-coordinate
    if isfield(criteria.initialXCoord,'max') && ~isempty(criteria.initialXCoord.max)
        maxCrit = criteria.initialXCoord.max;
    else
        maxCrit = max(comparisonVec) + 1;
    end
    
    %assign one to tracks that satisfy this criterion plus all criteria
    %above
    trackIndx = trackIndx & comparisonVec >= minCrit & comparisonVec <= maxCrit;
end

%% Initial y-coordinate
if isfield(criteria,'initialYCoord') && ~isempty(criteria.initialYCoord)
    
    %find initial y-coordinates
    comparisonVec = zeros(numTracks,1);
    if isstruct(trackedFeatureInfo)
        for i=1:numTracks
            comparisonVec(i) = trackedFeatureInfo(i).tracksCoordAmpCG(1,2);
        end
    else
        for i=1:numTracks
            comparisonVec(i) = trackedFeatureInfo(i,(trackSEL(i,1)-1)*8+2);
        end
    end
    
    %get minimum initial y-coordinate
    if isfield(criteria.initialYCoord,'min') && ~isempty(criteria.initialYCoord.min)
        minCrit = criteria.initialYCoord.min;
    else
        minCrit = min(comparisonVec) - 1;
    end
    
    %get maximum initial y-coordinate
    if isfield(criteria.initialYCoord,'max') && ~isempty(criteria.initialYCoord.max)
        maxCrit = criteria.initialYCoord.max;
    else
        maxCrit = max(comparisonVec) + 1;
    end
    
    %assign one to tracks that satisfy this criterion plus all criteria
    %above
    trackIndx = trackIndx & comparisonVec >= minCrit & comparisonVec <= maxCrit;
end

%% Initial z-coordinate
if isfield(criteria,'initialZCoord') && ~isempty(criteria.initialZCoord)
    
    %find initial z-coordinates
    comparisonVec = zeros(numTracks,1);
    if isstruct(trackedFeatureInfo)
        for i=1:numTracks
            comparisonVec(i) = trackedFeatureInfo(i).tracksCoordAmpCG(1,3);
        end
    else
        for i=1:numTracks
            comparisonVec(i) = trackedFeatureInfo(i,(trackSEL(i,1)-1)*8+3);
        end
    end
    
    %get minimum initial z-coordinate
    if isfield(criteria.initialZCoord,'min') && ~isempty(criteria.initialZCoord.min)
        minCrit = criteria.initialZCoord.min;
    else
        minCrit = min(comparisonVec) - 1;
    end
    
    %get maximum initial z-coordinate
    if isfield(criteria.initialZCoord,'max') && ~isempty(criteria.initialZCoord.max)
        maxCrit = criteria.initialZCoord.max;
    else
        maxCrit = max(comparisonVec) + 1;
    end
    
    %assign one to tracks that satisfy this criterion plus all criteria
    %above
    trackIndx = trackIndx & comparisonVec >= minCrit & comparisonVec <= maxCrit;
end

%% Number of segments
if isfield(criteria,'numSegments') && ~isempty(criteria.numSegments)
    
    %calculate number of segments in each track
    if isstruct(trackedFeatureInfo)
        comparisonVec = getNumSegmentsPerTrack(trackedFeatureInfo);
    else
        comparisonVec = ones(size(trackedFeatureInfo,1),1);
    end

    %get minimum number of segments
    if isfield(criteria.numSegments,'min') && ~isempty(criteria.numSegments.min)
        minCrit = criteria.numSegments.min;
    else
        minCrit = min(comparisonVec) - 1;
    end
    
    %get maximum number of segments
    if isfield(criteria.numSegments,'max') && ~isempty(criteria.numSegments.max)
        maxCrit = criteria.numSegments.max;
    else
        maxCrit = max(comparisonVec) + 1;
    end
    
    %assign one to tracks that satisfy this criterion plus all criteria
    %above
    trackIndx = trackIndx & comparisonVec >= minCrit & comparisonVec <= maxCrit;
end

%% Track type
if isfield(criteria,'trackType') && ~isempty(criteria.trackType)
    
    %get track types
    comparisonVec = getTrackType(trackedFeatureInfo,probDim);

    %assign one to tracks that satisfy this criterion plus all criteria
    %above
    trackIndx = trackIndx & comparisonVec == criteria.trackType;
    
end

%% Output

%keep only the indices of tracks that satisfy all input criteria
trackIndx = find(trackIndx);


%% ~~~ the end ~~~

