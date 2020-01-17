function [gapInfo,gapInfoSpace] = findTrackGaps(trackedFeatureInfo)
%FINDTRACKGAPS finds the gaps in each track and gives back their location and length
%
%SYNOPSIS [gapInfo,gapInfoSpace] = findTrackGaps(trackedFeatureInfo)
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
%
%OUTPUT gapInfo           : An array with 6 columns:
%                           1st column: (compound) track to which gap
%                           belongs;
%                           2nd column: segment within track to which gap
%                           belongs;
%                           3rd column: frame where gap starts; 
%                           4th column: gap length.
%                           5th column: ratio of gap length to length of
%                           segment before it.
%                           6th column: ratio of gap length to length of
%                           segment after it.
%       gapInfoSpace      : An array with 1 column, continuing gapInfo,
%                           storing net displacement during each gap.
%
%Khuloud Jaqaman, February 2007
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

gapInfo = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%check whether correct number of input arguments was used
if nargin < 1
    disp('--trackedFeatureInfo: Incorrect number of input arguments!');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pre-processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get number of tracks and number of time points
if isstruct(trackedFeatureInfo) %if tracks are in structure format
    numTracks = length(trackedFeatureInfo);
    tmp = vertcat(trackedFeatureInfo.seqOfEvents);
    if isempty(tmp)
        return
    end
    numTimePoints = max(tmp(:,1));
    clear tmp
else %if tracks are in matrix format
    [numTracks,numTimePoints] = size(trackedFeatureInfo);
    numTimePoints = numTimePoints/8;
    if numTimePoints == 0
        return
    end
end

%put tracks into big matrix if the input is in structure format
if isstruct(trackedFeatureInfo)

    %store the input structure as a variable with a different name
    inputStructure = trackedFeatureInfo;
    clear trackedFeatureInfo;
    
    %get number of segments making each track
    numSegments = zeros(numTracks,1);
    for i = 1 : numTracks
        numSegments(i) = size(inputStructure(i).tracksCoordAmpCG,1);
    end

    %locate the row of the first track of each compound track in the
    %big matrix of all tracks (to be constructed in the next step)
    trackStartRow = ones(numTracks,1);
    for iTrack = 2 : numTracks
        trackStartRow(iTrack) = trackStartRow(iTrack-1) + numSegments(iTrack-1);
    end

    %put all tracks together in a matrix
    trackedFeatureInfo = NaN*ones(trackStartRow(end)+numSegments(end)-1,8*numTimePoints);
    for i = 1 : numTracks
        startTime = inputStructure(i).seqOfEvents(1,1);
        endTime = inputStructure(i).seqOfEvents(end,1);
        trackedFeatureInfo(trackStartRow(i):trackStartRow(i)+...
            numSegments(i)-1,8*(startTime-1)+1:8*endTime) = ...
            inputStructure(i).tracksCoordAmpCG;
    end

    %clear memory
    clear inputStructure

else %if input was already in matrix format

    numSegments = ones(numTracks,1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Track information extraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get track start times and end times
trackSEL = getTrackSEL(trackedFeatureInfo);

%make new matrix which contains only one column per time point
trackedFeatureInfoOrig = trackedFeatureInfo;
trackedFeatureInfo = trackedFeatureInfo(:,1:8:end);

%alocate memory for output (assume that each track will have 10 gaps on average)
gapInfo = zeros(10*numTracks,6);
gapInfoSpace = zeros(10*numTracks,1);

%assign value of index showing last stored position in gapInfo.
indxLast = 0;

%look for gaps in each track in the big matrix
iBig = 0;
for i = 1 : numTracks    
    for j = 1 : numSegments(i)
        
        %update index enumerating tracks in the big matrix
        iBig = iBig + 1;
        
        %get current track, its starting point and its ending point
        track0 = trackedFeatureInfo(iBig,:);
        start0 = trackSEL(iBig,1);
        end0 = trackSEL(iBig,2);

        %find all time points where track is NaN between its beginning and its end
        missing = find(isnan(track0))';
        missing = missing(missing > start0 & missing < end0);

        %if there are gaps in track ...
        if ~isempty(missing)

            %find the time increment between one NaN and the next one
            missingDiff = diff(missing);

            %find all places where time increment is larger than 1 - this is the
            %start of a gap
            gapStart = missing([1;find(missingDiff~=1)+1]);

            %find places just before those with an increment larger than 1 - this
            %is the end of a gap
            gapEnd = missing([find(missingDiff~=1);end]);

            %calculate gap length
            gapLength = gapEnd - gapStart + 1;

            %get number of gaps in track
            numGaps = length(gapStart);
            
            %get lengths of segments before and after gaps
            segmentLength = [gapStart-1; end0] - [start0; gapEnd+1] + 1;
            
            %calculate ratio of gap length to segments before it and after
            %it
            lengthRatio = [gapLength ./ segmentLength(1:end-1) ...
                gapLength ./ segmentLength(2:end)];

            %place the gaps of current track in gapInfo
            gapInfo(indxLast+1:indxLast+numGaps,:) = [i*ones(numGaps,1) ...
                j*ones(numGaps,1) gapStart gapLength lengthRatio];
            
            %output displacement during each gap
            xCoord = trackedFeatureInfoOrig(iBig,1:8:end);
            yCoord = trackedFeatureInfoOrig(iBig,2:8:end);
            dispGap = sqrt((xCoord(gapEnd+1)-xCoord(gapStart-1)).^2 + ...
                (yCoord(gapEnd+1)-yCoord(gapStart-1)).^2);
            gapInfoSpace(indxLast+1:indxLast+numGaps,:) = dispGap;
            
            %update indxLast
            indxLast = indxLast + numGaps;

        end %(if ~isempty(missing))

    end %(for j = 1 : numSegment(i))
end %(for i = 1 : numTracks)

%remove any unused entries in gapInfo
indxKeep = find(gapInfo(:,1)~=0);
gapInfo = gapInfo(indxKeep,:);
gapInfoSpace = gapInfoSpace(indxKeep,:);


%%%%% ~~ the end ~~ %%%%%

