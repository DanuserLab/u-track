function [trackType,asymParam] = getTrackType(tracks,probDim)
%GETTRACKTYPE checks whether tracks are directed or Brownian
%
%SYNOPSIS [trackType,asymParam] = getTrackType(tracks,probDim)
%
%INPUT  tracks        : -- EITHER -- 
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
%                       -- OR -- 
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
%       probDim       :     2 for 2D, 3 for 3D. Optional. Default: 2.
%
%OUTPUT trackType     :     A 1D array indicating the type of each track.
%                           1 = directed, 0 = Brownian.
%       asymParam     :     Value of asymmetry parameter used to determine
%                           track type.
%
%
%Khuloud Jaqaman, December 2007
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

trackType = [];

%% Input

%check whether correct number of input arguments was used
if nargin < 1
    disp('--getTrackType: Incorrect number of input arguments!');
    return
end

%check dimensionality
if nargin < 2 || isempty(probDim)
    probDim = 2;
end

%get number of tracks and frames
if isstruct(tracks)
    
    numTracks = length(tracks);
    seqOfEvents = vertcat(tracks.seqOfEvents);
    numFrames = max(seqOfEvents(:,1));
    
else
    
    [numTracks,numFrames] = size(tracks);
    numFrames = numFrames / 8;
    
end

%assign the asymmetry parameter thresholds that indicate directed motion
%for different track lengths
if probDim == 2
    % % %         %80th percentile:
    % % %         asymThresh = [[NaN NaN 3.5 2 1.5 1.4 1.3 1.2 1.2 1.1 1.1 1.1 1.1 1.1 1.1 ...
    % % %             1.05 1.05 1.05 1.05 1.05]'; ones(numFrames-20,1)];
    %90th percentile:
    asymThresh = [[NaN NaN 5 2.7 2.1 1.8 1.7 1.6 1.5 1.45 1.45 1.4 1.4 ...
        1.4 1.4 1.4 1.4 1.35 1.35 1.35]'; 1.3*ones(numFrames-20,1)];
    % % %     %99th percentile:
    % % %     asymThresh = [[NaN NaN 9 5 3.5 3 2.7 2.5 2.4 2.4 2.3 2.2 ...
    % % %         2.2 2.2 2.2 2.2 2.2 2.2 2.2 2.2]'; 2.1*ones(numFrames-20,1)];
else
    %90th percentile:
    asymThresh = [[NaN NaN 2.9 1.9 1.5 1.4 1.3 1.3 1.2 1.2 1.2 1.2 ...
        1.2 1.2 1.2]'; 1.1*ones(numFrames-15,1)];
    % % %     %99th percentile:
    % % %     asymThresh = [[NaN NaN 5.5 3 2.5 2.3 2 2 1.9 1.9 1.9 1.8 1.8 ...
    % % %         1.8 1.7]'; 1.7*ones(numFrames-15,1)];
end

%% asymmetry calculation

trackType = NaN(numTracks,1);
asymParam = NaN(numTracks,1);

%go over all tracks ...
for iTrack = 1 : numTracks
    
    %get coordinates of all sements in current track
    if isstruct(tracks)
        trackCoord = tracks(iTrack).tracksCoordAmpCG';
    else
        trackCoord = tracks(iTrack,:)';
    end
    
    %reshape array to get an 8-by-n array where 1st column = x, 2nd column = y, ...
    trackCoord = (reshape(trackCoord(:),8,[]))';
    trackCoord = trackCoord(:,1:probDim);
    trackCoord = trackCoord(~isnan(trackCoord(:,1)),:);

    %calculate asymmetry in track
    if size(trackCoord,1)>2
        asymParam(iTrack) = asymDeterm2D3D(trackCoord);
    else
        asymParam(iTrack) = Inf;
    end
    
    %assign track type: 1 = directed, 0 = Brownian
    trackType(iTrack) = asymParam(iTrack) > asymThresh(min(size(trackCoord,1),numFrames));
    
end

%% ~~~ the end ~~~
    

