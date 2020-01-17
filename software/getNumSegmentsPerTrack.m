function numSegments = getNumSegmentsPerTrack(tracksFinal)
%GETNUMSEGMENTS gives back the number of segments in each compound track
%
%SYNPOSIS numSegments = getNumSegmentsPerTrack(tracksFinal)
%
%INPUT  tracksFinal: Output of trackCloseGapsKalman
%OUTPUT numSegments: Vector of number of segments in each compound track
%
%Khuloud Jaqaman, 8 September 2009
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

%find number of tracks
numTracks = length(tracksFinal);

%initialize numSegments
numSegments = NaN(numTracks,1);

%go over the tracks and get number of segments
for iTrack = 1 : numTracks
    numSegments(iTrack) = size(tracksFinal(iTrack).tracksCoordAmpCG,1);
end

