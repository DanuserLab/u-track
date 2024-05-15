function tracksFinal = convertMat2Struct2(tracksMat)

%get number of tracks
%
% Copyright (C) 2024, Danuser Lab - UTSouthwestern 
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
numTracks = size(tracksMat,1);

%reserve memory for structure storing tracks
tracksFinal = repmat(struct('tracksFeatIndxCG',[],...
    'tracksCoordAmpCG',[],'seqOfEvents',[]),numTracks,1);

%get the start and end time of tracks
trackSEL = getTrackSEL(tracksMat);

%go over all tracks and store information
for iTrack = 1 : numTracks
    
    %track start time and end time
    startTime = trackSEL(iTrack,1);
    endTime = trackSEL(iTrack,2);
    
    if ~isnan(startTime)
        
        %feature indices
        tracksFinal(iTrack).tracksFeatIndxCG = ones(1,endTime-startTime+1);
        
        %feature coordinates and amplitudes
        tracksFinal(iTrack).tracksCoordAmpCG = full(tracksMat(iTrack,...
            (startTime-1)*8+1:endTime*8));
        
        %sequence of events
        tracksFinal(iTrack).seqOfEvents = [startTime 1 1 NaN; endTime 2 1 NaN];
        
    else
        
        tracksFinal(iTrack).tracksFeatIndxCG = 1;
        tracksFinal(iTrack).tracksCoordAmpCG = NaN(1,8);
        tracksFinal(iTrack).seqOfEvents = [1 1 1 NaN; 1 2 1 NaN];
        
    end
    
end
