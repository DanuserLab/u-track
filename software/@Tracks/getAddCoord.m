function  newTracks=getAddCoord(obj,tracks)
   % tracks.addCoord(TracksToAdd)
   % add the coordinate when two tracks have the same lifetime.
   % the temporal window of /obj/ is kept if both tracks have different temporal window.
   %
   % Philippe Roudot 2017
%
% Copyright (C) 2019, Danuser Lab - UTSouthwestern 
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
   newTracks=obj.copy();
   newTracks.addCoord(tracks);
end
