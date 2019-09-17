function obj=addCoord(obj,tracks)
   % tracks.addCoord(TracksToAdd)
   % addition coordinate, lifetime must interset
   % TODO: 
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
   
   if((length(tracks)~=1)&&(length(tracks)~=length(obj)))
      error('Added tracks set must have the same size or unitary.')
   end

   if(length(obj)==1)
     tr=tracks.getOverlapping(obj);
     if(obj.lifetime~=tr.lifetime)
             error('tracks lifetime must overlap entirely with obj')
      end
          
     obj.x=obj.x+tr.x;
     obj.y=obj.y+tr.y;
     obj.z=obj.z+tr.z;
   else
    if(length(tracks)==1)
      arrayfun(@(o) o.addCoord(tracks),obj,'unif',0);
    else
      arrayfun(@(o,t) o.addCoord(t),obj,tracks,'unif',0);
   end
 end
