function obj=overlapping(obj,tr)
 % Does not work with merge and split
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
  if(length(obj)==1)
    [F,idxTr,idxObj] = intersect(tr.f,obj.f);
    % M=obj.tracksCoordAmpCG;
    % M=M(8*(min(F)-1)+(1:(8*numel(F))),:);
    % M=[nan(size(M,1),min(F)-1)  M];
    % obj=TracksHandle(M);
    obj.startFrame=min(F);
    obj.endFrame=max(F);
    obj.segmentStartFrame=min(F);
    obj.segmentEndFrame=max(F);
    obj.endFrame=max(F);
    obj.x=obj.x(idxObj);
    obj.y=obj.y(idxObj);
    obj.z=obj.z(idxObj);
  else
    arrayfun(@(o,t) o.overlapping(t),obj,tracks );
  end
end
