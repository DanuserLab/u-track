function  [vert,edges,frames,edgesLabel]=getGraph(obj)
   % Philippe Roudot 2018
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
   if(numel(obj)==1)
         obj=fillTrackGaps(obj);
   		vert=[obj.x' obj.y' obj.z'];
   		edges=[1:numel(obj.f)-1 ; 2:numel(obj.f)]';
   		frames=obj.f(2:end)';
         edgesLabel=ones(size(edges,1),1);
   else
   		[vert,edges,frames,edgesLabel]=arrayfun(@(t) t.getGraph,obj,'unif',0);
         NID=0;
   		for tIdx=1:numel(edges)
            edges{tIdx}=edges{tIdx}+NID;
            NID=NID+size(vert{tIdx},1);
            edgesLabel{tIdx}=tIdx*edgesLabel{tIdx};
   		end

   		vert=vertcat(vert{:});
   		edges=vertcat(edges{:});
   		frames=vertcat(frames{:});
         edgesLabel=vertcat(edgesLabel{:});
   end

end
