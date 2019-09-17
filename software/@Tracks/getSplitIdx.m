function idx = getSplitIdx(obj,msM)
    if(nargin < 2)
        msM = obj.getMergeSplitMatrix;
    end
    splitM = msM(msM(:,2) == 1,:);
    idx = sub2ind([obj.totalSegments obj.numTimePoints],splitM(:,3:4),[splitM(:,1) splitM(:,1)-1])';
end
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
