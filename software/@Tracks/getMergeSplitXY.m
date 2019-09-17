function [merge, split] = getMergeSplitXY(obj, matrix, msM)
    if(nargin < 2)
        matrix = obj.getMatrix;
    end
    if(nargin < 3)
        msM = obj.getMergeSplitMatrix;
    end
    X = matrix(:,1:8:end);
    Y = matrix(:,2:8:end);
    
    merge.idx = obj.getMergeIdx(msM);
    merge.X = X(merge.idx);
    merge.Y = Y(merge.idx);
    split.idx = obj.getSplitIdx(msM);
    split.X = X(split.idx);
    split.Y = Y(split.idx);
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
