function [ oldIdx ] = reindex( obj, idx )
%reindex Reindex array elements by their current array position or as
%specified
%
% INPUT
% obj - Tracks object array
% idx - (optional) new indices
%
% OUTPUT
% oldIdx - Previous indices or empty if not available
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

if(nargout > 0)
    oldIdx = [obj.index];
end

if(nargin < 2)
    idx = 1:numel(obj);
end

idx = num2cell(idx);
[obj.index] = idx{:};


end

