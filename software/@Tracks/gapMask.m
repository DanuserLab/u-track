function [ varargout ] = gapMask( obj )
%gapMask creates a binary mask of where gaps are in the compound track
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
    if(isscalar(obj))
        mask = obj.aliveMask;
        mask = mask & obj.tracksFeatIndxCG == 0;
        varargout{1} = mask;
    else
        varargout = cellfun(@gapMask,num2cell(obj(1:nargout)),'UniformOutput',false);
    end

end

