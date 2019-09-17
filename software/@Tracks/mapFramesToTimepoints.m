function [ obj ] = mapFramesToTimepoints( obj, timepoints, frames )
%mapFramesToTime map frames to timepoints for the collection of Tracks
%
% INPUT
% timepoints to set the .t property to based on frames property .f
% (optional) frames to map to timepoints
%            default: min(.startFrame):max(.endFrame)
%
% OUTPUT
% the object
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

if(nargin < 3)
    frames = min([obj.startFrame]):max([obj.endFrame]);
end

assert(numel(timepoints) == numel(frames), ...
    'Tracks:mapFramesToTime:invalidNumberOfTimepoints', ...
    'The number of timepoints must match the number of frames');

framesMap(frames) = timepoints;

for ii=1:numel(obj)
    obj(ii).t = framesMap(obj(ii).f);
end


end

