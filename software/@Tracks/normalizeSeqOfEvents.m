function [ obj ] = normalizeSeqOfEvents( obj )
%normalizeSeqOfEvents If seqOfEvents is empty, then create a default one.
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
    if(isa(obj,'TracksHandle'))
        % TracksHandle.seqOfEvents will always be populated
        % see TracksHandle.get.seqOfEvents
        return;
    end
    noSeqOfEvents = cellfun('isempty',{obj.seqOfEvents});
    tracksToFix = obj(noSeqOfEvents);
    if(isempty(tracksToFix))
        return;
    end
    seqOfEvents = cellfun( ...
        @(nF) [1 1 1 NaN; nF 2 1 NaN], ...
        {tracksToFix.numFrames}, ...
        'UniformOutput',false);
    [tracksToFix.seqOfEvents] = deal(seqOfEvents{:});
%     obj(noSeqOfEvents) = tracksToFix;
end

