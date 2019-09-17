function seqM = getSeqOfEventsMatrix(obj)
%getSeqOfEventsMatrix returns an aggregate seqOfEvents matrix combining all of the tracks
% 
% The third and forth columns are shifted to represent the segments in the segments matrix
% The fifth column represents the original compound track index number
%
% See also Tracks
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

%    seqM = cellfun( ...
%        @(s,iSeg,iTrack) [s iSeg(ones(size(s,1),1)) iTrack(ones(size(s,1),1)) ], ...
%        {obj.seqOfEvents}, ...
%        num2cell( ...
%            [0 ...
%                cumsum( ...
%                    cellfun('size',{obj(1:end-1).tracksFeatIndxCG},1) ...
%                ) ...
%            ] ...
%        ), ...
%        num2cell( 1:length(obj) ), ...
%        'UniformOutput',false);
    seqC = {obj.seqOfEvents};
    seqM = vertcat(seqC{:});
    seqI = [0 cumsum(cellfun('size',seqC(1:end-1),1))]+1;

    seqM(seqI,5) = [0 obj(1:end-1).numSegments]';
    seqM(:,5) = cumsum(seqM(:,5));

    seqM(seqI,6) = 1;
    seqM(:,6) = cumsum(seqM(:,6));
    
    seqM(:,3) = seqM(:,3) + seqM(:,5);
    seqM(:,4) = seqM(:,4) + seqM(:,5);
    seqM = seqM(:,[1:4 6]);
end
