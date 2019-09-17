function tracks=fillTrackGaps(tracks)
    if(~isempty(tracks))
    se=[zeros(1,tracks.numTimePoints()) 1 ones(1,tracks.numTimePoints())];
    for tIdx=1:length(tracks)
        gi=isnan(tracks(tIdx).x);
        if(any(gi))
            copyIdx=1:tracks(tIdx).lifetime;
            copyIdx(gi)=0;
            copyIdx=imdilate(copyIdx,se);
            tracks(tIdx).x=tracks(tIdx).x(copyIdx);
            tracks(tIdx).y=tracks(tIdx).y(copyIdx);
            tracks(tIdx).z=tracks(tIdx).z(copyIdx);
        end
    end
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
