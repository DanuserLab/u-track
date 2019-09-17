function tracksCoordAmpCG = getFeatFromIdx(tracks,movieInfo)
    offsets = cumsum(cellfun('size',{movieInfo.xCoord},1));
    offsets = [ 0 offsets ];
    X = vertcat(movieInfo.xCoord)';
    Y = vertcat(movieInfo.yCoord)';
    A = vertcat(movieInfo.amp)';
    is3D = false;
    if(isfield(movieInfo,'zCoord'))
        Z = vertcat(movieInfo.zCoord)';
        is3D = true;
    else
        Z = zeros(size(X));
    end
    clear movieInfo;
    trackSEL = getTrackSEL(tracks);
    tracksCoordAmpCG = cell(1,length(tracks));
    for iTrack = 1:length(tracks)
        frames = trackSEL(iTrack,1):trackSEL(iTrack,2);
        idx = tracks(iTrack).tracksFeatIndxCG;
        nz = idx~=0;
        idx = idx + repmat(offsets(frames),size(idx,1),1);
        idxnz = idx(nz);
        m = NaN(8,size(idx,1),length(frames));
        m([1 5],nz) = X(:,idxnz);
        m([2 6],nz) = Y(:,idxnz);
        m([3 7],nz) = Z(:,idxnz);
        m([4 8],nz) = A(:,idxnz);
        m = permute(m,[2 1 3]);
        tracksCoordAmpCG{iTrack} = m(:,:);
    end
%    trackSEL = num2cell(getTrackSEL(tracks),2)';
%    [s.X,s.dX] = cellfun(@(idx,se) extractFeat(idx,se,X),{tracks.tracksFeatIndxCG},trackSEL,'UniformOutput',false);
%    [s.Y,s.dY] = cellfun(@(idx,se) extractFeat(idx,se,Y),{tracks.tracksFeatIndxCG},trackSEL,'UniformOutput',false);
%    if(is3D)
%        [s.Z,s.dZ] = cellfun(@(idx,se) extractFeat(idx,se,Z),{tracks.tracksFeatIndxCG},trackSEL,'UniformOutput',false);
%    else
%        s.Z = cellfun(@(idx) zeros(size(idx)),{tracks.tracksFeatIndxCG},'UniformOutput',false);
%        s.dZ = s.Z;
%    end
%    [s.A,s.dA] = cellfun(@(idx,se) extractFeat(idx,se,A),{tracks.tracksFeatIndxCG},trackSEL,'UniformOutput',false);
%    function [out, dout] = extractFeat(idx,se,feat)
%        frames = se(1):se(2);
%        out = NaN(size(idx));
%        nz  = idx~=0;
%        idx = idx + repmat(offsets(frames),size(idx,1),1);
%        out(nz) = feat(idx(nz),1);
%        dout(nz) = feat(idx(nz),2);
%    end
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
end
