function [ movieInfo ] = getMovieInfo( obj )
%getMovieInfo Reconstruct movieInfo struct from Tracks
%
% See also setMovieInfo, setFeatFromIdx, getFeatFromIdx
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

maxFrame = max([obj.endFrame]);
maxIdx = max(cellfun(@(x) max(x(:)),{obj.tracksFeatIndxCG}));
movieInfo(maxFrame,1) = struct('xCoord',[],'yCoord',[],'zCoord',[],'amp',[]);

x = NaN(2,maxIdx,maxFrame);
y = x;
z = x;
A = x;

[linIdx,coords] = arrayfun(@getIndexFrameLinearCoordinates,obj,'UniformOutput',false);

coords = [coords{:}];

x(:,[linIdx{:}]) = [coords.x];
y(:,[linIdx{:}]) = [coords.y];
z(:,[linIdx{:}]) = [coords.z];
A(:,[linIdx{:}]) = [coords.A];

x = permute(x,[2 1 3]);
y = permute(y,[2 1 3]);
z = permute(z,[2 1 3]);
A = permute(A,[2 1 3]);

x = num2cell(x,[1 2]);
y = num2cell(y,[1 2]);
z = num2cell(z,[1 2]);
A = num2cell(A,[1 2]);

[movieInfo.xCoord] = x{:};
[movieInfo.yCoord] = y{:};
[movieInfo.zCoord] = z{:};
[movieInfo.amp] = A{:};

function [linIdx,coords] =  getIndexFrameLinearCoordinates(track)
    nzIdx = track.tracksFeatIndxCG ~= 0;
    nzIdx = nzIdx(:);
    idx = track.tracksFeatIndxCG(nzIdx);
    frame = repmat(track.f,track.numSegments,1);
    frame = frame(nzIdx);
    linIdx = sub2ind([maxIdx maxFrame],idx,frame);
    linIdx = linIdx(:)';
    
    coords.x = [ joinColumns(track.x(nzIdx))  joinColumns(track.dx(nzIdx)) ]';
    coords.y = [ joinColumns(track.y(nzIdx))  joinColumns(track.dy(nzIdx)) ]';
    coords.z = [ joinColumns(track.z(nzIdx))  joinColumns(track.dz(nzIdx)) ]';
    coords.A = [ joinColumns(track.A(nzIdx))  joinColumns(track.dA(nzIdx)) ]';
end

end
