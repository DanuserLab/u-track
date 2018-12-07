function patch=stamp3d(data,patchSize,center,symmetric)
%STAMP3D copy 3D sub-patch out of larger 3D data set
%
% SYNOPSIS patch=stamp3d(data,patchSize,center)
%
% INPUT data   : 3D data
%       patchSize  : size of patch
%       center : position of center pixel of patch in 3D data
%       symmetric : (opt) if part of the patch would fall outside of the 
%                    img: whether to cut accordingly on the other side 
%                    [{0}/1]
%
% OUTPUT patch : 3D patch 
%
% Copyright (C) 2018, Danuser Lab - UTSouthwestern 
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

% c: 18/6/01	dT

% test input, assign defaults
if nargin < 4 || isempty(symmetric)
    symmetric = 0;
end

% expand patchSize, center if necessary (as with data of lower
% dimensionality)
tmp = ones(1,3);
tmp(1:length(patchSize)) = patchSize;
patchSize = tmp;
tmp = ones(1,3);
tmp(1:length(center)) = center;
center = tmp;

% find patch in data
[ds1,ds2,ds3]=size(data);
hl=floor(patchSize/2);

% find extension towards 0
hx1=min([center(1)-1,hl(1)]);
hy1=min([center(2)-1,hl(2)]);
hz1=min([center(3)-1,hl(3)]);

% find extension towards inf
hx2 = min([hl(1),ds1-center(1)]);
hy2 = min([hl(2),ds2-center(2)]);
hz2 = min([hl(3),ds3-center(3)]);

% make symmetric, if necessary
if symmetric
   [hx1,hx2] = deal(min(hx1,hx2));
   [hy1,hy2] = deal(min(hy1,hy2));
   [hz1,hz2] = deal(min(hz1,hz2));
end

% take patch aout of img.
patch=data(center(1)-hx1:center(1)+hx2,...
    center(2)-hy1:center(2)+hy2,...
    center(3)-hz1:center(3)+hz2);

