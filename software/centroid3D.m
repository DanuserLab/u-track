function ce = centroid3D(img,exponent,matrixCoords)
% CENTROID compute the centroid of a gray value patch
%
% SYNOPSIS ce = centroid3D(img, exponent)
%
% INPUT img : an image 3D patch matrix
%       exponent: (opt) if the image is large and noisy, increase the exponent to
%                 get better results!
%       matrixCoords : (opt) if 1, centroid is returned in matrix
%                 coordinates. Default: 0
% 
% OUTPUT ce : vector with the centroid coordinates (in image coords!)
%
% REMARKS : image patches masked with NaNs will not be counted
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

% c 19/04/00

if nargin < 2 || isempty(exponent)
    exponent = 1;
end
if nargin < 3 || isempty(matrixCoords)
    matrixCoords = false;
end

[s1,s2,s3] = size(img);
% reshape image only once
%img = permute(img,[2,1,3]);
img = img(:).^exponent;

% use ndgrid so that we get xyz in matrix coordinates
[x,y,z] = ndgrid(1:s1,1:s2,1:s3);
% nansum in case there are masked regions
cx = nansum(x(:).*img);
cy = nansum(y(:).*img);
cz = nansum(z(:).*img);

% swap cy, cx so that we get image coordinates
if matrixCoords
    ce = [cx, cy, cz]/nansum(img);
else
    ce = [cy, cx, cz]/nansum(img);
end

% old code (replaced 10/05 by jonas)
%
% for l = 1 : s(2)
%    cx = cx + sum(sum(img(:,l,:).*l));
% end;
% for l = 1 : s(1)
%    cy = cy + sum(sum(img(l,:,:).*l));
% end;
% for l = 1 : s(3)
%    cz = cz + sum(sum(img(:,:,l).*l));
% end;
% sTot=nansum(img(:));
% ce=[cx cy cz]/sTot;