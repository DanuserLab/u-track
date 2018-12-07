function ce = centroid2D(img)
% CENTROID compute the centoid of a gray value patch
%
% SYNOPSIS ce = centroid2D(img)
%
% INPUT img : an image patch matrix
% 
% OUTPUT ce : vector with the centroid coordinates
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

s = size(img);
cx = 0;
cy = 0;
for l = 1 : s(2);
   cx = cx + sum(img(:,l).*l);
end
for l = 1 : s(1);
   cy = cy + sum(img(l,:).*l);
end
sTot=sum(img(:));
ce=[cx cy]/sTot;