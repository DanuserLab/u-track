%KDTREEBALLQUERY finds all of the points which are within the specified radius of the query points
% 
% [idx, dist] = KDTreeBallQuery(inPts,queryPts,radii)
% 
% This function returns the indices of the input points which are within
% the specified radii of the query points. Supports 1D, 2D or 3D point sets. 
% In other words, this returns all the indices of the input points which are
% contained in the spheres whose centers are given by the query points and
% whose radii are given by the input radii vector.
%
% Input:
% 
%     inPts - an MxK matrix specifying the input points to test for distance
%     from the query points, where M is the number of points and K is the
%     dimensionality of the points.
% 
%     queryPts - an NxK matrix specifying the query points, e.g. the centers of
%     the spheres within which input points will be found.
% 
%     radii - an Nx1 vector or a scalar specifying the distances from each 
%     query point to find input points, e.g. the radii of the spheres within 
%     which input points will be found. If scalar, the same radius is used
%     for all query points.
%     NOTE: This value should be of class double, or strange behavior may
%     occur.
% 
% 
% Output:
% 
%   idx - Nx1 cell array, the n-th element of which gives the indices of
%   the input points which are within the n-th radii of the n-th query
%   point.
% 
%   dist - Nx1 cell array, the n-th element of which gives the corresponding 
%   distances between the input points and the n-th query point.
%
%
% Copyright (C) 2020, Danuser Lab - UTSouthwestern 
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
