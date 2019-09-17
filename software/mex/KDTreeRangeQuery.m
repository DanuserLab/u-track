%KDTREERANGEQUERY finds all of the points which are within the specified radius of the query points
% 
% [idx, dist] = KDTreeRangeQuery(inPts,queryPts,ranges)
% 
% This function returns the indices of the input points which are within
% the specified range of the query points. Supports 2D or 3D point sets. In
% other words, this returns all the indices of the input points which are
% contained in the cuboids whose centers are given by the query points and
% whose dimensions are given by the input range vector.
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
%     ranges - an NxK matrix specifying the ranges from each query point to
%     find input points, e.g. the dimensions of the cuboids whithin which input
%     points will be found.
%     NOTE: This value should be of class double, or strange behavior may
%     occur.
% 
% 
% Output:
% 
%   idx - Nx1 cell array, the n-th element of which gives the indices of
%   the input points which are within the n-th range the n-th query
%   point.
% 
%   dist - Nx1 cell array, the n-th element of which gives the corresponding 
%   distances between the input points and the n-th query point.
%
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
