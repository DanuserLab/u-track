function D=createDistanceMatrix(M,N)
% createDistanceMatrix calculates the distance matrix for two sets of points
%
% SYNOPSIS   D=createDistanceMatrix(M,N)
%
% INPUT      M and N are the matrices containing the set of point coordinates.
%            M and N can represent point positions in 1, 2 and 3D, as follows.
%            
%            In 1D: M=[ x1        and   N=[ x1
%                       x2                  x2
%                       ...                ... 
%                       xm ]                xn ]
%
%            In 2D:
%                   M=[ y1 x1     and   N=[ y1 x1
%                       y2 x2              y2 x2
%                        ...                ...
%                       ym xm ]            yn xn ]
%
%            In 3D:
%                   M=[ y1 x1 z1  and   N=[ y1 x1 z1
%                       y2 x2 z2            y2 x2 z2
%                         ...                ...
%                       ym xm zm ]          yn xn zn ]
%
%
% OUTPUT   D : distance matrix D=(dij), i=1..m, j=1..n
% 
% REMARK   For 1D, both positive and negative distances are returned.
%
% C-MEX file - Aaron Ponti 28/08/15
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
