%CONV3FAST Fast 3D convolution with symmetric kernels
%
%  Usage:
%    [F] = conv3fast(volume, kernel)
%    [F] = conv3fast(volume, xKernel, yKernel, zKernel)
%
%  Outputs:
%     F : filtered volume
%
%  Example: convolution with a Gaussian kernel
%     s = 2;
%     w = ceil(4*s);
%     g = exp(-(0:w).^2/(2*s^2)); % symmetric kernel starts at '0'
%     F = conv3fast(data, g);
%
%  Note: NaNs in input are allowed
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

% Francois Aguet, 09/19/2013