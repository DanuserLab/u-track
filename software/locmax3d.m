%[lm] = locmax3d(img, wdims, varargin) searches for local maxima in 3D 
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

% Francois Aguet (last modified 12/17/2012)

function [lm] = locmax3d(img, wdims, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('img', @isnumeric);
ip.addRequired('wdims', @isnumeric);
ip.addParamValue('ClearBorder', true, @islogical);
ip.parse(img, wdims, varargin{:});

if numel(wdims)==1
    wx = wdims;
    wy = wdims;
    wz = wdims;
    if mod(wx,2)==0 || mod(wy,2)==0 || mod(wz,0)==0
        error('Mask dimensions must be odd integers');
    end
elseif numel(wdims)==3
    wx = wdims(1);
    wy = wdims(2);
    wz = wdims(3);
    if mod(wx,2)==0 || mod(wy,2)==0 || mod(wz,0)==0
        error('Mask dimensions must be odd integers');
    end    
end
mask = ones(wy,wx);

[ny,nx,nz] = size(img);
lm2D = zeros(ny,nx,nz);

for z = 1:nz
    lm2D(:,:,z) = ordfilt2(img(:,:,z), wx*wy, mask);  
end

lm = zeros(ny,nx,nz);

b = (wz-1)/2;
for z = 1:nz
    lm(:,:,z) = max(lm2D(:,:,max(1,z-b):min(nz,z+b)), [], 3);
end

% if max. filter response is equal to input, point is a local maximum
lm(lm~=img) = 0;

% set borders to zero
if ip.Results.ClearBorder
    b = (wx-1)/2;
    lm(:,[1:b end-b+1:end],:) = 0;
    b = (wy-1)/2;
    lm([1:b end-b+1:end],:,:) = 0;
    b = (wz-1)/2;
    lm(:,:,[1:b end-b+1:end]) = 0;
end