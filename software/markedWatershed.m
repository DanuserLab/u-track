function [label]=markedWatershed(vol,scales,thresh,varargin)
ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('vol', @isnumeric);
ip.addRequired('scales', @isnumeric);
ip.addRequired('thresh', @isnumeric);
ip.addParamValue('WindowSize', []); % Default: 2*sigma, see fitGaussians3D.
ip.parse(vol, scales,thresh, varargin{:});

if ~isa(vol, 'double')
    vol = double(vol);
end

if numel(scales)==1
    scales = [scales scales];
end

ws = ip.Results.WindowSize;
if isempty(ws)
    ws = ceil(2*scales);
elseif numel(ws)==1
    ws = [ws ws];
end

%-------------------------------------------------------------------------------------------
% Convolutions
%-------------------------------------------------------------------------------------------
% right-hand side of symmetric kernels
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
gx = exp(-(0:ws(1)).^2/(2*scales(1)^2));
gz = exp(-(0:ws(2)).^2/(2*scales(2)^2));
fg = conv3fast(vol, gx, gx, gz);
fu =  conv3fast(vol,    ones(1,ws(1)+1), ones(1,ws(1)+1), ones(1,ws(2)+1));
fu2 = conv3fast(vol.^2, ones(1,ws(1)+1), ones(1,ws(1)+1), ones(1,ws(2)+1));

% Laplacian of Gaussian-filtered input
gx2 = (0:ws(1)).^2 .*gx;
gz2 = (0:ws(2)).^2 .*gz;
fgx2 = conv3fast(vol, gx2, gx, gz);
fgy2 = conv3fast(vol, gx, gx2, gz);
fgz2 = conv3fast(vol, gx, gx, gz2);
imgLoG = (2/scales(1)^2+1/scales(2)^2)*fg - ((fgx2+fgy2)/scales(1)^4 + fgz2/scales(2)^4);
clear fgx2 fgy2 fgz2;

imseriesshow(imgLoG)

label=watershed(-imgLoG);
label(smooth3(vol,'gaussian',[3 3 3],scales(1))<ip.Results.thresh)=0;
pstruct=[];
