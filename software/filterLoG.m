% filterLoG : filters a signal with a Laplacian of Gaussian filter. Unlike
% the built-in Matlab function fspecial('Laplacian'), this function computes
% the exact convolution, using FFTs.
%
%    out = filterLoG(signal, sigma);
%
%    INPUT: signal: input (1-3 dimensional)
%           sigma : standard deviation of the Gaussian kernel
%
%    OUTPUT: y : LoG filtered signal
%
% Francois Aguet, last modified: 11/11/2010
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

function y = filterLoG(signal, sigma)
dims = ndims(signal);
if numel(signal)==max(size(signal))
    dims = 1;
end

switch dims
    case 1
        nx = length(signal);
        w1 = -nx/2:nx/2-1;
        w1 = fftshift(w1);
        w1 = w1*2*pi/nx;
        
        I = fft(signal);
        LoG = w1.^2 .* exp(-0.5*sigma^2.*w1.^2);
        y = real(ifft(I.*LoG));
    case 2
        [ny,nx] = size(signal);
        
        [w1,w2] = meshgrid(-nx/2:nx/2-1, -ny/2:ny/2-1);
        w1 = fftshift(w1);
        w2 = fftshift(w2);
        
        w1 = w1*2*pi/nx;
        w2 = w2*2*pi/ny;
        I = fft2(signal);
        LoG = (w1.^2 + w2.^2) .* exp(-0.5*sigma^2*(w1.^2 + w2.^2));
        y = real(ifft2(I.*LoG));
    case 3
        [ny,nx,nz] = size(signal);
        
        [w1,w2,w3] = meshgrid(-nx/2:nx/2-1, -ny/2:ny/2-1, -nz/2:nz/2-1);
        w1 = fftshift(w1);
        w2 = fftshift(w2);
        w3 = fftshift(w3);

        w1 = w1*2*pi/nx;
        w2 = w2*2*pi/ny;
        w3 = w3*2*pi/nz;
        I = fftn(signal);
        LoG = (w1.^2 + w2.^2 + w3.^2) .* exp(-0.5*sigma^2*(w1.^2 + w2.^2 + w3.^2));
        y = real(ifftn(I.*LoG));
end