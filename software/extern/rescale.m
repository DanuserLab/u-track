function [B, lims] = rescale(A, lims, out_lims)
%RESCALE  Linearly rescale values in an array
%
% Examples:
%   B = rescale(A)
%   B = rescale(A, lims)
%   B = rescale(A, lims, out_lims)
%   [B lims] = rescale(A)
%
% Linearly rescales values in an array, saturating values outside limits.
%
% IN:
%   A - Input array of any size and class.
%   lims - 1x2 array of saturation limits to be used on A. Default:
%          [min(A(:)) max(A(:))].
%   out_lims - 1x2 array of output limits the values in lims are to be
%              rescaled to. Default: [0 1].
%
% OUT:
%   B - size(A) double array.
%   lims - 1x2 array of saturation limits used on A. Equal to the input
%          lims, if given.
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

% Copyright: Oliver Woodford, 2009 - 2011

% 14/01/11 Fixed a bug brought to my attention by Ming Wu (Many thanks!).

if nargin < 3
    out_lims = [0 1];
end
if nargin < 2 || isempty(lims)
    M = isfinite(A);
    if ~any(reshape(M, numel(M), 1))
        % All NaNs, Infs or -Infs
        B = double(A > 0);
        lims = [0 1];
    else
        lims = [min(A(M)) max(A(M))];
        B = normalize(A, lims, out_lims);
        B = min(max(B, out_lims(1)), out_lims(2));
    end
    clear M
else
    B = normalize(A, lims, out_lims);
    B = min(max(B, out_lims(1)), out_lims(2));
end
return

function B = normalize(A, lims, out_lims)
if lims(2) == lims(1) || out_lims(1) == out_lims(2)
    B = zeros(size(A));
else
    B = double(A);
    if lims(1)
        B = B - lims(1);
    end
    v = (out_lims(2) - out_lims(1)) / (lims(2) - lims(1));
    if v ~= 1
        B = B * v;
    end
end
if out_lims(1)
    B = B + out_lims(1);
end
return