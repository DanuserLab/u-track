%STR2COLOR Convert a color character to an RGB vector
%
%    color = str2color(char)
%
% This function takes in a color character, e.g. 'k', 'r', 'g', 'b', etc.
% and returns the RGB vector for that color.
%
%IN:
%   char - One of 'k', 'b', 'c', 'r', 'm', 'y', 'w'.
%
%OUT:
%   color - 1x3 truecolor (double in range [0,1]) RGB vector for the input.
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

function x = str2color(x)
x = rem(floor((strfind('kbgcrmyw', x) - 1) * [0.25 0.5 1]), 2);
end