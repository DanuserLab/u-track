function ML = buildByRegexp(filter,outputDirectory)
% Build a MovieList using a regular expression.
%
% INPUT
% filter - regular expression for directories in the current directory
% which should contain .mat files of the same name as their directory
% outputDirectory - outputDirectory for the MovieList constructor
%
% OUTPUT
% ML - A MovieList containing MovieData objects that matches filter
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

% Mark Kittisopikul, March 2018
% Goldman Lab
% Northwestern University

    if(nargin < 2)
        outputDirectory = pwd;
    end
    D = dir;
    D = D([D.isdir]);
    D = D(~cellfun(@(x) isempty(regexp(x,filter, 'once')),{D.name}));
    movieDataFileNames = strcat(pwd,filesep,{D.name},filesep,{D.name},'.mat');
    ML = MovieList(movieDataFileNames,outputDirectory);
    ML.sanityCheck;
end