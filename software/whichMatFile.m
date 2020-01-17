function [ fullFileName ] = whichMatFile( filename )
%whichMatFile Finds the location of the filename on the path
%
% The function first searches in the current working directory (quick). If
% this fails, then it uses matfile to locate the file.
% 
% INPUT
% filename is an absolute or relative path for matfile
%
% OUTPUT
% fullFileName is the full file name of the file that would be loaded by
% the builtin load function
%
% See also matfile
%
% Mark Kittisopikul
% December 2014
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

% check if file exists in pwd or is absolute 
[s,attrib] = fileattrib(filename);
if(s)
    fullFileName = attrib.Name;
else
    % this is a potentially expensive operation
    matObj = matfile(filename);
    fullFileName = [];
    if(exist(matObj.Properties.Source,'file'))
        fullFileName = matObj.Properties.Source;
    end
    % free memory
    delete(matObj);
end

end

