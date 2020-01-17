function save( filename, varargin )
%cached.save saves a MAT file as per the builtin save, but also invalidates
%the cache for that file for cached.load
%
% See also save
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

% Inspired by Sebastien Besson

% Mark Kittisopikul
% December 2014

% execute the builtin save function in the caller workspace
expr = strjoin([filename,varargin],''',''');
expr = [ 'builtin(''save'',''' expr ''')'];
evalin('caller',expr);

% clear the cache for filename
cached.load(filename,'-clear');

end

