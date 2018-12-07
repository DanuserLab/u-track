function [ movieObject, process, processID ] = getOwnerAndProcess( movieObject, processClass, createProcessIfNoneExists, varargin )
%getMovieObjectAndProcess Get the MovieObject and indicated Process
%instance
%
% INPUT
% movieObjectOrProcess:      Either a MovieObject instance or a Process instance
% processClass:              String, Process, or meta.class indicating
%                             the class of the Process we are looking for
% createProcessIfNoneExists: logical indicating whether to create the process if it
%                             does not exist (optional, default: false)
% Extra arguments will be passed to the constructor if creating a new
% process
%
% OUTPUT
% movieObject: MovieObject handle that owns the Process
% process: Process handle of type processClass, empty if it does not exist
%          and createProcess is false
% processID: index value of the Process in movieObject
%
% See also MovieObject.getOwnerAndProcess, Process.getOwnerAndProcess
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

% Mark Kittisopikul, February 2016

narginchk(2,3);

if(~ischar(processClass))
    if(isa(processClass,'Process'))
        processClass = class(processClass);
    elseif(isa(processClass,'meta.class'))
        processClass = processClass.Name;
    else
        error('MovieObject:getProcessObjectAndProcess', ... 
            'The 2nd argument, processClass, is not a string, Process, or a meta.class');
    end
end

if(nargin < 3)
    createProcessIfNoneExists = false;
end

processID = movieObject.getProcessIndex(processClass,1,false);
if(~isempty(processID))
    process = movieObject.getProcess(processID);
elseif(createProcessIfNoneExists)
    constructor = str2func(processClass);
    process = constructor(movieObject,varargin{:});
    movieObject.addProcess(process);
    processID = numel(movieObject.processes_);
else
    emptyConstructor = str2func([processClass '.empty']);
    process = emptyConstructor();
end



end

