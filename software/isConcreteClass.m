function status = isConcreteClass(classname)
%ISSUBCLASS check if a classname is a concrete class
%
% Synopsis  status = isConcreteClass(classname)
%
% Input:
%   classname - a string giving the name of the class
%
% Output
%   status - boolean. Return true if class is concrete (can be instantiated)
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

% Sebastien Besson, Jan 2012

% Input check
ip =inputParser;
ip.addRequired('classname',@ischar);
ip.parse(classname);


% Check validity of child class
class_meta = meta.class.fromName(classname);
assert(~isempty(class_meta),' %s is not a valid class',classname);

if isprop(class_meta, 'MethodList')
    metaMethods = class_meta.MethodList;
else
    metaMethods = [class_meta.Methods{:}];
end
status = ~any([metaMethods.Abstract]);


