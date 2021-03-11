classdef HandleObject < handle
    %Class to create a handle object of any other object in order to
    %allow pass-by-reference. Create a reference object using syntax:
    %<name> = HandleObject(<thing-that-is-referred-to>). Access the 
    %underlying object by calling the <name>.Object property. 
%
% Copyright (C) 2021, Danuser Lab - UTSouthwestern 
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
    
    %Honestly, this should be a standard Matlab method, like with adding
    %the @ symbol to a function name to pass it as a handle.
    
% Updated in Jan 2020 to incorporate the changes made by Carmen Klein Herenbrink 
% and Brian Devree from Copenhagen University to reduce the tracking time.
% This is a new function provided by them.

    properties
      Object=[];
   end
 
   methods
      function obj=HandleObject(receivedObject)
         obj.Object=receivedObject;
      end
   end
end
