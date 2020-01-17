classdef FigDisplay < MovieDataDisplay
    %Concreate class to display general figure plot
    % Andrew R. Jamieson Mar 2017
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
    
    properties
        plotFunc = @plot; 
        plotFunParams = {};
    end

    methods
        function obj=FigDisplay(varargin)
            obj@MovieDataDisplay(varargin{:})
        end
        
        function h = initDraw(obj, data, tag, varargin)

            ip =inputParser;
            ip.addParameter('Parent', [], @ishandle);
            ip.parse(varargin{:})
            parent_h = ip.Results.Parent;
            if isempty(obj.plotFunParams)
                h = obj.plotFunc(data);    
            elseif ~isempty(data) && ~isstruct(data) && ~isa(data.obj,'MovieData') && ~isa(data.obj,'Process')
                h = obj.plotFunc(data, obj.plotFunParams{:});
            elseif isstruct(data) && (isa(data.obj,'MovieData') || isa(data.obj,'Process'))
                h = obj.plotFunc(data, obj.plotFunParams{:},'figHandleIn',parent_h);
            end
            set(h,'Tag', tag);
        end
        function updateDraw(obj, h, data)   
        end  
    end 

    methods (Static)
        function params=getParamValidators()
            params(1).name = 'plotFunc';
            params(1).validator = @(A)validateattributes(A,{'function_handle'},{'nonempty'});
            params(2).name = 'plotFunParams';
            params(2).validator = @iscell;
        end
        function f=getDataValidator()
            f=@(x)isstruct(x) || isa(x,'MovieData') || isa(x,'Process');% (A)validateattributes(A,{'struct'},{'nonempty'});
        end
    end    
end