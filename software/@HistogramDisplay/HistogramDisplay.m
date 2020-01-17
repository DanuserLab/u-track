classdef HistogramDisplay < MovieDataDisplay
    %Concreate display class for displaying points or lines
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
        Marker = 'none';
        Linewidth = 2;
        XLabel = '';
        YLabel = '' ;
        sfont = {'FontName', 'Helvetica', 'FontSize', 18};
        lfont = {'FontName', 'Helvetica', 'FontSize', 22};
        tfont = {'FontName', 'Helvetica', 'FontSize', 14, 'FontAngle', 'italic'};
    end
    methods
        function obj = HistogramDisplay(varargin)
            nVarargin = numel(varargin);
            if nVarargin > 1 && mod(nVarargin,2)==0
                for i=1 : 2 : nVarargin-1
                    obj.(varargin{i}) = varargin{i+1};
                end
            end
        end
        function h=initDraw(obj,data,tag,varargin)
            
            % Generate plot
            hold on;
            if ~isempty(data)
                [n,x]=hist(data, 1:max(data));
                h=bar(x, n, 'Linewidth', obj.Linewidth, 'Tag', tag);
            else
                h = text(.3,.4,'No data found', obj.lfont{:});
            end
            xlabel(obj.XLabel,obj.lfont{:});
            ylabel(obj.YLabel,obj.lfont{:});
            set(gca, 'LineWidth', 1.5, obj.sfont{:});
        end
        function updateDraw(obj,h,data)
            tag = get(h(1),'Tag');
            cla(get(get(h(1),'Parent')))
            obj.initDraw(data,tag);
        end
    end    
    
    methods (Static)
        function params=getParamValidators()
            params(1).name='Marker';
            params(1).validator=@ischar;
            params(2).name='Linewidth';
            params(2).validator=@isscalar;
            params(3).name='XLabel';
            params(3).validator=@ischar;
            params(4).name='YLabel';
            params(4).validator=@ischar;
            params(5).name='sfont';
            params(5).validator=@iscell;
            params(6).name='lfont';
            params(6).validator=@iscell;
            params(7).name='tfont';
            params(7).validator=@iscell;
        end
        
        function f=getDataValidator()
            f=@isnumeric;
        end
    end    
end