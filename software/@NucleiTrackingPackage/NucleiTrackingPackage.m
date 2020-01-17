classdef NucleiTrackingPackage < TrackingPackage
    % A concrete process for UTrack Package
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
    
    methods (Access = public)
        function obj = NucleiTrackingPackage(varargin)
            
            % Call the superclass constructor
            obj = obj@TrackingPackage(varargin{:});
        end
        
    end
    methods (Static)

        function procConstr = getDefaultProcessConstructors(index)
            procConstr = {
                @NucleiDetectionProcess,...
                @(x,y)TrackingProcess(x,y, NucleiTrackingPackage.getDefaultTrackingParams(x,y)),...
                @MotionAnalysisProcess};
            if nargin==0, index=1:numel(procConstr); end
            procConstr=procConstr(index);
        end
        
        function funParams = getDefaultTrackingParams(owner,outputDir)
            funParams = TrackingProcess.getDefaultParams(owner,outputDir);

            % Set default gap length and allow splitting only
            funParams.gapCloseParam.timeWindow = 1;
            funParams.gapCloseParam.mergeSplit = 3;

            % Set default kalman functions (Brownian & linear motions)
            funParams.kalmanFunctions = TrackingProcess.getKalmanFunctions(1);

            % Set default cost matrices
            funParams.costMatrices(1) = TrackingProcess.getDefaultLinkingCostMatrices(owner, funParams.gapCloseParam.timeWindow,1);
            funParams.costMatrices(2) = TrackingProcess.getDefaultGapClosingCostMatrices(owner, funParams.gapCloseParam.timeWindow,1);
            
            % Set default linking parameters
            parameters = funParams.costMatrices(1).parameters;
            parameters.minSearchRadius = 5; 
            parameters.maxSearchRadius = 50; 
            funParams.costMatrices(1).parameters = parameters;
            
            % Set default gap closing parameters
            parameters = funParams.costMatrices(2).parameters;
            parameters.minSearchRadius = 5; 
            parameters.maxSearchRadius = 50; 
            parameters.ampRatioLimit = [];
            parameters.maxAngleVV = 45; 
            parameters.gapPenalty = []; 
            parameters.resLimit = 10; 
            funParams.costMatrices(2).parameters = parameters;
        end
        
        function varargout = GUI(varargin)
            % Start the package GUI
            varargout{1} = nucleiTrackingPackage(varargin{:});
        end
    end
    
end