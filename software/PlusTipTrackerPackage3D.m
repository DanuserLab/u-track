classdef PlusTipTrackerPackage3D < TrackingPackage
    % A concrete package for tracking microtubules
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
        function obj = PlusTipTrackerPackage3D(varargin)
            
            % Call the superclass constructor
            obj = obj@TrackingPackage(varargin{:});
        end
        
        function [status, processExceptions] = sanityCheck(obj,varargin)
            
            % Check that the time interval is correctly setup
            missingMetadataMsg = ['Missing %s! The %s is necessary to analyze '...
                'microtubule plus-end tracks. Please edit the movie and fill the %s.'];
            errorMsg = @(x) sprintf(missingMetadataMsg, x, x, x);
            
            assert(~isempty(obj.owner_.pixelSize_), errorMsg('pixel size'));
            assert(~isempty(obj.owner_.pixelSizeZ_), errorMsg('pixel Z size defined!'));
            assert(~isempty(obj.owner_.timeInterval_), errorMsg('time interval'));
            assert(~isempty(obj.owner_.camBitdepth_), errorMsg('camera bit depth'));
            
            [status, processExceptions] = sanityCheck@Package(obj, varargin{:});
        end
    end
    methods (Static)
        
        function name = getName()
            name = '+TipTrackerPackage3D';
        end

        function m = getDependencyMatrix(i,j)   
            m = [0 0 0 0;  %1 ComputeMIPProcess
                 0 0 0 0;  %2 DetectionProcess
                 0 1 0 0;  %3 TrackingProcess
                 0 1 1 0;];%4 PostTrackingProcess
            if nargin<2, j=1:size(m,2); end
            if nargin<1, i=1:size(m,1); end
            m=m(i,j);
        end

        function varargout = GUI(varargin)
            % Start the package GUI
            varargout{1} = plusTipTrackerPackageGUI(varargin{:});
        end

        function classes = getProcessClassNames(index)
            classes = {
                'ComputeMIPProcess', ...
                'PointSourceDetectionProcess3D',...
                'TrackingProcess',...
                'MotionAnalysisProcess'};
            if nargin==0, index=1:numel(classes); end
            classes=classes(index);
        end
        
        function procConstr = getDefaultProcessConstructors(index)
            procConstr = {
                @ComputeMIPProcess,...
                @(x,y)PointSourceDetectionProcess3D(x,y,PlusTipTrackerPackage3D.getDefaultDetectionParams(x,y)),...
                @(x,y)TrackingProcess(x,y,PlusTipTrackerPackage3D.getDefaultTrackingParams(x,y)),...
                @MotionAnalysisProcess};
            if nargin == 0
                index = 1:numel(procConstr); 
            end
            procConstr = procConstr(index);
        end
        
        % Check this...        
        % function procConstr = getAlternateDetector()
        %     procConstr = @AnisoGaussianDetectionProcess;
        % end
        function funParams = getDefaultDetectionParams(owner, outputDir)

            funParams = PointSourceDetectionProcess3D.getDefaultParams(owner, outputDir);
            % Set default parameters
            funParams.OutputDirectory = [outputDir  filesep 'detectPlusTip_3D'];
            
            % types = PointSourceDetectionProcess3D.getDetectionTypeOptions;
            %% EB3 Detection parameters (from Philippe)
            % The detection method must be chosen depending on the data at hand
            % - 'pointSourceAutoSigmaLM' is parameter-free and suitable for raw 
            %   (but deskewed) data
            % - 'bandPassWatershed' is more suitable for deconvolve data, requires 
            %   input threshold parameter <waterThresh>.
            funParams.algorithmType = {'pointSourceAutoSigmaLM'};
            % funParams.algorithmType = {'bandPassWatershed'};

            funParams.alpha=.05;
            funParams.Mode = {'xyzAc'};
            funParams.FitMixtures = false;
            funParams.MaxMixtures = 5;
            funParams.RemoveRedundant = true;
            funParams.RedundancyRadius = .25;
            funParams.RefineMaskLoG = false;
            funParams.RefineMaskValid = false;

            % DetectComets3D & watershed params
            funParams.waterThresh = 100;
            funParams.waterStep = 10;
            funParams.lowFreq = 3;
            funParams.highFreq = 1;
            funParams.showAll = 0;
            funParams.isoCoord = true;

            funParams = prepPerChannelParams(funParams, numel(owner.channels_));
        end        

        function funParams = getDefaultTrackingParams(owner,outputDir)
            funParams = TrackingProcess.getDefaultParams(owner, outputDir);

            % Set default kalman functions
            kalmanFunctions = TrackingProcess.getKalmanFunctions(2); 
            fields = fieldnames(kalmanFunctions);
            validFields = {'reserveMem','initialize','calcGain','timeReverse'};
            kalmanFunctions = rmfield(kalmanFunctions,fields(~ismember(fields,validFields)));
            funParams.kalmanFunctions = kalmanFunctions;


            % --------------- gapCloseParam ----------------
            %% general gap closing parameters
            funParams.gapCloseParam.timeWindow = 2; %maximum allowed time gap (in frames) %between a track segment end and a track segment start that allows linking them.
            funParams.gapCloseParam.mergeSplit = 0; %1 if merging and splitting are to be considered, 2 if only merging is to be considered, 3 if only splitting is to be considered, 0 if no merging or splitting are to be considered.
            funParams.gapCloseParam.minTrackLen = 2; %minimum length of track segments from linking to be used in gap closing.
            %optional input:
            funParams.gapCloseParam.diagnostics = 0; %1 to plot a histogram of gap lengths in the end; 0 or empty otherwise.

            % --------------- LinkingParams ----------------
            funParams.costMatrices(1) = TrackingProcess.getDefaultLinkingCostMatrices(owner, funParams.gapCloseParam.timeWindow, 2);

            %% cost matrix for frame-to-frame linking
            %function name
            funParams.costMatrices(1).parameters.linearMotion = 1; %use linear motion Kalman filter.
            funParams.costMatrices(1).parameters.minSearchRadius = 2; 
            funParams.costMatrices(1).parameters.maxSearchRadius = 8; 
            funParams.costMatrices(1).parameters.brownStdMult = 3;  %searchRadiusMult
            funParams.costMatrices(1).parameters.useLocalDensity = 1; 
            funParams.costMatrices(1).parameters.kalmanInitParam.initVelocity = []; %Kalman filter initialization parameters.
            funParams.costMatrices(1).parameters.kalmanInitParam.convergePoint = []; %Kalman filter initialization parameters.
            funParams.costMatrices(1).parameters.kalmanInitParam.searchRadiusFirstIteration = 10; %Kalman filter initialization parameters.    
            funParams.costMatrices(1).parameters.diagnostics = []; 
            
            % --------------- gapClosing Params ----------------

            funParams.costMatrices(2) = TrackingProcess.getDefaultGapClosingCostMatrices(owner, funParams.gapCloseParam.timeWindow, 2);

            funParams.costMatrices(2).parameters.maxFAngle = 10; %use linear motion Kalman filter.
            funParams.costMatrices(2).parameters.maxBAngle = 10; %use linear motion Kalman filter.
            funParams.costMatrices(2).parameters.backVelMultFactor = 1.5;
            funParams.costMatrices(2).parameters.fluctRad = 1.0;
            funParams.costMatrices(2).parameters.breakNonLinearTracks = false;

            % ----------------------------------------------------
            %verbose state
            funParams.verbose = 1;
            funParams.probDim = 3;
        
            %% TODO -- make export options for Amira/Ariviis/Imarus...
            funParams.saveResults.export = 0; %FLAG allow additional export of the tracking results into matrix

        end
        
        function tools = getTools(index)
            plusTipTools(1).name = '+TIP Group Analysis';
            plusTipTools(1).funHandle = @plusTipGroupAnalysisGUI;
            if nargin==0, index=1:numel(plusTipTools); end
            tools = plusTipTools(index);
        end
    end
end