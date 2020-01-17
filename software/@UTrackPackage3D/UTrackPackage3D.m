classdef UTrackPackage3D < TrackingPackage
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
        function obj = UTrackPackage3D(varargin)

            % Call the superclass constructor
            obj = obj@TrackingPackage(varargin{:});
        end        

        function [status, processExceptions] = sanityCheck(obj, varargin) % throws Exception Cell Array
            
            %% TODO - add more to sanitycheck
            disp('TODO: SanityCheck!');
            missingMetadataMsg = ['Missing %s! The %s is necessary to analyze '...
            '3D Tracking Movies. Please edit the movie and fill the %s.'];
            errorMsg = @(x) sprintf(missingMetadataMsg, x, x, x);
            
            assert(obj.owner_.is3D, errorMsg('MovieData is not 3D!'));
            assert(~isempty(obj.owner_.pixelSize_), errorMsg('pixel size not defined!'));
            assert(~isempty(obj.owner_.pixelSizeZ_), errorMsg('pixel Z size defined!'));
            assert(~isempty(obj.owner_.timeInterval_), errorMsg('time interval defined!'));
            [status, processExceptions] = sanityCheck@Package(obj, varargin{:});

            % possible PSF sanity check?

            % psfSigmaCheck = arrayfun(@(x) ~isempty(x.psfSigma_) || ~isempty(x.psfSigma_),obj.owner_.channels_);
            % assert(any(psfSigmaCheck),...
            %     ['Missing standard deviation of the theoretical point-spread function! '...
            %     'Please fill the numerical aperture, pixel size and'...
            %     ' emission wavelengths of all channels!']);

        end
    end

    methods (Static)

        function name = getName()
            name = 'UTrackPackage3D';
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
            varargout{1} = uTrackPackageGUI(varargin{:});
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
                @(x,y)PointSourceDetectionProcess3D(x,y,UTrackPackage3D.getDefaultDetectionParams(x,y)),...
                @(x,y)TrackingProcess(x,y,UTrackPackage3D.getDefaultTrackingParams(x,y)),...
                @MotionAnalysisProcess};
            if nargin==0, index=1:numel(procConstr); end
            procConstr=procConstr(index);
        end
        
        function funParams = getDefaultDetectionParams(owner, outputDir)

            funParams = PointSourceDetectionProcess3D.getDefaultParams(owner, outputDir);
            % Set default parameters
            funParams.OutputDirectory = [outputDir  filesep 'pointsource3D_detect'];
            funParams.alpha = .01;
            funParams = prepPerChannelParams(funParams, numel(owner.channels_));
        end

        function funParams = getDefaultTrackingParams(owner,outputDir)
            funParams = TrackingProcess.getDefaultParams(owner, outputDir);

            % Set default kalman functions
            kalmanFunctions = TrackingProcess.getKalmanFunctions(1); 
            
            fields = fieldnames(kalmanFunctions);
            validFields = {'reserveMem','initialize','calcGain','timeReverse'};
            kalmanFunctions = rmfield(kalmanFunctions,fields(~ismember(fields,validFields)));
            funParams.kalmanFunctions = kalmanFunctions;


            % --------------- gapCloseParam ----------------
            %% general gap closing parameters
            funParams.gapCloseParam.timeWindow = 2; %maximum allowed time gap (in frames) %between a track segment end and a track segment start that allows linking them.
            funParams.gapCloseParam.mergeSplit = 0; %1 if merging and splitting are to be considered, 2 if only merging is to be considered, 3 if only splitting is to be considered, 0 if no merging or splitting are to be considered.
            funParams.gapCloseParam.minTrackLen = 3; %minimum length of track segments from linking to be used in gap closing.
            %optional input:
            funParams.gapCloseParam.diagnostics = 0; %1 to plot a histogram of gap lengths in the end; 0 or empty otherwise.

            % --------------- LinkingParams ----------------
            funParams.costMatrices(1) = TrackingProcess.getDefaultLinkingCostMatrices(owner, funParams.gapCloseParam.timeWindow, 1);

            %% cost matrix for frame-to-frame linking
            %function name
            funParams.costMatrices(1).parameters.linearMotion = 1; %use linear motion Kalman filter.
            funParams.costMatrices(1).parameters.minSearchRadius = 2; 
            funParams.costMatrices(1).parameters.maxSearchRadius = 5; 
            funParams.costMatrices(1).parameters.brownStdMult = 6; 
            funParams.costMatrices(1).parameters.useLocalDensity = 0; 
            funParams.costMatrices(1).parameters.kalmanInitParam.searchRadiusFirstIteration = 4; 
            funParams.costMatrices(1).parameters.diagnostics = []; 
            
            % --------------- gapClosing Params ----------------
            funParams.costMatrices(2) = TrackingProcess.getDefaultGapClosingCostMatrices(owner, funParams.gapCloseParam.timeWindow, 1);
            funParams.costMatrices(2).parameters.linearMotion = 0; %use linear motion Kalman filter.

            funParams.costMatrices(2).parameters.minSearchRadius = 2; %minimum allowed search radius.
            funParams.costMatrices(2).parameters.maxSearchRadius = 5; %maximum allowed search radius.
            funParams.costMatrices(2).parameters.brownStdMult = 6*ones(funParams.gapCloseParam.timeWindow, 1); %multiplication factor to calculate Brownian search radius from standard deviation.
            %power for scaling the Brownian search radius with time, before and
            %after timeReachConfB (next parameter). Note that it is only the gap
            %value which is powered, then we have brownStdMult*powered_gap*sig*sqrt(dim)
            funParams.costMatrices(2).parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the gap closing and merging/splitting step.
            funParams.costMatrices(2).parameters.nnWindow = funParams.gapCloseParam.timeWindow; %number of frames before/after the current one where you want to look for a track's nearest neighbor at its end/start (in the gap closing step).
            funParams.costMatrices(2).parameters.brownScaling = [0.25 0.01];
            funParams.costMatrices(2).parameters.timeReachConfB = funParams.gapCloseParam.timeWindow; %before timeReachConfB, the search radius grows with time with the power in brownScaling(1); after timeReachConfB it grows with the power in brownScaling(2).
            funParams.costMatrices(2).parameters.ampRatioLimit = [0.7 4]; %for merging and splitting. Minimum and maximum ratios between the intensity of a feature after merging/before splitting and the sum of the intensities of the 2 features that merge/split.
            
            funParams.costMatrices(2).parameters.lenForClassify = 5; %minimum track segment length to classify it as linear or random.
            funParams.costMatrices(2).parameters.linStdMult = 1*ones(funParams.gapCloseParam.timeWindow,1); %multiplication factor to calculate linear search radius from standard deviation.
            funParams.costMatrices(2).parameters.linScaling = [0.25 0.01]; %power for scaling the linear search radius with time (similar to brownScaling).
            funParams.costMatrices(2).parameters.timeReachConfL = funParams.gapCloseParam.timeWindow; %similar to timeReachConfB, but for the linear part of the motion.
            funParams.costMatrices(2).parameters.maxAngleVV = 30; %maximum angle between the directions of motion of two tracks that allows linking them (and thus closing a gap). Think of it as the equivalent of a searchRadius but for angles.
            %optional; if not input, 1 will be used (i.e. no penalty)
            
            funParams.costMatrices(2).parameters.gapPenalty = 1.5; %penalty for increasing temporary disappearance time (disappearing for n frames gets a penalty of gapPenalty^n).
            %optional; to calculate MS search radius
            %if not input, MS search radius will be the same as gap closing search radius
            funParams.costMatrices(2).parameters.resLimit = []; %resolution limit, which is generally equal to 3 * point spread function sigma.

            % ----------------------------------------------------
            %verbose state
            funParams.verbose = 1;
            funParams.probDim = 3;
        
            %% TODO -- make export options for Amira/Ariviis/Imarus...
            funParams.saveResults.export = 0; %FLAG allow additional export of the tracking results into matrix
        end
        
    end
    
end