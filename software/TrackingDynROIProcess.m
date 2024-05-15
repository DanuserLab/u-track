classdef TrackingDynROIProcess < TrackingProcess
    % Process Class for tracking for dynamic ROI
    % trackMovie.m is the wrapper function
    % TrackingDynROIProcess is part of New Utrack 3D package
    %
    % This process class is a subclass of TrackingProcess class.
    % So different getName and funParams.OutputDirectory can be put for this process.
    %
    % Qiongjing (Jenny) Zou, July 2019
%
% Copyright (C) 2024, Danuser Lab - UTSouthwestern 
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
    
    methods(Access = public)
        
        function obj = TrackingDynROIProcess(owner, varargin)
            if nargin == 0
                super_args = {};
            else
                % Input check
                ip = inputParser;
                ip.addRequired('owner',@(x) isa(x,'MovieData'));
                ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
                ip.addOptional('funParams',[],@isstruct);
                ip.parse(owner,varargin{:});
                outputDir = ip.Results.outputDir;
                funParams = ip.Results.funParams;
                
                % Define arguments for superclass constructor
                super_args{1} = owner;
                super_args{2} = TrackingDynROIProcess.getName;
                super_args{3} = @trackMovie;
                if isempty(funParams)
                    funParams = TrackingDynROIProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            obj = obj@TrackingProcess(super_args{:});
            obj.is3Dcompatible_ = true;
        end
                
        function varargout = loadChannelOutput(obj, iChan, varargin)
            % This functions was adapted from TrackingProcess.loadChannelOutput
            % The changes are (1) tracks results for Dyn ROI are saved in tracksFinalDynROIRef_oldFormat instead of in tracksFinal for TrackingProcess.
            % (2) the struct of tracksFinalDynROIRef_oldFormat are in a transposed format as the struct of tracksFinal, 
            % so need to transpose back as s.tracksFinalDynROIRef_oldFormat'
            
            % Input check
            outputList = {'tracksFinal', 'gapInfo', 'staticTracks',...
                          'plottracks3d','tracksFinalDynROIRef_oldFormat'};
            ip =inputParser;
            ip.addRequired('obj');
            ip.addRequired('iChan', @(x) obj.checkChanNum(x));
            ip.addOptional('iFrame', [] ,@(x) obj.checkFrameNum(x));
            ip.addParameter('useCache', obj.useCache_, @islogical);
            ip.addParameter('iZ',[], @(x)ismember(x,1:obj.owner_.zSize_));
            ip.addParameter('output', outputList{1}, @(x) all(ismember(x,outputList)));
            ip.addParameter('projectionAxis3D','Z', @(x) ismember(x,{'Z','X','Y','three'}));
            ip.parse(obj,iChan,varargin{:})
            output = ip.Results.output;
            iFrame = ip.Results.iFrame;
            projAxis3D = ip.Results.projectionAxis3D;
            iZ = ip.Results.iZ;
            if ischar(output),output={output}; end
            
            % Data loading
            if ~strcmp(projAxis3D,'Z') %&& obj.owner_.is3D % && ~obj.funParams_.isoOutput
                s = cached.load(obj.funParams_.processBuildDynROI.outFilePaths_{3}, '-useCache', true);
                s1 = cached.load(s.movieDataDynROICell{1}, '-useCache', true); % s1.MD is a movieData built based on DynROI raw images.
                ZXRatio = s1.MD.pixelSizeZ_/s1.MD.pixelSize_;
            end

            s = cached.load(obj.outFilePaths_{4,iChan}, '-useCache', ip.Results.useCache, 'tracksFinalDynROIRef_oldFormat');

            varargout = cell(numel(output), 1);
        
            
            for i = 1:numel(output)
                switch output{i}
                    case {'tracksFinalDynROIRef_oldFormat', 'staticTracks', 'plottracks3d'}
                        varargout{i} = s.tracksFinalDynROIRef_oldFormat';
                    case 'gapInfo'
                        varargout{i} = findTrackGaps(s.tracksFinalDynROIRef_oldFormat');
                end
                if strcmp(output{i}, 'tracksFinalDynROIRef_oldFormat') && ~isempty(iFrame)
                    % Filter tracks existing in input frame
                    trackSEL=getTrackSEL(s.tracksFinalDynROIRef_oldFormat');
                    validTracks = (iFrame>=trackSEL(:,1) &iFrame<=trackSEL(:,2));
                    [varargout{i}(~validTracks).tracksCoordAmpCG]=deal([]);
                    
                    nFrames = iFrame-trackSEL(validTracks,1)+1;
                    nCoords = nFrames*8;
                    validOut = varargout{i}(validTracks);
                    
                    if strcmp(projAxis3D, 'three')
                        validOut_zx = validOut;
                        validOut_zy = validOut;
                        varargout_zx = varargout{i};%,varargout{i},varargout{i});
                        varargout_zy = varargout{i};%,varargout{i},varargout{i});
                        varargout_xy = varargout{i};%,varargout{i},varargout{i});
                    end
                    
                    for j=1:length(validOut)
                        % validOut(j).tracksCoordAmpCG = validOut(j).tracksCoordAmpCG(:,1:nCoords(j));
                        tT = validOut(j).tracksCoordAmpCG(:,1:nCoords(j));
                        % a1 = tT(:,4:8:end);
                        % dx1 = tT(:,5:8:end);
                        % dy1 = tT(:,6:8:end);
                        % dz1 = tT(:,7:8:end);
                        % da1 = tT(:,8:8:end);                        
                        switch projAxis3D
                            % [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
                            case 'Z'

                            case 'Y' %(ZX)
                                xCoords = tT(:,1:8:end);
                                yCoords = tT(:,2:8:end);
                                zCoords = tT(:,3:8:end);
                                
                                tT(:,1:8:end) = zCoords;
                                tT(:,2:8:end) = xCoords;
                                tT(:,3:8:end) = yCoords*ZXRatio;

                            case 'X' %(ZY)
                                xCoords = tT(:, 1:8:end);
                                yCoords = tT(:, 2:8:end);
                                zCoords = tT(:, 3:8:end);

                                tT(:, 1:8:end) = zCoords;
                                tT(:, 2:8:end) = yCoords;
                                tT(:, 3:8:end) = xCoords*ZXRatio;

                            case 'three'
                                zx_tT = tT;
                                zy_tT = tT;
                                
                                xCoords = tT(:,1:8:end);
                                yCoords = tT(:,2:8:end);
                                zCoords = tT(:,3:8:end);
                                
                                %ZX - BAD
                                zx_tT(:,1:8:end) = xCoords;
                                zx_tT(:,2:8:end) = zCoords+4+s1.MD.imSize_(1); %%%
                                zx_tT(:,3:8:end) = yCoords*ZXRatio;

                                %ZY - GOOD
                                zy_tT(:,1:8:end) = zCoords+4+s1.MD.imSize_(2);
                                zy_tT(:,2:8:end) = yCoords;
                                zy_tT(:,3:8:end) = xCoords*ZXRatio;

                                validOut_zy(j).tracksCoordAmpCG = zy_tT;
                                validOut_zx(j).tracksCoordAmpCG = zx_tT;

                            otherwise
                        end
                        validOut(j).tracksCoordAmpCG = tT;
                    end
                    if strcmp(projAxis3D,'three')
                        varargout_xy(validTracks) = validOut;
                        varargout_zx(validTracks) = validOut_zx;
                        varargout_zy(validTracks) = validOut_zy;
                        varargout{i} = vertcat(varargout_xy,varargout_zy,varargout_zx);                        
                    else
                        varargout{i}(validTracks) = validOut;
                    end
                end
                if  strcmp(output{i}, 'staticTracks')
                    validOut = varargout{i};
                    if strcmp(projAxis3D, 'three')
                        validOut_zx = validOut;
                        validOut_zy = validOut;
                        varargout_zx = varargout{i};%,varargout{i},varargout{i});
                        varargout_zy = varargout{i};%,varargout{i},varargout{i});
                        varargout_xy = varargout{i};%,varargout{i},varargout{i});
                    end
                    for j=1:length(validOut)
                        tT = validOut(j).tracksCoordAmpCG;
                        switch projAxis3D
                            % [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
                            case 'Z'
                            case 'Y'
                                xCoords = tT(:,1:8:end);
                                yCoords = tT(:,2:8:end);
                                zCoords = tT(:,3:8:end);
                                
                                tT(:,1:8:end) = zCoords;
                                tT(:,2:8:end) = xCoords;
                                tT(:,3:8:end) = yCoords*ZXRatio;

                            case 'X'
                                xCoords = tT(:, 1:8:end);
                                yCoords = tT(:, 2:8:end);
                                zCoords = tT(:, 3:8:end);

                                tT(:, 1:8:end) = zCoords;
                                tT(:, 2:8:end) = yCoords;
                                tT(:, 3:8:end) = xCoords*ZXRatio;

                            case 'three'
                                zx_tT = tT;
                                zy_tT = tT;
                                
                                xCoords = tT(:,1:8:end);
                                yCoords = tT(:,2:8:end);
                                zCoords = tT(:,3:8:end);
                                
                                %ZX - BAD
                                zx_tT(:,1:8:end) = xCoords;
                                zx_tT(:,2:8:end) = zCoords+4+s1.MD.imSize_(1); %%%
                                zx_tT(:,3:8:end) = yCoords*ZXRatio;

                                %ZY - GOOD
                                zy_tT(:,1:8:end) = zCoords+4+s1.MD.imSize_(2);
                                zy_tT(:,2:8:end) = yCoords;
                                zy_tT(:,3:8:end) = xCoords*ZXRatio;

                                validOut_zy(j).tracksCoordAmpCG = zy_tT;
                                validOut_zx(j).tracksCoordAmpCG = zx_tT;

                            otherwise
                        end
                        validOut(j).tracksCoordAmpCG = tT;
                    end
                    if strcmp(projAxis3D,'three')
                        varargout_xy = validOut;
                        varargout_zx = validOut_zx;
                        varargout_zy = validOut_zy;
                        varargout{i} = vertcat(varargout_xy,varargout_zy,varargout_zx);                        
                    else
                        varargout{i} = validOut;
                    end                    
%                     varargout{i} = validOut;
                end
            end
        end
        
        function output = getDrawableOutput(obj)
            colors = lines(numel(obj.owner_.channels_));
            output(1).name='Tracks';
            output(1).var='tracksFinalDynROIRef_oldFormat';
            output(1).type='overlay';
            output(1).defaultDisplayMethod = @(x) TracksDisplay(...
                'Color',colors(x,:));
            output(2).name='Gap length histogram';
            output(2).var='gapInfo';
            output(2).formatData=@(x) x(:,4);
            output(2).type='graph';
            output(2).defaultDisplayMethod=@(x)HistogramDisplay('XLabel',...
                'Gap length',...
                'YLabel','Counts');
            output(3).name='Static tracks';
            output(3).var='staticTracks';
            output(3).type='overlay';
            output(3).defaultDisplayMethod = @(x) TracksDisplay(...
                'Color',colors(x,:), 'useDragtail', false);                        

            output(1).formatData=@TrackingProcess.formatTracks2D;
            output(3).formatData=@TrackingProcess.formatTracks2D;
            
            if obj.funParams_.probDim == 3

                output(4).name='PlotTracks3D';
                output(4).var='plottracks3d';
                output(4).formatData=[];
                output(4).type='graph';
                % output(4).defaultDisplayMethod=@(x)plotTracks3DFigDisplay('plotFunc', @plotTracks3D);
                output(4).defaultDisplayMethod=@(x)FigDisplay('plotFunc', @plotTracks3D,...
                                                              'plotFunParams', {[], []});
            end
            
        end
    end
    
    methods(Static)
        function name = getName()
            name = 'Tracking in Dynamic ROI';
        end
        
        function h = GUI()
            h= @trackingProcessGUI;
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            nChan = numel(owner.channels_);
            funParams.ChannelIndex = 1:numel(owner.channels_);
            funParams.EstimateTrackability=false;
            funParams.processBuildDynROI=[]; % DynROI used for computation; % make first available DynROIProc selected&set on the GUI, even default is empty in the process class. edited on 2021-01-04
            funParams.buildDynROIProcessChannel=1; % Added for the setting GUI, but not used in the wrapper func.
            
            % should detect for which channels a detection process output exists.
            funParams.DetProcessIndex = []; % perhaps tag by process & channel
            
            funParams.OutputDirectory = [outputDir  filesep 'tracksInROI'];
            % --------------- time range ----------------
            funParams.timeRange = []; % empty implies entier movie
            % --------------- gapCloseParam ----------------
            funParams.gapCloseParam.timeWindow = 5; %IMPORTANT maximum allowed time gap (in frames) between a track segment end and a track segment start that allows linking them.
            funParams.gapCloseParam.mergeSplit = 0; % (SORT OF FLAG: 4 options for user) 1 if merging and splitting are to be considered, 2 if only merging is to be considered, 3 if only splitting is to be considered, 0 if no merging or splitting are to be considered.
            funParams.gapCloseParam.minTrackLen = 3; %minimum length of track segments from linking to be used in gap closing.
            funParams.gapCloseParam.diagnostics = 1; %FLAG 1 to plot a histogram of gap lengths in the end; 0 or empty otherwise.
            % --------------- kalmanFunctions ----------------
            kalmanFunctions = TrackingProcess.getKalmanFunctions(1);
            fields = fieldnames(kalmanFunctions);
            validFields = {'reserveMem','initialize','calcGain','timeReverse'};
            kalmanFunctions = rmfield(kalmanFunctions,fields(~ismember(fields,validFields)));
            funParams.kalmanFunctions = kalmanFunctions;
            % --------------- saveResults ----------------
            funParams.saveResults.export = 0; %FLAG allow additional export tracking results into matrix
            funParams.saveResults.exportTrackabilityData = 1; %FLAG allow exporting Kalman filter variable
            
            
            % --------------- Others ----------------
            funParams.verbose = 1;
            
            if owner.is3D
                funParams.probDim = 3;
            else
                funParams.probDim = 2;
            end
            
            funParams.costMatrices(1) = TrackingProcess.getDefaultLinkingCostMatrices(owner, funParams.gapCloseParam.timeWindow,1);
            funParams.costMatrices(2) = TrackingProcess.getDefaultGapClosingCostMatrices(owner, funParams.gapCloseParam.timeWindow,1);
            
            
            %list of parameters which can be specified at a per-channel
            %level. If specified as scalar these will  be replicated
            %             funParams.PerChannelParams = {'gapCloseParam','kalmanFunctions','costMatrices'};
            % funParams = prepPerChannelParams(funParams, nChan);
        end
    end
    
end