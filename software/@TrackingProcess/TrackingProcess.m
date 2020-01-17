classdef TrackingProcess < DataProcessingProcess & NonSingularProcess
    % A class definition for a generic tracking process.
    %
    % Chuangang Ren, 11/2010
    % Sebastien Besson (last modified Dec 2011)
    % Mark Kittisopikul, Nov 2014, Added channelOutput cache
    % Andrew R. Jamieson, Dec 2016, updated parameters for getDefaultGapClosingCostMatrices and GUI
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
        useCache_ = false;
    end
    
    methods(Access = public)
        
        function obj = TrackingProcess(owner, varargin)
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
                super_args{2} = TrackingProcess.getName;
                super_args{3} = @trackMovie;
                if isempty(funParams)
                    funParams = TrackingProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            obj = obj@DataProcessingProcess(super_args{:});
            obj.is3Dcompatible_ = true;
        end

        function h=draw(obj,iChan,varargin)
            h = obj.draw@DataProcessingProcess(iChan,varargin{:},'useCache',true);
        end
        
        function varargout = loadChannelOutput(obj, iChan, varargin)
            
            % Input check
            outputList = {'tracksFinal', 'gapInfo', 'staticTracks',...
                          'plottracks3d'};
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
                ZXRatio = obj.owner_.pixelSizeZ_/obj.owner_.pixelSize_;
            end

            s = cached.load(obj.outFilePaths_{1,iChan}, '-useCache', ip.Results.useCache, 'tracksFinal');

            varargout = cell(numel(output), 1);
        
            
            for i = 1:numel(output)
                switch output{i}
                    case {'tracksFinal', 'staticTracks', 'plottracks3d'}
                        varargout{i} = s.tracksFinal;
                    case 'gapInfo'
                        varargout{i} = findTrackGaps(s.tracksFinal);
                end
                if strcmp(output{i}, 'tracksFinal') && ~isempty(iFrame)
                    % Filter tracks existing in input frame
                    trackSEL=getTrackSEL(s.tracksFinal);
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
                                zx_tT(:,2:8:end) = zCoords+4+obj.owner_.imSize_(1); %%%
                                zx_tT(:,3:8:end) = yCoords*ZXRatio;

                                %ZY - GOOD
                                zy_tT(:,1:8:end) = zCoords+4+obj.owner_.imSize_(2);
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
                                zx_tT(:,2:8:end) = zCoords+4+obj.owner_.imSize_(1); %%%
                                zx_tT(:,3:8:end) = yCoords*ZXRatio;

                                %ZY - GOOD
                                zy_tT(:,1:8:end) = zCoords+4+obj.owner_.imSize_(2);
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
            output(1).var='tracksFinal';
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
            name = 'Tracking';
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
            
            
            % should detect for which channels a detection process output exists.
            funParams.DetProcessIndex = []; % perhaps tag by process & channel

            funParams.OutputDirectory = [outputDir  filesep 'tracks'];
            % --------------- time range ----------------
            funParams.timeRange = []; % empty implies entier movie
            % --------------- gapCloseParam ----------------
            funParams.gapCloseParam.timeWindow = 5; %IMPORTANT maximum allowed time gap (in frames) between a track segment end and a track segment start that allows linking them.
            funParams.gapCloseParam.mergeSplit = 0; % (SORT OF FLAG: 4 options for user) 1 if merging and splitting are to be considered, 2 if only merging is to be considered, 3 if only splitting is to be considered, 0 if no merging or splitting are to be considered.
            funParams.gapCloseParam.minTrackLen = 1; %minimum length of track segments from linking to be used in gap closing.
            funParams.gapCloseParam.diagnostics = 1; %FLAG 1 to plot a histogram of gap lengths in the end; 0 or empty otherwise.
            % --------------- kalmanFunctions ----------------
            kalmanFunctions = TrackingProcess.getKalmanFunctions(1);
            fields = fieldnames(kalmanFunctions);
            validFields = {'reserveMem','initialize','calcGain','timeReverse'};
            kalmanFunctions = rmfield(kalmanFunctions,fields(~ismember(fields,validFields)));
            funParams.kalmanFunctions = kalmanFunctions;            
            % --------------- saveResults ----------------
            funParams.saveResults.export = 0; %FLAG allow additional export of the tracking results into matrix
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
        
        function kalmanFunctions = getKalmanFunctions(index)
            % Brownian + Directed motion models
            kalmanFunctions(1).name = 'Brownian + Directed motion models';
            kalmanFunctions(1).reserveMem = func2str(@kalmanResMemLM);
            kalmanFunctions(1).initialize  = func2str(@kalmanInitLinearMotion);
            kalmanFunctions(1).initializeGUI  = @kalmanInitializationGUI;
            kalmanFunctions(1).calcGain    = func2str(@kalmanGainLinearMotion);
            kalmanFunctions(1).timeReverse = func2str(@kalmanReverseLinearMotion);
            
            % Microtubule plus-end dynamics
            kalmanFunctions(2).name = 'Microtubule plus-end dynamics';
            kalmanFunctions(2).reserveMem = func2str(@kalmanResMemLM);
            kalmanFunctions(2).initialize  = func2str(@plusTipKalmanInitLinearMotion);
            kalmanFunctions(2).initializeGUI  = @kalmanInitializationGUI;
            kalmanFunctions(2).calcGain    = func2str(@plusTipKalmanGainLinearMotion);
            kalmanFunctions(2).timeReverse = func2str(@kalmanReverseLinearMotion);

            if nargin>0
                assert(all(ismember(index, 1:numel(kalmanFunctions))));
                kalmanFunctions=kalmanFunctions(index);
            end
        end
        
        function costMatrix = getDefaultLinkingCostMatrices(owner,timeWindow,varargin)
            
            % Brownian + Directed motion models
            costMatrices(1).name = 'Brownian + Directed motion models';
            costMatrices(1).funcName = func2str(@costMatRandomDirectedSwitchingMotionLink);
            costMatrices(1).GUI = @costMatRandomDirectedSwitchingMotionLinkGUI;
            costMatrices(1).parameters.linearMotion = 0; % use linear motion Kalman filter.
            costMatrices(1).parameters.minSearchRadius = 2; %minimum allowed search radius. The search radius is calculated on the spot in the code given a feature's motion parameters. If it happens to be smaller than this minimum, it will be increased to the minimum.
            costMatrices(1).parameters.maxSearchRadius = 5; %IMPORTANT maximum allowed search radius. Again, if a feature's calculated search radius is larger than this maximum, it will be reduced to this maximum.
            costMatrices(1).parameters.brownStdMult = 3; %multiplication factor to calculate search radius from standard deviation.
            costMatrices(1).parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the linking (initial tracking) step.
            costMatrices(1).parameters.nnWindow = timeWindow; %number of frames before the current one where you want to look to see a feature's nearest neighbor in order to decide how isolated it is (in the initial linking step).
            costMatrices(1).parameters.kalmanInitParam = []; %Kalman filter initialization parameters.
            costMatrices(1).parameters.diagnostics = owner.nFrames_-1;
            
            % plusTip markers
            plusTipCostMatrix.name = 'Microtubule plus-end dynamics';
            plusTipCostMatrix.funcName = func2str(@plusTipCostMatLinearMotionLink);
            plusTipCostMatrix.GUI = @plusTipCostMatLinearMotionLinkGUI;
            plusTipCostMatrix.parameters.linearMotion = 1; % use linear motion Kalman filter.
            plusTipCostMatrix.parameters.minSearchRadius = 2; %minimum allowed search radius. The search radius is calculated on the spot in the code given a feature's motion parameters. If it happens to be smaller than this minimum, it will be increased to the minimum.
            plusTipCostMatrix.parameters.maxSearchRadius = 10; %IMPORTANT maximum allowed search radius. Again, if a feature's calculated search radius is larger than this maximum, it will be reduced to this maximum.
            plusTipCostMatrix.parameters.brownStdMult = 3; %multiplication factor to calculate search radius from standard deviation.
            plusTipCostMatrix.parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the linking (initial tracking) step.
            plusTipCostMatrix.parameters.nnWindow = timeWindow; %number of frames before the current one where you want to look to see a feature's nearest neighbor in order to decide how isolated it is (in the initial linking step).
            plusTipCostMatrix.parameters.kalmanInitParam.initVelocity = []; %Kalman filter initialization parameters.
            plusTipCostMatrix.parameters.kalmanInitParam.convergePoint = []; %Kalman filter initialization parameters.
            plusTipCostMatrix.parameters.kalmanInitParam.searchRadiusFirstIteration = 20; %Kalman filter initialization parameters.
            plusTipCostMatrix.parameters.diagnostics = [];
            costMatrices(2)=plusTipCostMatrix;

            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addRequired('timeWindow',@isscalar);
            ip.addOptional('index',1:length(costMatrices),@isvector);
            ip.parse(owner,timeWindow,varargin{:});
            index = ip.Results.index;
            costMatrix=costMatrices(index);
        end
        
        function costMatrix = getDefaultGapClosingCostMatrices(owner,timeWindow,varargin)
            
            % Linear motion
            costMatrices(1).name = 'Brownian + Directed motion models';
            costMatrices(1).funcName = func2str(@costMatRandomDirectedSwitchingMotionCloseGaps);
            costMatrices(1).GUI = @costMatRandomDirectedSwitchingMotionCloseGapsGUI;
            costMatrices(1).parameters.linearMotion = 0; %use linear motion Kalman filter.
            
            costMatrices(1).parameters.minSearchRadius = 2; %minimum allowed search radius.
            costMatrices(1).parameters.maxSearchRadius = 5; %maximum allowed search radius.
            costMatrices(1).parameters.brownStdMult = 3*ones(timeWindow,1); %multiplication factor to calculate Brownian search radius from standard deviation.
            
            costMatrices(1).parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the gap closing and merging/splitting step.
            costMatrices(1).parameters.nnWindow = timeWindow; %number of frames before/after the current one where you want to look for a track's nearest neighbor at its end/start (in the gap closing step).
            costMatrices(1).parameters.brownScaling = [0.5 0.01]; %power for scaling the Brownian search radius with time, before and after timeReachConfB (next parameter).
            costMatrices(1).parameters.timeReachConfB = timeWindow; %before timeReachConfB, the search radius grows with time with the power in brownScaling(1); after timeReachConfB it grows with the power in brownScaling(2).
            costMatrices(1).parameters.ampRatioLimit = [0.5 2]; % (FLAG + VALUES small-big value) for merging and splitting. Minimum and maximum ratios between the intensity of a feature after merging/before splitting and the sum of the intensities of the 2 features that merge/split.
            
            % If parameters.linearMotion = 1
            costMatrices(1).parameters.lenForClassify = 5; %minimum track segment length to classify it as linear or random.
            costMatrices(1).parameters.linStdMult = 3*ones(timeWindow,1); %multiplication factor to calculate linear search radius from standard deviation.
            costMatrices(1).parameters.linScaling = [1 0.01]; %power for scaling the linear search radius with time (similar to brownScaling).
            costMatrices(1).parameters.timeReachConfL = timeWindow; %similar to timeReachConfB, but for the linear part of the motion.
            costMatrices(1).parameters.maxAngleVV = 30; %maximum angle between the directions of motion of two tracks that allows linking them (and thus closing a gap). Think of it as the equivalent of a searchRadius but for angles.
            % ---------------------------------
            
            costMatrices(1).parameters.gapPenalty = 1.5; %penalty for increasing temporary disappearance time (disappearing for n frames gets a penalty of gapPenalty^n).
            costMatrices(1).parameters.resLimit = []; % text field resolution limit, which is generally equal to 3 * point spread function sigma.
            % ---------------------------------
            % NEW PARAMETERS by Khuloud - Dec 2016
            costMatrices(1).parameters.gapExcludeMS = 0; %(1)%flag to allow gaps to exclude merges and splits
            costMatrices(1).parameters.strategyBD = 0; %(-1)%strategy to calculate birth and death cost
            % ---------------------------------

            % Linear motion
            plusTipCostMatrix.name = 'Microtubule plus-end dynamics';
            plusTipCostMatrix.funcName = func2str(@plusTipCostMatCloseGaps);
            plusTipCostMatrix.GUI = @plusTipCostMatCloseGapsGUI;
            plusTipCostMatrix.parameters.maxFAngle = 30; %use linear motion Kalman filter.
            plusTipCostMatrix.parameters.maxBAngle = 10; %use linear motion Kalman filter.
            plusTipCostMatrix.parameters.backVelMultFactor = 1.5;
            plusTipCostMatrix.parameters.fluctRad = 1.0;
            plusTipCostMatrix.parameters.breakNonLinearTracks = false;
            costMatrices(2)=plusTipCostMatrix;
            
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addRequired('timeWindow',@isscalar);
            ip.addOptional('index',1:length(costMatrices),@isvector);
            ip.parse(owner,timeWindow,varargin{:});
            index = ip.Results.index;
            costMatrix=costMatrices(index);
        end
        
        function displayTracks = formatTracks(tracks)
            displayTracks = TrackingProcess.formatTracks2D(tracks);
        end
        
        function displayTracks = formatTracks2D(tracks)
            % Format tracks structure into compound tracks for display
            % purposes
            
            % Determine the number of compound tracks
            nCompoundTracks = cellfun('size',{tracks.tracksCoordAmpCG},1)';
            trackIdx = 1:length(tracks);
            
            % Filter by the tracks that are nonzero
            filter = nCompoundTracks > 0;
            tracks = tracks(filter);
            trackIdx = trackIdx(filter);
            nCompoundTracks = nCompoundTracks(filter);
            
            % Get the track lengths (nFrames x 8)
            trackLengths = cellfun('size',{tracks.tracksCoordAmpCG},2)';
            % Unique track lengths for batch processing later
            uTrackLengths = unique(trackLengths);
            % Running total of displayTracks for indexing
            nTracksTot = [0 cumsum(nCompoundTracks(:))'];
            % Total number of tracks
            nTracks = nTracksTot(end);
            
            % Fail fast if no track
            if nTracks == 0
                displayTracks = struct.empty(1,0);
                return
            end
            
            % Number of events in each seqOfEvents for indexing
            % Each will normally have 2 events, beginning and end
            nEvents = cellfun('size',{tracks.seqOfEvents},1);

            % Initialize displayTracks structure
            % xCoord: x-coordinate of simple track
            % yCoord: y-coordinate of simple track
            % events (deprecated): split or merge
            % number: number corresponding to the original input compound
            % track number
            % splitEvents: times when splits occur
            % mergeEvents: times when merges occur
            displayTracks(nTracks,1) = struct('xCoord', [], 'yCoord', [], 'number', [], 'splitEvents', [], 'mergeEvents', []);
            
            hasSeqOfEvents = isfield(tracks,'seqOfEvents');
            hasLabels = isfield(tracks, 'label');
            
            if(hasLabels)
                labels = vertcat(tracks.label);
                if(size(labels,1) == nTracks)
                    hasPerSegmentLabels = true;
                end
            end

            % Batch by unique trackLengths
            for trackLength = uTrackLengths'
                %% Select original track numbers
                selection = trackLengths == trackLength;
                sTracks = tracks(selection);
                sTrackCoordAmpCG = vertcat(sTracks.tracksCoordAmpCG);
                % track number relative to input struct array
                sTrackIdx = trackIdx(selection);
                % index of selected tracks
                siTracks = nTracksTot(selection);
                
                % runs in current selection
                snCompoundTracks = nCompoundTracks(selection);
                snTracksTot = [0 ; cumsum(snCompoundTracks)];
                
                % decode run lengths
                snTracks = snTracksTot(end);
                idx = zeros(snTracks,1);
                idx(snTracksTot(1:end-1) + 1) = ones(size(snCompoundTracks));
                idx = cumsum(idx);

                % absolute index of tracks
                iTracks = ones(snTracks,1);
                iTracks(snTracksTot(1:end-1) + 1) = [siTracks(1); diff(siTracks)'] - [0 ; snCompoundTracks(1:end-1)-1];
                iTracks = cumsum(iTracks) + 1;
                
                % grab x and y coordinate matrices
                xCoords = sTrackCoordAmpCG(:,1:8:end);
                yCoords = sTrackCoordAmpCG(:,2:8:end);
                
                %% Process sequence of events
                if(hasSeqOfEvents)
                    % make sequence of events matrix
                    seqOfEvents = vertcat(sTracks.seqOfEvents);
                    nSelectedEvents = nEvents(selection);
                    iStartEvents = [1 cumsum(nSelectedEvents(1:end-1))+1];

                    % The fifth column is the start frame for each track
                    seqOfEvents(iStartEvents,5) = [seqOfEvents(1) ; diff(seqOfEvents(iStartEvents,1))];
                    seqOfEvents(:,5) = cumsum(seqOfEvents(:,5));

                    % The sixth column is the offset for the current selected
                    % tracks
                    seqOfEvents(iStartEvents,6) = [0; snCompoundTracks(1:end-1)];
                    seqOfEvents(:,6) = cumsum(seqOfEvents(:,6));

                    % Isolate merges and splits
                    seqOfEvents = seqOfEvents(~isnan(seqOfEvents(:,4)),:);

                    % Apply offset 
                    seqOfEvents(:,3) = seqOfEvents(:,3) + seqOfEvents(:,6);
                    seqOfEvents(:,4) = seqOfEvents(:,4) + seqOfEvents(:,6);

                    % Number of Frames
                    nFrames = trackLength/8;

                    %% Splits
                    % The 2nd column indicates split (1) or merge(2)
                    splitEvents = seqOfEvents(seqOfEvents(:,2) == 1,:);
                    % Evaluate time relative to start of track
                    splitEventTimes = splitEvents(:,1) - splitEvents(:,5);
                    % Time should not exceed the number of coordinates we have
                    splitEvents = splitEvents(splitEventTimes < nFrames,:);
                    splitEventTimes = splitEventTimes(splitEventTimes < nFrames,:);
                    iTrack1 = splitEvents(:,3);
                    iTrack2 = splitEvents(:,4);

                    % Use accumarray to gather the splitEventTimes into cell
                    % arrays
                    if(~isempty(splitEventTimes))
                        splitEventTimeCell = accumarray([iTrack1 ; iTrack2], [splitEventTimes ; splitEventTimes ],[size(xCoords,1) 1],@(x) {x'},{});
                    else
                        splitEventTimeCell = cell(size(xCoords,1),1);
                    end

                   %  leftIdx = sub2ind(size(xCoords),iTrack1,splitEventTimes);
                   %  rightIdx = sub2ind(size(xCoords),iTrack2,splitEventTimes);
                    leftIdx = (splitEventTimes-1)*size(xCoords,1) + iTrack1;
                    rightIdx = (splitEventTimes-1)*size(xCoords,1) + iTrack2;
                    xCoords(leftIdx) = xCoords(rightIdx);
                    yCoords(leftIdx) = yCoords(rightIdx);

                    %% Merges
                    % The 2nd column indicates split (1) or merge(2)
                    mergeEvents = seqOfEvents(seqOfEvents(:,2) == 2,:);
                    mergeEventTimes = mergeEvents(:,1) - mergeEvents(:,5);
                    mergeEvents = mergeEvents(mergeEventTimes < nFrames,:);
                    mergeEventTimes = mergeEventTimes(mergeEventTimes < nFrames,:);
                    mergeEventTimes = mergeEventTimes + 1;
                    iTrack1 = mergeEvents(:,3);
                    iTrack2 = mergeEvents(:,4);

                    if(~isempty(mergeEventTimes))
                        mergeEventTimeCell = accumarray([iTrack1 ; iTrack2], [mergeEventTimes ; mergeEventTimes ],[size(xCoords,1) 1],@(x) {x'},{});
                    else
                        mergeEventTimeCell = cell(size(xCoords,1),1);
                    end

                    % leftIdx = sub2ind(size(xCoords),iTrack1,mergeEventTimes);
                    % rightIdx = sub2ind(size(xCoords),iTrack2,mergeEventTimes);
                    leftIdx = (mergeEventTimes-1)*size(xCoords,1) + iTrack1;
                    rightIdx = (mergeEventTimes-1)*size(xCoords,1) + iTrack2;
                    xCoords(leftIdx) = xCoords(rightIdx);
                    yCoords(leftIdx) = yCoords(rightIdx);
                end
                          
                %% Load cells into struct fields
                for i=1:length(iTracks)
                    iTrack = iTracks(i);
                    displayTracks(iTrack).xCoord = xCoords(i,:);
                    displayTracks(iTrack).yCoord = yCoords(i,:);
                    displayTracks(iTrack).number = sTrackIdx(idx(i));
                end            
                
                if(hasSeqOfEvents)
                    [displayTracks(iTracks).splitEvents] = splitEventTimeCell{:};
                    [displayTracks(iTracks).mergeEvents] = mergeEventTimeCell{:};
                end
                
                if hasLabels
                    if hasPerSegmentLabels
                        for i=1:length(iTracks)
                            iTrack = iTracks(i);
                            displayTracks(iTrack).label = labels(iTrack);
                        end    
                    else
                        [displayTracks(iTracks).label] = sTracks(idx).label;
                    end
                end
            end
        end
        
        function displayTracks = formatTracks3D(tracks)
            % Format tracks structure into compound tracks for display
            % purposes
            
            % Determine the number of compound tracks
            nCompoundTracks = cellfun('size',{tracks.tracksCoordAmpCG},1)';
            trackIdx = 1:length(tracks);
            
            % Filter by the tracks that are nonzero
            filter = nCompoundTracks > 0;
            tracks = tracks(filter);
            trackIdx = trackIdx(filter);
            nCompoundTracks = nCompoundTracks(filter);
            
            % Get the track lengths (nFrames x 8)
            trackLengths = cellfun('size',{tracks.tracksCoordAmpCG},2)';
            % Unique track lengths for batch processing later
            uTrackLengths = unique(trackLengths);
            % Running total of displayTracks for indexing
            nTracksTot = [0 cumsum(nCompoundTracks(:))'];
            % Total number of tracks
            nTracks = nTracksTot(end);
            
            % Fail fast if no track
            if nTracks == 0
                displayTracks = struct.empty(1,0);
                return
            end
            
            % Number of events in each seqOfEvents for indexing
            % Each will normally have 2 events, beginning and end
            nEvents = cellfun('size',{tracks.seqOfEvents},1);

            % Initialize displayTracks structure
            % xCoord: x-coordinate of simple track
            % yCoord: y-coordinate of simple track
            % events (deprecated): split or merge
            % number: number corresponding to the original input compound
            % track number
            % splitEvents: times when splits occur
            % mergeEvents: times when merges occur
            displayTracks(nTracks,1) = struct('xCoord', [], 'yCoord', [], 'number', [], 'splitEvents', [], 'mergeEvents', []);
            
            hasSeqOfEvents = isfield(tracks,'seqOfEvents');
            hasLabels = isfield(tracks, 'label');
            
            if(hasLabels)
                labels = vertcat(tracks.label);
                if(size(labels,1) == nTracks)
                    hasPerSegmentLabels = true;
                end
            end

            % Batch by unique trackLengths
            for trackLength = uTrackLengths'
                %% Select original track numbers
                selection = trackLengths == trackLength;
                sTracks = tracks(selection);
                sTrackCoordAmpCG = vertcat(sTracks.tracksCoordAmpCG);
                % track number relative to input struct array
                sTrackIdx = trackIdx(selection);
                % index of selected tracks
                siTracks = nTracksTot(selection);
                
                % runs in current selection
                snCompoundTracks = nCompoundTracks(selection);
                snTracksTot = [0 ; cumsum(snCompoundTracks)];
                
                % decode run lengths
                snTracks = snTracksTot(end);
                idx = zeros(snTracks,1);
                idx(snTracksTot(1:end-1) + 1) = ones(size(snCompoundTracks));
                idx = cumsum(idx);

                % absolute index of tracks
                iTracks = ones(snTracks,1);
                iTracks(snTracksTot(1:end-1) + 1) = [siTracks(1); diff(siTracks)'] - [0 ; snCompoundTracks(1:end-1)-1];
                iTracks = cumsum(iTracks) + 1;
                
                % grab x and y coordinate matrices
                xCoords = sTrackCoordAmpCG(:,1:8:end);
                yCoords = sTrackCoordAmpCG(:,2:8:end);
                
                %% Process sequence of events
                if(hasSeqOfEvents)
                    % make sequence of events matrix
                    seqOfEvents = vertcat(sTracks.seqOfEvents);
                    nSelectedEvents = nEvents(selection);
                    iStartEvents = [1 cumsum(nSelectedEvents(1:end-1))+1];

                    % The fifth column is the start frame for each track
                    seqOfEvents(iStartEvents,5) = [seqOfEvents(1) ; diff(seqOfEvents(iStartEvents,1))];
                    seqOfEvents(:,5) = cumsum(seqOfEvents(:,5));

                    % The sixth column is the offset for the current selected
                    % tracks
                    seqOfEvents(iStartEvents,6) = [0; snCompoundTracks(1:end-1)];
                    seqOfEvents(:,6) = cumsum(seqOfEvents(:,6));

                    % Isolate merges and splits
                    seqOfEvents = seqOfEvents(~isnan(seqOfEvents(:,4)),:);

                    % Apply offset 
                    seqOfEvents(:,3) = seqOfEvents(:,3) + seqOfEvents(:,6);
                    seqOfEvents(:,4) = seqOfEvents(:,4) + seqOfEvents(:,6);

                    % Number of Frames
                    nFrames = trackLength/8;

                    %% Splits
                    % The 2nd column indicates split (1) or merge(2)
                    splitEvents = seqOfEvents(seqOfEvents(:,2) == 1,:);
                    % Evaluate time relative to start of track
                    splitEventTimes = splitEvents(:,1) - splitEvents(:,5);
                    % Time should not exceed the number of coordinates we have
                    splitEvents = splitEvents(splitEventTimes < nFrames,:);
                    splitEventTimes = splitEventTimes(splitEventTimes < nFrames,:);
                    iTrack1 = splitEvents(:,3);
                    iTrack2 = splitEvents(:,4);

                    % Use accumarray to gather the splitEventTimes into cell
                    % arrays
                    if(~isempty(splitEventTimes))
                        splitEventTimeCell = accumarray([iTrack1 ; iTrack2], [splitEventTimes ; splitEventTimes ],[size(xCoords,1) 1],@(x) {x'},{});
                    else
                        splitEventTimeCell = cell(size(xCoords,1),1);
                    end

                    % leftIdx = sub2ind(size(xCoords),iTrack1,splitEventTimes);
                    % rightIdx = sub2ind(size(xCoords),iTrack2,splitEventTimes);
                    leftIdx = (splitEventTimes-1)*size(xCoords,1) + iTrack1;
                    rightIdx = (splitEventTimes-1)*size(xCoords,1) + iTrack2;
                    xCoords(leftIdx) = xCoords(rightIdx);
                    yCoords(leftIdx) = yCoords(rightIdx);

                    %% Merges
                    % The 2nd column indicates split (1) or merge(2)
                    mergeEvents = seqOfEvents(seqOfEvents(:,2) == 2,:);
                    mergeEventTimes = mergeEvents(:,1) - mergeEvents(:,5);
                    mergeEvents = mergeEvents(mergeEventTimes < nFrames,:);
                    mergeEventTimes = mergeEventTimes(mergeEventTimes < nFrames,:);
                    mergeEventTimes = mergeEventTimes + 1;
                    iTrack1 = mergeEvents(:,3);
                    iTrack2 = mergeEvents(:,4);

                    if(~isempty(mergeEventTimes))
                        mergeEventTimeCell = accumarray([iTrack1 ; iTrack2], [mergeEventTimes ; mergeEventTimes ],[size(xCoords,1) 1],@(x) {x'},{});
                    else
                        mergeEventTimeCell = cell(size(xCoords,1),1);
                    end

                    % leftIdx = sub2ind(size(xCoords),iTrack1,mergeEventTimes);
                    % rightIdx = sub2ind(size(xCoords),iTrack2,mergeEventTimes);
                    leftIdx = (mergeEventTimes-1)*size(xCoords,1) + iTrack1;
                    rightIdx = (mergeEventTimes-1)*size(xCoords,1) + iTrack2;
                    xCoords(leftIdx) = xCoords(rightIdx);
                    yCoords(leftIdx) = yCoords(rightIdx);
                end
                          
                %% Load cells into struct fields
                for i=1:length(iTracks)
                    iTrack = iTracks(i);
                    displayTracks(iTrack).xCoord = xCoords(i,:);
                    displayTracks(iTrack).yCoord = yCoords(i,:);
                    displayTracks(iTrack).number = sTrackIdx(idx(i));
                end            
                
                if(hasSeqOfEvents)
                    [displayTracks(iTracks).splitEvents] = splitEventTimeCell{:};
                    [displayTracks(iTracks).mergeEvents] = mergeEventTimeCell{:};
                end
                
                if hasLabels
                    if hasPerSegmentLabels
                        for i=1:length(iTracks)
                            iTrack = iTracks(i);
                            displayTracks(iTrack).label = labels(iTrack);
                        end    
                    else
                        [displayTracks(iTracks).label] = sTracks(idx).label;
                    end
                end
            end
        end
    end
end

