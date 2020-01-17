classdef MotionAnalysisProcess < PostTrackingProcess
    % A concrete class for analyzing tracks diffusion
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
    
    % Sebastien Besson, March 2012
    % Mark Kittisopikul, Nov 2014, Added channelOutput cache
    
    methods (Access = public)
        function obj = MotionAnalysisProcess(owner, varargin)
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
                
                super_args{1} = owner;
                super_args{2} = MotionAnalysisProcess.getName;
                super_args{3} = @analyzeMovieMotion;
                if isempty(funParams)  % Default funParams
                    funParams = MotionAnalysisProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            
            obj = obj@PostTrackingProcess(super_args{:});
        end

        function h=draw(obj,iChan,varargin)
            h = obj.draw@PostTrackingProcess(iChan,varargin{:},'useCache',true);
        end

        function varargout = loadChannelOutput(obj,iChan,varargin)
            
            % Input check
            outputList = {'diffAnalysisRes', 'tracks'};
            ip =inputParser;
            ip.addRequired('iChan',@(x) isscalar(x) && obj.checkChanNum(x));
            ip.addOptional('iFrame',[],@(x) isempty(x) || isscalar(x) && obj.checkFrameNum(x));
            ip.addParamValue('iZ',[], @(x)ismember(x,1:obj.owner_.zSize_));
            ip.addParamValue('useCache',false,@islogical);
            ip.addParamValue('output',outputList,@(x) all(ismember(x,outputList)));
            ip.addParameter('projectionAxis3D','Z', @(x) ismember(x,{'Z','X','Y'}));
            ip.parse(iChan,varargin{:})
            iFrame = ip.Results.iFrame;
            output = ip.Results.output;
            projAxis3D = ip.Results.projectionAxis3D;
            if ischar(output),output={output}; end
            nOutput = numel(output);
            
            % Data loading

            if ~strcmp(projAxis3D,'Z') %&& obj.owner_.is3D % && ~obj.funParams_.isoOutput
                ZXRatio = obj.owner_.pixelSizeZ_/obj.owner_.pixelSize_;
            end


            % load outFilePaths_{1,iChan}
            s = cached.load(obj.outFilePaths_{1,iChan}, '-useCache', ip.Results.useCache, output{:});
            
            varargout = cell(nOutput);
            for i = 1:nOutput
                switch output{i}
                    case 'tracks'
                        tracksFinal = s.(output{i});
                        if ~isempty(iFrame),
                            % Filter tracks existing in input frame
                            trackSEL=getTrackSEL(tracksFinal);
                            validTracks = (iFrame>=trackSEL(:,1) &iFrame<=trackSEL(:,2));
                            [tracksFinal(~validTracks).tracksCoordAmpCG]=deal([]);
                            
                            for j=find(validTracks)'

                                tracksFinal(j).tracksCoordAmpCG = tracksFinal(j).tracksCoordAmpCG(:,1:8*(iFrame-trackSEL(j,1)+1));
                                tT = tracksFinal(j).tracksCoordAmpCG;

                                switch projAxis3D
                                    % [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
                                    case 'Z'
                                    case 'Y'
                                        xCoords = tT(:,1:8:end);
                                        yCoords = tT(:,2:8:end);
                                        zCoords = tT(:,3:8:end)*ZXRatio;
                                        tT(:,1:8:end) = zCoords;
                                        tT(:,2:8:end) = xCoords;
                                        tT(:,3:8:end) = yCoords;
                                    case 'X'
                                        xCoords = tT(:,1:8:end);
                                        yCoords = tT(:,2:8:end);
                                        zCoords = tT(:,3:8:end)*ZXRatio;
                                        tT(:,1:8:end) = zCoords;
                                        tT(:,2:8:end) = yCoords;
                                        tT(:,3:8:end) = xCoords;
                                    otherwise
                                end
                                tracksFinal(j).tracksCoordAmpCG = tT;                            

                            end
                            varargout{i} = tracksFinal;
                        else
                            varargout{i} = tracksFinal;
                        end
                    case 'diffAnalysisRes'
                        varargout{i} = s.(output{i});
                end
            end
        end
        
        function output = getDrawableOutput(obj)
            types = MotionAnalysisProcess.getTrackTypes();
            colors = vertcat(types.color);
            output(1).name='Classified tracks';
            output(1).var='tracks';
            output(1).formatData=@MotionAnalysisProcess.formatTracks;
            output(1).type='overlay';
            output(1).defaultDisplayMethod=@(x)TracksDisplay('Color', colors);
        end
        
    end
    methods (Static)
        
        function name = getName()
            name = 'Track Analysis';
        end
        function h = GUI()
            h = @motionAnalysisProcessGUI;
        end
        
        function alpha = getAlphaValues()
            alpha=[0.01 0.05 0.1 0.2];
        end
        
        function methods = getConfinementRadiusMethods()
            methods(1).type = 0;
            methods(1).name = 'Mean positional standard deviation';
            methods(2).type = 1;
            methods(2).name = 'Minimum positional standard deviation';
            methods(3).type = 2;
            methods(3).name = 'Rectangle approximation';
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelIndex = 1 : numel(owner.channels_);
            funParams.OutputDirectory = [outputDir  filesep 'MotionAnalysis'];

            if owner.is3D
                funParams.probDim = 3;
            else
                funParams.probDim = 2;
            end
            
            funParams.checkAsym = 0;
            funParams.alphaValues = [0.05 0.1];
            funParams.confRadMin=0;
            funParams.driftCorrect = 0;
        end
        
        function displayTracks = formatTracks(tracks)
            % Format classified tracks into structure for display
            
            % Read track types and classification matrix
            types = MotionAnalysisProcess.getTrackTypes();
            % Matrix of the number of track segments by 3
            track_class = vertcat(tracks.classification);

            % Number of labels per track needed
            nLabels = cellfun('size',{tracks.classification},1);
            % Labels is a cell array of indices corresponding to to track segment types,
            %   see getTracksTypes
            % Initialize all labels as unlabeled
            labels = arrayfun(@(x) ones(x,1)*numel(types),nLabels,'UniformOutput',false);
            % Map indicates position of last label for each track
            map = cumsum(nLabels);
            
            % logical array per track index of if there is more than one label per track
            nLabels_gt_1 = nLabels > 1;
            % labels index where the labels for each track starts
            %  if there is more than one label per track
            start = map(nLabels_gt_1)-nLabels(nLabels_gt_1)+1;
            % labels index where the labels for each track ends
            %  if there is more than one label per track
            finish = map(nLabels_gt_1);

            % Set labels as needed
            for i = 1 : numel(types) - 1
                % idx is a logical array if each track _segment_ belongs to types(i)
                idx = types(i).f(track_class);
                % idx2 is a cell array of logical arrays
                %  the number of cells corresponds to each track index
                %  the index of the logical array in each cell refers to each track segment
                % Setup logical arrays. This works for when nLabels == 1
                idx2 = num2cell(idx(map));
                % deal with nLabels > 1 separately, grab range of labels for each track segment
                %  corresponding to each track        
                idx2(nLabels_gt_1) = arrayfun(@(s,e) idx(s:e),start,finish,'UniformOutput',false);
                % idx is now a cell array of logical arrays marking for each track
                %  if each track segment belongs to types(i)
                idx = idx2;
                % Select only the indices where at least one segment belongs to type(i)
                any_idx = cellfun(@any,idx);
                % Assign label as i for each track segment belonging to types(i)
                % Merge with previous labels
                labels(any_idx) = cellfun(@(idx_i,labels_i) labels_i.*~idx_i + i.*idx_i, ...
                    idx(any_idx), labels(any_idx)','UniformOutput',false);
            end
            % Assign labels to each track
            [tracks.label] = deal(labels{:});

            % Format tracks using TrackingProcess utility function
            displayTracks = TrackingProcess.formatTracks(tracks);
        end
        
        function types = getTrackTypes()
            % Get the color map for classified tracks
            %
            % see also: plotTracksDiffAnalysis2D
            % Immobile: brown
            types(1).name = 'immobile';
            types(1).f = @(x) x(:, 1) ~= 1 & x(:, 2) == 0;
            types(1).color = [0.5 0.3 0];
            % Linear 1D confined: orange
            types(2).name = 'linear & 1D confined diffusion';
            types(2).f = @(x) x(:, 1) == 1 & x(:, 3) == 1;
            types(2).color = [1 0.7 0];
            % Linear 1D normal: bright red 
            types(3).name = 'linear & 1D normal diffusion';
            types(3).f = @(x) x(:, 1) == 1 & x(:, 3) == 2;
            types(3).color = [1 0 0];
            % Linear 1D super: bright green
            types(4).name = 'linear & 1D super diffusion';
            types(4).f = @(x) x(:, 1) == 1 & x(:, 3) == 3;
            types(4).color = [0 1 0];
            % Linear 1D too short: yellow
            types(5).name = 'linear & too short to analyze 1D diffusion';
            types(5).f = @(x) x(:, 1) == 1 & isnan(x(:, 3));
            types(5).color = [1 1 0];
            % Random/unclassified & 2D confined: blue
            types(6).name = 'random/unclassified & 2D confined diffusion';
            types(6).f = @(x) x(:, 1) ~= 1 & x(:, 2) == 1;
            types(6).color = [0 0 1];
            % Random/unclassified & 2D normal: cyan
            types(7).name = 'random/unclassified & 2D normal diffusion';
            types(7).f = @(x) x(:, 1) ~= 1 & x(:, 2) == 2;
            types(7).color = [0 1 1];
            % Random/unclassified & 2D super: magenta
            types(8).name = 'random/unclassified & 2D super diffusion';
            types(8).f = @(x) x(:, 1) ~= 1 & x(:, 2) == 3;
            types(8).color = [1 0 1];
            % Random & 2D too short: purple
            types(9).name = 'random & too short to analyze 2D diffusion';
            types(9).f = @(x) x(:, 1) == 0 & isnan(x(:, 2));
            types(9).color = [.6 0 1];
            % Too short: grey
            types(10).name = 'too short for any analysis';
            types(10).f = @(x) 1;
            types(10).color = [.7 .7 .7];
        end
        
        function hIm = showTrackTypes(newFigure)
            % MotionAnalysisProcess.showTypes: Show types and their colors
            
            if(nargin < 1)
                newFigure = false;
            end
            if(newFigure)
                figure;
            end
            
            types = MotionAnalysisProcess.getTrackTypes;
            
            % Create a grey background
            hIm = imshow(ones(300)*0.3);
            % Display type names in their color from top to bottom
            for i=1:length(types)
                text(1,length(types)-i+1,types(i).name, ...
                    'Color',types(i).color, ...
                    'Units','characters');
            end
        end
        
    end
end
