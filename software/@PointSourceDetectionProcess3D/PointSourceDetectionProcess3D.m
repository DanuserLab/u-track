classdef PointSourceDetectionProcess3D < DetectionProcess
%PointSourceDetectionProcess3D is a concrete class of a point source
%detection process for 3d
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
        function obj = PointSourceDetectionProcess3D(owner, outputDir, funParams)
            % Constructor of the SubResolutionProcess
            super_args{1} = owner;
            super_args{2} = PointSourceDetectionProcess3D.getName;
            super_args{3} = @detectMoviePointSources3D;
            
            if nargin < 3 || isempty(funParams)  % Default funParams
                if nargin <2, outputDir = owner.outputDirectory_; end
                funParams = PointSourceDetectionProcess3D.getDefaultParams(owner,outputDir);
            end
            
            super_args{4} = funParams;
            
            obj = obj@DetectionProcess(super_args{:});
            obj.is3Dcompatible_ = true;

        end
        
        function varargout = loadChannelOutput(obj,iChan,varargin)
            
            % Input check
            outputList = {'movieInfo', 'detect3D','detect3Dall', 'detectionsLabRef'};
            ip = inputParser;
            ip.addRequired('iChan',@(x) ismember(x,1:numel(obj.owner_.channels_)));
            ip.addOptional('iFrame',1:obj.owner_.nFrames_,...
                @(x) ismember(x,1:obj.owner_.nFrames_));
            ip.addParameter('useCache',true,@islogical);
            ip.addParameter('iZ',[], @(x) ismember(x,1:obj.owner_.zSize_)); 
            ip.addParameter('output', outputList{1}, @(x) all(ismember(x,outputList)));
            ip.addParameter('projectionAxis3D','Z', @(x) ismember(x,{'Z','X','Y','three'}));
            ip.parse(iChan, varargin{:})
            output = ip.Results.output;
            iFrame = ip.Results.iFrame;
            projAxis3D = ip.Results.projectionAxis3D;
            iZ = ip.Results.iZ;
            varargout = cell(numel(output), 1);
            ZXRatio = obj.owner_.pixelSizeZ_/obj.owner_.pixelSize_;              

            if ischar(output),output={output}; end
            
            for iout = 1:numel(output)
                switch output{iout}             
                    case 'detect3D'
                        s = cached.load(obj.outFilePaths_{1, iChan}, '-useCache', ip.Results.useCache, 'movieInfo');

                        if numel(ip.Results.iFrame)>1
                            v1 = s.movieInfo;
                        else
                            v1 = s.movieInfo(iFrame);
                        end
                        if ~isempty(v1.xCoord) && ~isempty(iZ)
                            % Only show Detections in Z. 
                            zThick = 1;
                            tt = table(v1.xCoord(:,1), v1.yCoord(:,1), v1.zCoord(:,1), 'VariableNames', {'xCoord','yCoord','zCoord'});
                            valid_states = ((tt.zCoord/ZXRatio)>=(iZ-zThick) & (tt.zCoord/ZXRatio)<=(iZ+zThick));
                            dataOut = tt{valid_states, :};

                            if isempty(dataOut) || numel(dataOut) <1 || ~any(valid_states)
                                dataOut = [];
                            end
                        else
                            dataOut = [];
                        end
                        dataOutz = obj.convertProjection3D(dataOut, projAxis3D, ZXRatio);
                        varargout{iout} = dataOutz;

                    case 'detect3Dall'
                        s = cached.load(obj.outFilePaths_{1, iChan}, '-useCache', ip.Results.useCache, 'movieInfo');

                        if numel(ip.Results.iFrame)>1
                            v1 = s.movieInfo;
                        else
                            v1 = s.movieInfo(iFrame);
                        end
                        if ~isempty(v1.xCoord) && ~isempty(iZ)
                            % Only show Detections in Z. 
%                             zThick = 1;
                            tt = table(v1.xCoord(:,1), v1.yCoord(:,1), v1.zCoord(:,1), 'VariableNames', {'xCoord','yCoord','zCoord'});
                            valid_states = ((tt.zCoord/ZXRatio)>=1 & (tt.zCoord/ZXRatio)<=obj.owner_.zSize_);
                            dataOut = tt{:, :};

                            if isempty(dataOut) || numel(dataOut) <1 || ~any(valid_states)
                                dataOut = [];
                            end
                        else
                            dataOut = [];
                        end
                        dataOutz = obj.convertProjection3D(dataOut, projAxis3D, ZXRatio);
                        varargout{iout} = dataOutz;            
                    case 'movieInfo'
                        varargout{iout} = obj.loadChannelOutput@DetectionProcess(iChan, varargin{:});
                    case 'detectionsLabRef'
                        varargout{iout} = load(obj.outFilePaths_{2, iChan}, 'detectionLabRef');
                    otherwise
                        error('Incorrect Output Var type');
                end
            end 
        end
        
        function output = getDrawableOutput(obj)
            output = getDrawableOutput@DetectionProcess(obj);
            output(1).name='Detected Objects by zSlice';
            output(1).var = 'detect3D';
            output(1).formatData=@DetectionProcess.formatOutput3D;
            
            output(2) = getDrawableOutput@DetectionProcess(obj);
            output(2).name='Detected Objects';
            output(2).var = 'detect3Dall';
            output(2).formatData=@DetectionProcess.formatOutput3D;
        end
        

        %% TODO 
        %% draw Amira? function

        
        %% TODO WIP -- TO FINISH (using Philippe'method)
        function scales = getEstSigmaPSF3D(obj, channel, varargin)

            volList=[];
            for i=1:5
                volList = [volList double(obj.owner_.getChannel(channel).loadStack(i))];
            end
            scales = getGaussianPSFsigmaFrom3DData(volList, varargin{:});
            disp(['Estimated scales: ' num2str(scales)]);

        end

    end
    methods (Static)
        
        function name = getName()
            name = 'Point source detection in 3D'; % name need to contain 'Point source detection', otherwise won't see it in the list of DetectionProcess on abstractProcessGUI.
        end
        
        function h = GUI()
            h = @pointSourceDetectionProcessGUI3D;
        end
        
        function funParams = getDefaultParams(owner, varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir = ip.Results.outputDir;
            
            nChan = numel(owner.channels_);
            
            % Set default parameters
            funParams.ChannelIndex = 1:nChan;
            funParams.InputImageProcessIndex = 0; % ?? (can we add some way to check what is availble.)
            funParams.MaskChannelIndex = []; %1:numel(owner.channels_);
            funParams.MaskProcessIndex = [];            
            funParams.OutputDirectory = [outputDir  filesep 'detect3D'];
            funParams.frameRange=[1 owner.nFrames_];
            
%             types = PointSourceDetectionProcess3D.getDetectionTypeOptions;
            funParams.algorithmType = {'pointSourceAutoSigmaFit'};

            funParams.alpha=.05;
            funParams.Mode = {'xyzAc'};
            funParams.FitMixtures = false;
            funParams.MaxMixtures = 5;
            funParams.RemoveRedundant = true;
            funParams.RedundancyRadius = .25;
            funParams.RefineMaskLoG = false;
            funParams.RefineMaskValid = false;
            funParams.ClearMaskBorder = true;
            
            % DetectComets3D & watershed params
            funParams.waterThresh = 120;
            funParams.waterStep = 10;
            funParams.lowFreq = 3;
            funParams.highFreq = 1;
            funParams.showAll = 0;
            funParams.isoCoord = true;
            
            % sigma estimation            
            nChan = numel(owner.channels_);
            funParams.filterSigma = 1.2*ones(1, nChan); %Minimum numerically stable sigma is ~1.2 pixels.
            hasPSFSigma = arrayfun(@(x) ~isempty(x.psfSigma_), owner.channels_);
            
            funParams.filterSigma(hasPSFSigma) = [owner.channels_(hasPSFSigma).psfSigma_];            
            funParams.filterSigma(funParams.filterSigma<1.2) = 1.2; %Make sure default isn't set to too small.
            funParams.filterSigma = repmat(funParams.filterSigma,[2 1]); %TEMP - use z-PSF estimate as well!!
            
            % For the GUI
            funParams.filterSigmaXY = funParams.filterSigma(1,1);
            funParams.filterSigmaZ = funParams.filterSigma(2,1);
            funParams.ConfRadius = arrayfun(@(x)(2*x), funParams.filterSigma);
            funParams.WindowSize = arrayfun(@(x)(ceil(4*x)), funParams.filterSigma);                       

            %list of parameters which can be specified at a per-channel
            %level. If specified as scalar these will  be replicated
            funParams.PerChannelParams = {'alpha','Mode','FitMixtures','MaxMixtures','RedundancyRadius',...
                'ConfRadius','WindowSize','RefineMaskLoG','filterSigma','InputImageProcessIndex',...
                'algorithmType', 'waterThresh', 'waterStep', 'lowFreq', 'highFreq'};
            funParams = prepPerChannelParams(funParams, nChan);
        end
        
        function validTypes =  getValidAlgorithmTypes()
            validTypes = {'watershedApplegateAuto', ...
                          'watershedApplegate',...
                          'bandPassWatershed',...
                          'watershedMatlab',...
                          'markedWatershed',...
                          'pointSourceLM',...
                          'pointSource',...
                          'pointSourceAutoSigma',...
                          'pointSourceAutoSigmaFit',...
                          'pSAutoSigmaMarkedWatershed',...
                          'pointSourceAutoSigmaMixture',... 
                          'pointSourceAutoSigmaLM',...     
                          'pointSourceAutoSigmaFitSig',... 
                          'pSAutoSigmaWatershed'};
        end
    end    
end