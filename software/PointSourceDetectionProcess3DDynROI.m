classdef PointSourceDetectionProcess3DDynROI < PointSourceDetectionProcess3D
    % Process Class for detection in 3D for dynamic ROI
    % detectMoviePointSources3D.m is the wrapper function (same as its
    % superclass, so not called here.)
    % PointSourceDetectionProcess3DDynROI is part of New Utrack 3D package
    %
    % This process class is a subclass of PointSourceDetectionProcess3D class.
    % So different getName and funParams.OutputDirectory can be put for this process.
    %
    % Qiongjing (Jenny) Zou, July 2019
%
% Copyright (C) 2025, Danuser Lab - UTSouthwestern 
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

        function obj = PointSourceDetectionProcess3DDynROI(owner, varargin)
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
                super_args{2} = outputDir;
                if isempty(funParams)
                    funParams = PointSourceDetectionProcess3DDynROI.getDefaultParams(owner,outputDir);
                end
                super_args{3} = funParams;
                super_args{4} = PointSourceDetectionProcess3DDynROI.getName;
            end
            obj = obj@PointSourceDetectionProcess3D(super_args{:});
            obj.is3Dcompatible_ = true;
        end

    end
    
    methods (Static)
        
        function name = getName()
            name = 'Point source detection for Dynamic ROI in 3D'; % name need to contain 'Point source detection', otherwise won't see it in the list of DetectionProcess on abstractProcessGUI.
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
            funParams.OutputDirectory = [outputDir  filesep 'detect3DInROI'];
            funParams.frameRange=[1 owner.nFrames_];
            
%             types = PointSourceDetectionProcess3DDynROI.getDetectionTypeOptions;
            % funParams.algorithmType = {'pointSourceAutoSigmaFit'};
            funParams.algorithmType = {'multiscaleDetectionDebug'}; % Change the default for NewUtrack3DPackage, on 2020-11-25

            % funParams.alpha=.05;
            funParams.alpha=.001; % Change the default for NewUtrack3DPackage, on 2020-11-25
            funParams.Mode = {'xyzAc'};
            funParams.FitMixtures = false;
            funParams.MaxMixtures = 5;
            funParams.RemoveRedundant = true;
            funParams.RedundancyRadius = .25;
            funParams.RefineMaskLoG = false;
            funParams.RefineMaskValid = false;
            funParams.ClearMaskBorder = true;
            funParams.processBuildDynROI=[]; % make first available DynROIProc selected&set on the GUI, even default is empty in the process class. edited on 2021-01-04
            funParams.buildDynROIProcessChannel=1;
            funParams.saveMaskFilePattern=[];
            funParams.samplePos=[];


            %% multiscale detector
            funParams.debug=false;
            % funParams.scales=[2:0.5:4];
            funParams.scales=[1.25:0.5:2.25]; % Change the default for NewUtrack3DPackage, on 2020-11-25
            funParams.version='';
            funParams.verbosity=1; % New param added 2022-05, and added again 2022-07
            


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
        
    end    
end