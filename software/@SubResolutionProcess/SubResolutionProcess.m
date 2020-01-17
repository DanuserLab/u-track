classdef SubResolutionProcess < DetectionProcess
    % A concrete class for detecting objects using Gaussian mixture-model fitting
    % Chuangang Ren 11/2010
    % Sebastien Besson (last modified Dec 2011)
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
        function obj = SubResolutionProcess(owner, varargin)
            % Input check
            ip = inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.addOptional('funParams',[],@isstruct);
            ip.parse(owner,varargin{:});
            outputDir = ip.Results.outputDir;
            funParams = ip.Results.funParams;
            
            
            % Constructor of the SubResolutionProcess
            
            super_args{1} = owner;
            super_args{2} = SubResolutionProcess.getName;
            super_args{3} = @detectMovieSubResFeatures;
            if isempty(funParams)  % Default funParams
                funParams = SubResolutionProcess.getDefaultParams(owner,outputDir);
            end
            super_args{4} = funParams;
            
            obj = obj@DetectionProcess(super_args{:});            
        end
      
 
        function output = getDrawableOutput(obj)
            % Rename default detection output
            output = getDrawableOutput@DetectionProcess(obj);
            output(1).name='Sub-resolution objects';
        end
        
    end
    methods (Static)
        
        function name = getName()
            name = 'Gaussian Mixture-Model Fitting';
        end
        function h = GUI()
            h = @subResolutionProcessGUI;
        end
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            % moviePara  
            funParams.ChannelIndex =1:numel(owner.channels_);
            funParams.OutputDirectory = [outputDir  filesep 'GaussianMixtureModels'];
            funParams.firstImageNum = 1;
            funParams.lastImageNum = owner.nFrames_;
            
            % detectionParam
            if ~isempty(owner.channels_(1).psfSigma_)
                funParams.detectionParam.psfSigma = mean(horzcat(owner.channels_.psfSigma_));%%%% (1) for testing
            else
                funParams.detectionParam.psfSigma=[];
            end
            if ~isempty(owner.camBitdepth_)
                funParams.detectionParam.bitDepth = owner.camBitdepth_;
            else
                funParams.detectionParam.bitDepth = [];
            end
            funParams.detectionParam.alphaLocMax = .05;
            funParams.detectionParam.integWindow = 0;
            funParams.detectionParam.doMMF = 0;
            funParams.detectionParam.testAlpha = struct('alphaR', .05,'alphaA', .05, 'alphaD', .05,'alphaF',0);
            funParams.detectionParam.numSigmaIter = 0;
            funParams.detectionParam.visual = 0;
            funParams.detectionParam.background = [];
            
        end

    end
    
end