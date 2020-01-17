classdef CometDetectionProcess < DetectionProcess
    % A concrete class associated to the comet detection
    %
    % Sebastien Besson, Sep 2011
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
        function obj = CometDetectionProcess(owner, varargin)
            % Constructor of the CometDetectionProcess
            
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
                super_args{2} = CometDetectionProcess.getName;
                super_args{3} = @detectMovieComets;
                if isempty(funParams)  % Default funParams
                    funParams = CometDetectionProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
                
                
            end
            
            obj = obj@DetectionProcess(super_args{:});
        end
        
        function plusTipImport(obj,path)
            obj.outFilePaths_{1}=path;
            obj.procChanged_=false;
            obj.success_=true;
        end
        
        function output = getDrawableOutput(obj)
            output=getDrawableOutput@DetectionProcess(obj);
            output(1).name='Comets';
        end
        
    end
    methods (Static)
        function name = getName()
            name = 'Comet Detection';
        end
        function h = GUI()
            h = @cometDetectionProcessGUI;
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
            funParams.OutputDirectory = [outputDir  filesep 'comets'];
            funParams.MaskProcessIndex = [];
            funParams.MaskChannelIndex =  1 : numel(owner.channels_);
            
            % Detection parameters
            funParams.firstFrame=1;
            funParams.lastFrame=owner.nFrames_;
            funParams.sigma1 = 1;
            funParams.sigma2 = 4;
            funParams.multFactorStepSize=1;
            funParams.multFactorThresh=3;
        end
    end    
end