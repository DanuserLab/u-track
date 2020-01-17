classdef NucleiDetectionProcess < DetectionProcess
    % A concrete class of the detection process
    % Sebastien Besson Oct 2011
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
        function obj = NucleiDetectionProcess(owner, varargin)
            % Input check
            ip = inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.addOptional('funParams',[],@isstruct);
            ip.parse(owner,varargin{:});
            outputDir = ip.Results.outputDir;
            funParams = ip.Results.funParams;
            
            % Define arguments for superclass constructor
            
            % Constructor of the SubResolutionProcess
            super_args{1} = owner;
            super_args{2} = NucleiDetectionProcess.getName;
            super_args{3} = @detectMovieNuclei;            
            if isempty(funParams)  % Default funParams
                funParams=NucleiDetectionProcess.getDefaultParams(owner, outputDir);
            end
            super_args{4} = funParams;
            
            obj = obj@DetectionProcess(super_args{:});
        end
        
   
        function output = getDrawableOutput(obj)
            output = getDrawableOutput@DetectionProcess(obj);
            output(1).name='Nuclei';
        end
        
    end
    methods (Static)
 
        function name = getName()
            name = 'Nuclei Detection';
        end
        function h = GUI()
            h = @nucleiDetectionProcessGUI;
        end
        function filters=getFilters
            filters{1}='sobel';
            filters{2}='canny';
            filters{3}='prewitt';
            
        end
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            % movieParam
            funParams.ChannelIndex = 1:numel(owner.channels_);
            funParams.OutputDirectory = [outputDir  filesep 'detected_nuclei'];
            funParams.ProcessIndex = [];%Default is to use raw images
            funParams.radius = 5;
            funParams.confluent = false;
            funParams.edgeFilter = 'sobel';
            funParams.sigma = 2;
            funParams.useDblLog = true;
            funParams.p = .01;
            funParams.firstFrame = 1;
            funParams.lastFrame = owner.nFrames_;
        end
    end
    
end