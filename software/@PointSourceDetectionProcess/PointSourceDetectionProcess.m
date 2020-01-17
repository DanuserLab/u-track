classdef PointSourceDetectionProcess < DetectionProcess
    % A concrete class of a point source detection process
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
        function obj = PointSourceDetectionProcess(owner, outputDir, funParams )
            % Constructor of the SubResolutionProcess
            super_args{1} = owner;
            super_args{2} = PointSourceDetectionProcess.getName;
            super_args{3} = @detectMoviePointSources;
            
            if nargin < 3 || isempty(funParams)  % Default funParams
                if nargin <2, outputDir = owner.outputDirectory_; end
                funParams=PointSourceDetectionProcess.getDefaultParams(owner,outputDir);
            end
            
            super_args{4} = funParams;
            
            obj = obj@DetectionProcess(super_args{:});
        end
        
        function movieInfo = loadChannelOutput(obj,iChan,varargin)
            
            % Input check
            ip =inputParser;
            ip.addRequired('iChan',@(x) ismember(x,1:numel(obj.owner_.channels_)));
            ip.addOptional('iFrame',1:obj.owner_.nFrames_,...
                @(x) ismember(x,1:obj.owner_.nFrames_));
            ip.addParamValue('output', 'movieInfo',@(x) strcmp(x, 'movieInfo'));
            ip.KeepUnmatched = true;
            ip.parse(iChan,varargin{:})
            
            % Data loading
            s = load(obj.outFilePaths_{1,iChan}, 'movieInfo');
            movieInfo = s.movieInfo(ip.Results.iFrame);
        end
        
        function output = getDrawableOutput(obj)
            output = getDrawableOutput@DetectionProcess(obj);
            output(1).name='Point sources';
            output(1).formatData=@PointSourceDetectionProcess.formatOutput;
        end
        
    end
    methods (Static)
        
        function name = getName()
            name = 'Point source detection';
        end
        function h = GUI()
            h = @pointSourceDetectionProcessGUI;
        end
        
        function funParams = getDefaultParams(owner,varargin)
            % Input check
            ip=inputParser;
            ip.addRequired('owner',@(x) isa(x,'MovieData'));
            ip.addOptional('outputDir',owner.outputDirectory_,@ischar);
            ip.parse(owner, varargin{:})
            outputDir=ip.Results.outputDir;
            
            % Set default parameters
            funParams.ChannelIndex=1;
            funParams.MaskChannelIndex = []; %1:numel(owner.channels_);
            funParams.MaskProcessIndex = [];            
            funParams.OutputDirectory = [outputDir  filesep 'point_sources'];
            funParams.alpha=.05;
            funParams.maskRadius=40;
            funParams.Mode = {'xyAc'};
            funParams.FitMixtures = false;
            funParams.MaxMixtures = 5;
            funParams.RedundancyRadius = .25;
            funParams.UseIntersection = true;            
            funParams.PreFilter = true;
            %list of parameters which can be specified at a per-channel
            %level. If specified as scalar these will  be replicated
            funParams.PerChannelParams = {'alpha','Mode','FitMixtures','MaxMixtures','RedundancyRadius','filterSigma','PreFilter','ConfRadius','WindowSize'};
            
            nChan = numel(owner.channels_);
            funParams.filterSigma = 1.2*ones(1,nChan);%Minimum numerically stable sigma is ~1.2 pixels.
            hasPSFSigma = arrayfun(@(x) ~isempty(x.psfSigma_), owner.channels_);
            funParams.filterSigma(hasPSFSigma) = [owner.channels_(hasPSFSigma).psfSigma_];            
            funParams.filterSigma(funParams.filterSigma<1.2) = 1.2;%Make sure default isn't set to too small.
            
            funParams.ConfRadius = arrayfun(@(x)(2*x),funParams.filterSigma);
            funParams.WindowSize = arrayfun(@(x)(ceil(4*x)),funParams.filterSigma);
            
            funParams = prepPerChannelParams(funParams,nChan);

            
        end
        
        function positions = formatOutput(pstruct)
            positions = formatOutput@DetectionProcess(pstruct);
            %positions = positions(pstruct.isPSF, :);
        end
    end    
end