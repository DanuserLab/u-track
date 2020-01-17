classdef CometPostTrackingProcess < PostTrackingProcess
    % A concrete class for classifying comet tracks
    %
    % Sebastien Besson, March 2012
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
        function obj = CometPostTrackingProcess(owner, varargin)
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
                super_args{2} = CometPostTrackingProcess.getName;
                super_args{3} = @postProcessMovieComets;
                if isempty(funParams)  % Default funParams
                    funParams = CometPostTrackingProcess.getDefaultParams(owner,outputDir);
                end
                super_args{4} = funParams;
            end
            
            obj = obj@PostTrackingProcess(super_args{:});
        end
        function varargout = loadChannelOutput(obj,iChan,varargin)
            
            % Input check
            outputList = {'projData','tracks'};
            ip =inputParser;
            ip.addRequired('iChan',@(x) isscalar(x) && obj.checkChanNum(x));
            ip.addOptional('iFrame',1:obj.owner_.nFrames_,@(x) all(obj.checkFrameNum(x)));
            ip.addParameter('useCache',true,@islogical);
            ip.addParamValue('output',outputList,@(x) all(ismember(x,outputList)));
            ip.parse(iChan,varargin{:})
            iFrame = ip.Results.iFrame;
            output = ip.Results.output;
            if ischar(output),output={output}; end
            
            % Data loading
            s = cached.load(obj.outFilePaths_{1,iChan}, '-useCache', ip.Results.useCache,'projData');
            for j=1:numel(output)
                if isequal(output{j},'projData')
                    varargout{1}=s.(output{j});
                else
                    
                    trackData=s.projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix;
                    fullIdx=trackData(:,1);
                    trackType=trackData(:,5);
                    sF=trackData(:,2);
                    [xMat,yMat]=plusTipGetSubtrackCoords(s.projData,[]);
                    
                    correspFullIdx=fullIdx(any(~isnan(xMat(:,iFrame)),2));
                    if ~isempty(correspFullIdx) && max(iFrame)>1
                        subtracks2keep=find(ismember(fullIdx,correspFullIdx));
                        varargout{1}.x=xMat(subtracks2keep,1:max(iFrame));
                        varargout{1}.y=yMat(subtracks2keep,1:max(iFrame));
                        varargout{1}.fullIdx=fullIdx(subtracks2keep);
                        varargout{1}.trackType=trackType(subtracks2keep);
                        varargout{1}.sF=sF(subtracks2keep);
                    else
                        varargout{1}.x=[];
                        varargout{1}.y=[];
                    end
                end
            end
        end
    end
    methods (Static)
        function output = getDrawableOutput()
            output(1).name='Classified tracks';
            output(1).var='tracks';
            output(1).formatData=[];
            output(1).type='overlay';
            output(1).defaultDisplayMethod=@(x) MTTracksDisplay();
        end
        
        function name = getName()
            name = 'Microtubule dynamics classification';
        end
        function h = GUI()
            h = @cometPostTrackingProcessGUI;
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
            funParams.OutputDirectory = [outputDir  filesep 'mtTracks'];
            funParams.makeHist = true;
            funParams.remBegEnd = true;
            funParams.fgapReclassScheme = 1;
            funParams.bgapReclassScheme = 1;
        end
        
        function schemes = getFgapReclassificationSchemes()
            
            schemes{1} = 'Using 2-3 frames before forward gap (local)';
            schemes{2} = 'Using full growth subtrack velocity (local)';
            schemes{3} = 'Unimodal thresholding';
            schemes{4} = 'No reclassification';            
        end
        
        function schemes = getBgapReclassificationSchemes()
            schemes{1} = '95th percentile of forward gap speed distribution';
            schemes{2} = 'Unimodal thresholding with comet latency correction';
            schemes{3} = 'Unimodal thresholding without comet latency correction';
            schemes{4} = 'Using fluctuation radius';
        end
        
    end
end