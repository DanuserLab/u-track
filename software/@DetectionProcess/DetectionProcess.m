classdef DetectionProcess < ImageAnalysisProcess
    % An abstract class for all detection processes with an output
    % structure compatible with the tracker
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
    
    % Chuangang Ren, 11/2010
    % Sebastien Besson (last modified May 2012)
    % Mark Kittisopikul, Nov 2014, Added channelOutput cache
    % Andrew R. Jamieson, mar 2017, adding 3D support.
    
    methods(Access = public)
        
        function obj = DetectionProcess(owner, name, funName, funParams )
            
            if nargin == 0
                super_args = {};
            else
                super_args{1} = owner;
                super_args{2} = name;
            end
            
            obj = obj@ImageAnalysisProcess(super_args{:});
            
            if nargin > 2
                obj.funName_ = funName;
            end
            if nargin > 3
                obj.funParams_ = funParams;
            end
        end

        % function overlayCell=displayBSBSBS(obj,iChan)

        %     %% Retrieve projection
        %     displayProjsProcess = obj.getOwner().searchProcessName('RenderDynROI');
        %     displayProjs = cell(1,numel(displayProjsProcess));
            
        %     for pIdx = 1:numel(displayProjsProcess)
        %         displayProjs{pIdx}=displayProjsProcess(pIdx).loadChannelOutput();
        %     end
            
        %     displayProjs=[displayProjs{:}];
        %     overlayCell=cell(1,numel(displayProjs));



        %     for rIdx=1:numel(overlayCell)
        %         overlay = ProjectDynROIRendering(displayProjs{rIdx},['ROI-detection']);
        %         overlay.ZRight = false;
        %         overlayCell{rIdx} = overlay;
        %     end

        %     oDetections = Detections(obj.loadChannelOutput(iChan));
            
        %     for rIdx=1:numel(overlayCell)
        %         myColormap=256*prism(256);
                
        %         overlayProjDetectionMovie(displayProjs{rIdx},'detections',oDetections, 'process',overlayCell{rIdx}, ... 
        %             'cumulative',false,'colormap',myColormap,'radius',5,'detectionBorderDisplay',1.5); 
        %     end


        %     for rIdx=1:numel(displayProjs)
        %         overlayCell{rIdx}.cachedOrtho.imdisp(); %% project data MIPs + detections
        %         overlayCell{rIdx}.draw(1);
        %         drawnow;
        %     end
        % end

        function status = checkChannelOutput(obj,iChan)
            
            %Checks if the selected channels have valid output files
            nChan = numel(obj.owner_.channels_);
            if nargin < 2 || isempty(iChan), iChan = 1:nChan; end
            
            status=  ismember(iChan,1:nChan) & ....
                arrayfun(@(x) exist(obj.outFilePaths_{1,x},'file'),iChan);
        end

        function h=draw(obj,iChan,varargin)
            h = obj.draw@ImageAnalysisProcess(iChan,varargin{:},'useCache',true);
        end
        
        function varargout = loadChannelOutput(obj,iChan,varargin)
            
            % Input check
            outputList = {'movieInfo'};
            ip =inputParser;
            ip.addRequired('iChan',@(x) isscalar(x) && obj.checkChanNum(x));
            ip.addOptional('iFrame',1:obj.owner_.nFrames_,@(x) all(obj.checkFrameNum(x)));
            ip.addParameter('useCache',false,@islogical);
            ip.addParameter('output',outputList,@(x) all(ismember(x,outputList)));
            ip.parse(iChan,varargin{:})
            iFrame = ip.Results.iFrame;
            output = ip.Results.output;
            if ischar(output),output={output}; end
            
            % Data loading
            % load outFilePaths_{1,iChan}
            %
            s = cached.load(obj.outFilePaths_{1,iChan}, '-useCache', ip.Results.useCache, output{:});
           
            if numel(ip.Results.iFrame)>1,
                varargout{1}=s.(output{1});
            else
                varargout{1}=s.(output{1})(iFrame);
            end
        end
        function output = getDrawableOutput(obj)
            colors = parula(numel(obj.owner_.channels_));
            output(1).name='Objects';
            output(1).var='movieInfo';
            output(1).formatData=@DetectionProcess.formatOutput;
            output(1).type='overlay';
            output(1).defaultDisplayMethod=@(x) LineDisplay('Marker','o',...
                'LineStyle','none','Color', colors(x,:));
        end  

        function drawImaris(obj,iceConn)
            
            dataSet = iceConn.mImarisApplication.GetDataSet;
            
            nChan = numel(obj.owner_.channels_);
            
            %Create data container
            spotFolder = iceConn.mImarisApplication.GetFactory.CreateDataContainer;
            spotFolder.SetName(obj.name_);
            iceConn.mImarisApplication.GetSurpassScene.AddChild(spotFolder,-1);
            for iChan = 1:nChan
                
                
                if obj.checkChannelOutput(iChan)
                    %TEMP - doesn't support timelapse yet!!
                    pts = obj.loadChannelOutput(iChan);
                    if ~isempty(pts)
                        chanCol = iceConn.mapRgbaScalarToVector(dataSet.GetChannelColorRGBA(iChan-1));
                        nPts = numel(pts.x);
                        ptXYZ = [pts.y' pts.x' pts.z'];%Swap xy for imaris display
                        ptXYZ(:,1:2) =ptXYZ(:,1:2) * obj.owner_.pixelSize_ / 1e3;
                        ptXYZ(:,3) =ptXYZ(:,3) * obj.owner_.pixelSizeZ_ / 1e3;
                        ptRad = pts.s(1,:)' * obj.owner_.pixelSize_ / 1e3;
                        ptObj = iceConn.createAndSetSpots(ptXYZ,zeros(nPts,1),...
                            ptRad,['Spots ' char(dataSet.GetChannelName(iChan-1))],chanCol,spotFolder);                    

                        if ismethod(ptObj,'SetRadiiXYZ')
                            %Make sure we have a version of imaris which supports anisotropic points                        
                            ptRad = zeros(nPts,3);
                            ptRad(:,1:2) = repmat(pts.s(1,:)' * obj.owner_.pixelSize_ / 1e3,[1 2]);
                            ptRad(:,3) = pts.s(2,:) * obj.owner_.pixelSizeZ_ / 1e3;
                            ptObj.SetRadiiXYZ(ptRad);                                            
                        end
                    end
                end
            end            
        end    
          function y = convertProjection3D(obj, x, zAxis, ZXRatio)
            switch zAxis
                case 'Y' %(ZX)
%                     y = horzcat(x(:,3)*ZXRatio,x(:,1),x(:,2));
                    y = horzcat(x(:,3), x(:,1), x(:,2)*ZXRatio);
                case 'X' % (ZY)
%                     y = horzcat(x(:,3)*ZXRatio,x(:,2),x(:,1));
                    y = horzcat(x(:,3), x(:,2), x(:,1)*ZXRatio);
                case 'Z'
                    y = x;
                case 'three'
                    yT = x;
                    zx = horzcat(x(:,1),x(:,3)+4+obj.owner_.imSize_(1), x(:,2)*ZXRatio);
                    zy = horzcat(x(:,3)+obj.owner_.imSize_(2)+4, x(:,2), x(:,1)*ZXRatio);
                    y = vertcat(yT, zx, zy);
            end
        end
    end
    methods(Static)

        function name = getName()
            name = 'Detection';
        end

        function h = GUI()
            h = @abstractProcessGUI;
        end

        function procClasses = getConcreteClasses(varargin)

            procClasses = ...
                {@SubResolutionProcess;
                 @CometDetectionProcess;
                 @AnisoGaussianDetectionProcess;
                 @NucleiDetectionProcess;
                 @PointSourceDetectionProcess;
                 @ExternalDetectionProcess; % both 2D and 3D
                 @PointSourceDetectionProcess3D;
                };

            % If input, check if 2D or 3D movie(s).
            ip =inputParser;
            ip.addOptional('MO', [], @(x) isa(x,'MovieData') || isa(x,'MovieList'));
            ip.parse(varargin{:});
            MO = ip.Results.MO;
            
            if ~isempty(MO)
                if isa(MO,'MovieList')
                    MD = MO.getMovie(1);
                elseif length(MO) > 1
                    MD = MO(1);
                else
                    MD = MO;
                end                
            end

            if isempty(MD)
               warning('MovieData properties not specified (2D vs. 3D)');
               disp('Displaying both 2D and 3D Detection processes');
            elseif MD.is3D
                disp('Detected 3D movie');
                disp('Displaying 3D Detection processes only');
                procClasses(1:5) = [];
            elseif ~MD.is3D
                disp('Detected 2D movie');
                disp('Displaying 2D Detection processes only');
                procClasses(7) = [];
            end
            procClasses = cellfun(@func2str, procClasses, 'Unif', 0);
        end
        
        function y = formatOutput(x)
            % Format output in xy coordinate system (default-backwards compatible)
            y = DetectionProcess.formatOutput2D(x);
        end

        function y = formatOutput2D(x)
            % Format output in xy coordinate system
            if isempty(x.xCoord)
                y = NaN(1,2);
            else
                y = horzcat(x.xCoord(:,1),x.yCoord(:,1));
            end
        end


        function y = formatOutput3D(x)
            if isempty(x)
                y = NaN(1,3);
            else
                y = x;
            end
        end

    end
end
