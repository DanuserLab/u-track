classdef Animation < hgsetget & matlab.mixin.Copyable & handle
%% This class encapsulate a 2D dynamic or static animation, it can be:
%% - A MIP
%% - An overlay
%% - Graph
%
% Copyright (C) 2019, Danuser Lab - UTSouthwestern 
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
    methods (Abstract)
        img=loadView(obj,fIdx)
        n=getFrameNb(obj)
    end

    
    methods
        
        function saveVideo(obj,pathToVideoFile,varargin)
            ip=inputParser();
            ip.CaseSensitive = false;
            ip.KeepUnmatched = true;
            ip.addParamValue('frameRate',24, @isnumeric);
            ip.addParamValue('quality',100, @isnumeric);
            ip.parse( varargin{:});
            p=ip.Results;

            mkdirRobust(fileparts(pathToVideoFile));
            video = VideoWriter(pathToVideoFile);
            video.FrameRate =p.frameRate;  % Default 30
            video.Quality = p.quality;    % Default 75
            open(video)
            for fIdx=1:obj.getFrameNb()
                img=obj.loadView(fIdx);
                writeVideo(video,img);
            end
            close(video)
        end

        % function imdisp(obj)
        %     imdisp(obj.loadView(1));
        % end

        function hImAnimation=buildImAnimation(obj,filePathTemplate)
            mkdirRobust(fileparts(filePathTemplate));
            for fIdx=1:obj.getFrameNb();
                img=obj.loadView(fIdx);
                imwrite(img, sprintfPath(filePathTemplate, fIdx), 'Compression', 'none');
            end
            hImAnimation=ImAnimation(filePathTemplate,obj.getFrameNb())
        end

    end
end